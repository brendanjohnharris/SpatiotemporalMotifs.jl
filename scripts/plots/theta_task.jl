#! /bin/bash
#=
exec julia -t auto "${BASH_SOURCE[0]}" "$@"
=#
# ? Expected execution time: 15 mins
using DrWatson
@quickactivate "SpatiotemporalMotifs"

import TimeseriesTools: freqs
using FileIO
using Unitful
import DimensionalData: metadata
using MultivariateStats
using USydClusters
using SpatiotemporalMotifs
structurecolors, structures, commondepths, parselayernum,
load_calculations, unify_calculations,
plotlayerints!, plotlayermap!,
@preamble
set_theme!(foresight(:physics))
Random.seed!(32)

stimulus = r"Natural_Images"
vars = [:x]

session_table = load(datadir("session_table.jld2"), "session_table")
oursessions = session_table.ecephys_session_id

path = datadir("calculations")
Q = calcquality(path)[structure = At(structures)]
quality = mean(Q[stimulus = At(stimulus)])

config = @strdict stimulus vars
data, file = produce_or_load(produce_out(Q), config, datadir(); filename = savepath,
                             prefix = "out")
out = data["out"]

begin # * Set up main figure
    fig = FourPanel()
    gs = subdivide(fig, 2, 2)
end

# ---------------------------------- Hit/miss prediction --------------------------------- #
begin # * Filter out poor-performing mice and regenerate order parameters
    performance = load_performance(; path)
    newsessions = performance.sessionid[performance.mean_dprime .> 1]
    out = [filter(x -> x.sessionid ∈ newsessions, o) for o in out]

    Og = map(out) do o
        O = map(o) do p
            x = p[:x]
            trials = p[:trials]
            trialtimes = trials.change_time_with_display_delay
            @assert all(isapprox.(ustripall(lookup(x, :changetime)), trialtimes,
                                  atol = 0.05)) # These won't match exactly, because the data change times have been adjusted for rectification. This is OK.
            x = set(x, Dim{:changetime} => Dim{:trial})
            set(x, DimensionalData.format(Dim{:trial}(trials.hit), lookup(x, :trial)))
        end
    end

    Og_h = [[o[:, :, lookup(o, :trial) .== true] for o in O] for O in Og]
    Og_m = [[o[:, :, lookup(o, :trial) .== false] for o in O] for O in Og]
    sessionids = [o[:sessionid] for o in out[1]]
    GC.gc()
end

function unflatten!(T, x) # T is a template to write into
    n = 0
    for i in eachindex(T)
        T[i] .= x[(n + 1):(n + length(T[i]))]
        n += length(T[i])
    end
end

begin # * LDA input data and downsampling
    H = [getindex.(Og, s) for s in eachindex(Og[1])] # Grouped by subject

    # * Make sure temporal dims same size
    ts = intersect([intersect(lookup.(h, Ti)...) for h in H]...)
    H = [[_h[Ti = At(ts)] for _h in h] for h in H]
    H = map(H) do h
        map(h) do _h
            _h[1:10:end, 1:3:end, :] # * Downsample
        end
    end
    H_pred = map(H) do h # How well can we predict the mouses response from data prior to the stimulus?
        map(h) do _h
            _h[Ti = -0.25u"s" .. 0.0u"s"]
        end
    end
    H_pre = map(H) do h # Or using data prior to the stimulus offset?
        map(h) do _h
            _h[Ti = -0.25u"s" .. 0.25u"s"]
        end
    end
    H_post = map(H) do h # Or using data following stimulus offset?
        map(h) do _h
            _h[Ti = 0.25u"s" .. 0.75u"s"]
        end
    end
    T_pre = map(H_pre) do h # Container for weights later on
        map(h) do _h
            similar(_h[:, :, 1])
        end
    end
    T_post = map(H_post) do h # Container for weights later on
        map(h) do _h
            similar(_h[:, :, 1])
        end
    end
    function toflat(H)
        H = map(H) do h
            map(h) do _h
                x = stack(map(Iterators.flatten, eachslice(_h, dims = 3)), dims = 1)
                x = DimArray(x', (DimensionalData.AnonDim(axes(x, 2)), dims(_h, :trial)))
            end
        end
        H = map(H) do h
            cat(h..., dims = 1)
        end
        return H
    end
    H_pred = toflat(H_pred)
    H_pre = toflat(H_pre)
    H_post = toflat(H_post)
end

begin
    procs = USydClusters.Physics.addprocs(32; mem = 15, ncpus = 2,
                                          project = projectdir())
    @everywhere using SpatiotemporalMotifs
end

begin # * Get a rough estimate for a good hyperparameter. Currently on pre-offset data. This gives us ~0.5 as a good estimate
    @everywhere begin
        using Optim
        function tuneclassifier(H; r0 = 0.5, repeats = 1, kwargs...)
            @info "Tuning classifier"
            r0 = log10(r0)
            objective = r -> -classify_kfold(H;
                                             regcoef = exp10(only(r)), repeats, kwargs...)
            o = optimize(objective, [r0], ParticleSwarm(), Optim.Options(; iterations = 10))
            return exp10(only(o.minimizer)), -o.minimum
        end
    end
    hfile = datadir("hyperparameters", "theta_task.jld2")
    if haskey(ENV, "JULIA_DISTRIBUTED")
        hyperr = pmap(tuneclassifier, H_pre) # Takes about .. minutes over 96 cores
        tagsave(hfile, @strdict hyperr)
    end
    if isfile(hfile)
        hyperr = load(hfile, "hyperr")
        f = Figure()
        ax = Axis(f[1, 1]; xlabel = "Regularization coefficient",
                  ylabel = "Balanced accuracy",
                  title = "Hyperparameter tuning")
        scatter!(ax, first.(hyperr), last.(hyperr))
        f
        save(datadir("hyperparameters", "theta_task.pdf"), f)
    end
end

begin # * Single-subject classifications, returning 5-fold balanced accuracy. Takes around 30 minutes across 30 cores.
    regcoef = 1
    folds = 5
    repeats = 10
    bac_pred = pmap(H_pred; distributed = true) do h
        @info "Commencing calculations"
        bac = classify_kfold(h; regcoef, k = folds)
    end
    bac_pre = pmap(H_pre; distributed = true) do h
        bac = classify_kfold(h; regcoef, k = folds)
    end
    bac_pre = pmap(H_pre; distributed = true) do h
        bac = classify_kfold(h; regcoef, k = folds)
    end

    bac_post = p(H_post; distributed = true) do h
        bac = classify_kfold(h; regcoef, k = folds)
    end

    bac_sur = p(H_post; distributed = true) do h
        idxs = randperm(size(h, :trial))
        h = set(h, Dim{:trial} => lookup(h, :trial)[idxs])
        bac = classify_kfold(h; regcoef, k = folds)
    end
end

begin # * Plot classification performance
    f = Figure()
    ax = Axis(f[1, 1]; xlabel = "Balanced accuracy", ylabel = "Density",
              title = "Hit/miss classification")
    ziggurat!(ax, bac_sur; bins = 10, label = "Null", normalizaton = :pdf, color = crimson)
    ziggurat!(ax, bac_post; bins = 10, label = "Post-offset", normalizaton = :pdf,
              color = cucumber)
    ziggurat!(ax, bac_pre; bins = 10, label = "Pre-offset", normalizaton = :pdf,
              color = cornflowerblue)
    axislegend(ax, position = :lt, framevisible = true, labelsize = 12)
    tightlimits!(ax)
    f
end

begin # * Map of region-wise weightings
    W = map(H) do h
        h = h[Ti = -0.25u"s" .. 0.75u"s"] # !!!
        N, M = classifier(h; regcoef = 0.5) # !!!
        W = projection(M)
        W = reshape(W, size(h)[1:2])
        return DimArray(W, dims(h)[1:2])
    end
    W = stack(Dim{:sessionid}(sessionids), W, dims = 3)
    W = W ./ maximum(abs.(W))
end

begin # * Plot region-wise weightings
    ax = Axis(gs[2]; xlabel = "Time (s)", ylabel = "Normalized LDA weight",
              title = "Regional classification weights")
    vlines!(ax, [0, 0.25], color = (:black, 0.2), linestyle = :dash)
    hlines!(ax, [0], color = (:black, 0.5), linewidth = 2)
    for structure in reverse(collect(lookup(W, :structure)))
        ws = W[structure = At(structure)]
        ws = upsample(ws, 5, 1)
        ts = ustripall(lookup(ws, Ti))
        μ = mean(ws, dims = :sessionid) |> vec
        σ = std(ws, dims = :sessionid) ./ sqrt(size(ws, :sessionid)) |> ustripall |> vec
        band!(ax, ts, μ - σ, μ + σ; color = (structurecolormap[structure], 0.3))
        lines!(ax, ts, μ; color = structurecolormap[structure], label = structure)
    end
    tightlimits!(ax)
    l = axislegend(ax, position = :lt, nbanks = 2, framevisible = true, labelsize = 12)
    reverselegend!(l)
end

# # -------------------------------- Rection-time prediction ------------------------------- #
# begin # * K-fold validation with Huber (outlier-robust) regression
#     reaction_times = [o[:trials].lick_latency for o in first(out)]
#     regcoef = 1e5
#     bac = pmap(eachindex(H)) do i
#         regress_kfold(H[i][Ti = -0.0u"s" .. 0.25u"s"], reaction_times[i]; regcoef,
#                       k = 5)
#     end
#     bac_sur = pmap(eachindex(H)) do i
#         regress_kfold(H[i][Ti = -0.0u"s" .. 0.25u"s"],
#                       reaction_times[i][randperm(length(reaction_times[i]))]; regcoef,
#                       k = 5)
#     end
# end
# begin # * Plot regression scores
#     ax = Axis(gs[4]; xlabel = "Spearman correlation", ylabel = "Density",
#               title = "Reaction time prediction")
#     ziggurat!(ax, bac_sur; bins = 10, label = "Null", normalizaton = :pdf, color = crimson)
#     # ziggurat!(ax, bac_post; bins = 10, label = "Post-stimulus", normalizaton = :pdf,
#     #           color = cucumber)
#     ziggurat!(ax, bac; bins = 10, label = "Pre-offset", normalizaton = :pdf,
#               color = cornflowerblue)
#     axislegend(ax, position = :lt, framevisible = true, labelsize = 12)
#     tightlimits!(ax)
# end

begin # * Save
    addlabels!(fig)
    wsave(plotdir("theta_waves_task", "theta_waves_task.pdf"), fig)
end
