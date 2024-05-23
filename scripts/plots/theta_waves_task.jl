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
using SpatiotemporalMotifs
@preamble
set_theme!(foresight(:physics))
Random.seed!(32)

stimulus = r"Natural_Images"
vars = [:k]
datafile = datadir("theta_waves_task.jld2")

session_table = load(datadir("session_table.jld2"), "session_table")
oursessions = session_table.ecephys_session_id

path = datadir("calculations")
Q = calcquality(path)[structure = At(structures)]
quality = mean(Q[stimulus = At(stimulus)])

config = @strdict stimulus vars
data, file = produce_or_load(produce_out(Q), config, datadir(); filename = savepath,
                             prefix = "out")
out = data["out"]

begin # * Supplemental material: average phase velocity maps in each region
    f = Figure(size = (720, 1440))

    for i in eachindex(uni)
        k = uni[i][:k][:, 2:end, :]
        k = uconvert.(u"mm^-1", k)

        # * Hit
        ax = Axis(f[i, 1], yreversed = true)
        structure = metadata(k)[:structure]
        ax.title = structure * ": hit"
        m = median(k[:, :, lookup(k, :trial) .== true], dims = :trial)
        m = dropdims(m, dims = :trial)
        colorrange = maximum(abs.(ustripall(m))) * [-1, 1]
        ints = uni[i][:layerints]
        p = plotlayermap!(ax, m, ints; arrows = true, colorrange) |> first
        if i > 5
            ax.xlabel = "Time (s)"
        end

        # * Miss
        ax = Axis(f[i, 2], yreversed = true)
        ax.title = structure * ": miss"
        m = median(k[:, :, lookup(k, :trial) .== false], dims = :trial)
        m = dropdims(m, dims = :trial)
        ints = uni[i][:layerints]
        p = plotlayermap!(ax, m, ints; arrows = true, colorrange) |> first
        c = Colorbar(f[i, 3], p)
        c.label = "k̃ ($(unit(eltype(k))))"
    end
    display(f)
    wsave(plotdir("theta_waves_task", "supplemental_wavenumber.pdf"), f)
end

begin # * Set up main figure
    fig = FourPanel()
    gs = subdivide(fig, 2, 2)
end

begin # * Calculate a global order parameter at each time point
    Og = map(out) do o
        O = map(o) do p
            k = p[:k]
            ret = dropdims(mean(sign.(k), dims = Dim{:depth}), dims = Dim{:depth})
            trials = p[:trials][1:size(k, :changetime), :]
            trialtimes = trials.change_time_with_display_delay
            @assert all(isapprox.(ustripall(lookup(k, :changetime)), trialtimes,
                                  atol = 0.05)) # These won't match exactly, because the data change times have been adjusted for rectification. This is OK.
            set(ret, Dim{:changetime} => Dim{:trial}(trials.hit))
        end
    end

    Og_h = [[o[:, lookup(o, :trial) .== true] for o in O] for O in Og]
    Og_m = [[o[:, lookup(o, :trial) .== false] for o in O] for O in Og]
end

begin # * Plot the mean order parameter across time, and the hit/miss contrast
    function mf(Og)
        map(Og) do O
            o = dropdims.(mean.(O; dims = Dim{:trial}), dims = Dim{:trial})
            sessionids = [metadata(_o)[:sessionid] for _o in O]
            o = [_o[(end - minimum(length.(o)) + 1):end] for _o in o]
            stack(Dim{:sessionid}(sessionids), o)
        end
    end
    Ō = mf(Og)
    Ō_h = mf(Og_h)
    Ō_m = mf(Og_m)

    ax = Axis(gs[1]; xlabel = "Time (s)", ylabel = "Mean order parameter",
              xautolimitmargin = (0, 0), xminorticksvisible = true,
              xminorticks = IntervalsBetween(5), yminorticksvisible = true,
              yminorticks = IntervalsBetween(5))
    hlines!(ax, [0]; color = (:black, 0.5), linestyle = :dash, linewidth = 2)
    vlines!(ax, [0, 0.25]; color = (:black, 0.5), linestyle = :dash, linewidth = 2)
    for (i, O) in reverse(collect(enumerate(Ō)))
        structure = metadata(O)[:structure]
        μ = dropdims(mean(O, dims = Dim{:sessionid}), dims = Dim{:sessionid})
        σ = dropdims(std(O, dims = Dim{:sessionid}), dims = Dim{:sessionid})
        σ = σ ./ sqrt(size(O, 2)) # SEM
        bargs = [times(μ), μ .- σ, μ .+ σ] .|> ustripall .|> collect
        band!(ax, bargs..., color = (structurecolors[i], 0.3))
        lines!(ax, μ |> ustripall, color = (structurecolors[i], 0.7), label = structure)
    end
    l = axislegend(ax, position = :lt, nbanks = 2, framevisible = true, labelsize = 12)
    reverselegend!(l)

    f = Figure()
    ax = Axis(f[1, 1]; xlabel = "Time (s)", ylabel = "Mean order parameter (hit - miss)",
              xautolimitmargin = (0, 0), xminorticksvisible = true,
              xminorticks = IntervalsBetween(5), yminorticksvisible = true,
              yminorticks = IntervalsBetween(5))
    hlines!(ax, [0]; color = (:black, 0.5), linestyle = :dash, linewidth = 2)
    vlines!(ax, [0, 0.25]; color = (:black, 0.5), linestyle = :dash, linewidth = 2)
    for (i, (O_h, O_m)) in enumerate(zip(Ō_h, Ō_m))
        structure = metadata(O_h)[:structure]
        x = O_h .- O_m
        μ = dropdims(mean(x, dims = Dim{:sessionid}), dims = Dim{:sessionid})
        σ = dropdims(std(x, dims = Dim{:sessionid}), dims = Dim{:sessionid})
        σ = σ ./ sqrt(size(O_h, 2)) # SEM
        bargs = [times(μ), μ .- σ, μ .+ σ] .|> ustripall .|> collect
        band!(ax, bargs..., color = (structurecolors[i], 0.3))
        lines!(ax, μ |> ustripall, color = (structurecolors[i], 0.7), label = structure)
    end
    axislegend(ax, position = :rb, framevisible = true, labelsize = 12)
    display(f)
    wsave(plotdir("theta_waves_task", "regional_orderparameter.pdf"), f)
end

# ---------------------------------- Hit/miss prediction --------------------------------- #
begin # * Filter out poor-performing mice and regenerate order parameters
    performance = load_performance(; path)
    newsessions = performance.sessionid[performance.mean_dprime .> 1]
    out = [filter(x -> x.sessionid ∈ newsessions, o) for o in out]

    Og = map(out) do o
        O = map(o) do p
            k = p[:k]
            ret = dropdims(mean(sign.(k), dims = Dim{:depth}), dims = Dim{:depth})
            trials = p[:trials]
            trialtimes = trials.change_time_with_display_delay
            @assert all(isapprox.(ustripall(lookup(k, :changetime)), trialtimes,
                                  atol = 0.05)) # These won't match exactly, because the data change times have been adjusted for rectification. This is OK.
            ret = set(ret, Dim{:changetime} => Dim{:trial})
            set(ret, DimensionalData.format(Dim{:trial}(trials.hit), lookup(ret, :trial)))
        end
    end

    Og_h = [[o[:, lookup(o, :trial) .== true] for o in O] for O in Og]
    Og_m = [[o[:, lookup(o, :trial) .== false] for o in O] for O in Og]
    sessionids = [o[:sessionid] for o in out[1]]
    GC.gc()
end

begin # * LDA input data and downsampling
    H = [getindex.(Og, s) for s in eachindex(Og[1])] # Grouped by subject

    # * Make sure temporal dims same size
    ts = intersect([intersect(lookup.(h, Ti)...) for h in H]...)
    H = [[_h[Ti = At(ts)] for _h in h] for h in H]
    H = stack.([Dim{:structure}(structures)], H, dims = 3)
    H = permutedims.(H, [(1, 3, 2)])

    H = [h[1:10:end, :, :] for h in H]
end

if !isfile(datadir("hyperparameters", "theta_waves_task.jld2")) ||
   !isfile(datafile) # Run calculations; needs to be on a cluster
    if haskey(ENV, "JULIA_DISTRIBUTED")
        procs = USydClusters.Physics.addprocs(32; mem = 15, ncpus = 2,
                                              project = projectdir())
        @everywhere using SpatiotemporalMotifs
    else
        error("Calculations must be run on a cluster, set ENV[\"JULIA_DISTRIBUTED\"] to confirm this.")
    end
    begin # * Get a rough estimate for a good hyperparameter. Currently on pre-offset data. This gives us ~0.5 as a good estimate
        function tuneclassifier(H; r0 = 0.5, repeats = 1, kwargs...)
            r0 = log10(r0)
            objective = r -> -classify_kfold(H[Ti = -0.25u"s" .. 0.25u"s"];
                                             regcoef = exp10(only(r)), repeats, kwargs...)
            o = optimize(objective, [r0], ParticleSwarm(), Optim.Options(; iterations = 10))
            return exp10(only(o.minimizer)), -o.minimum
        end
        hfile = datadir("hyperparameters", "theta_waves_task.jld2")
        hyperr = pmap(tuneclassifier, H)
        tagsave(hfile, @strdict hyperr)
        if isfile(hfile)
            hyperr = load(hfile, "hyperr")
            f = Figure()
            ax = Axis(f[1, 1]; xlabel = "Regularization coefficient",
                      ylabel = "Balanced accuracy",
                      title = "Hyperparameter tuning")
            scatter!(ax, first.(hyperr), last.(hyperr))
            f
            save(datadir("hyperparameters", "theta_waves_task.pdf"), f)
        end
    end

    begin # * Single-subject classifications, returning 5-fold balanced accuracy
        regcoef = 0.5
        folds = 5
        repeats = 10
        bac_pred = pmap(H) do h
            h = h[Ti = -0.25u"s" .. 0.0u"s"]
            bac = classify_kfold(h; regcoef, k = folds, repeats)
        end

        bac_pre = pmap(H) do h
            h = h[Ti = -0.25u"s" .. 0.25u"s"]
            bac = classify_kfold(h; regcoef, k = folds, repeats)
        end

        bac_post = pmap(H) do h
            h = h[Ti = 0.25u"s" .. 0.75u"s"]
            bac = classify_kfold(h; regcoef, k = folds, repeats)
        end

        bac_sur = pmap(H) do h
            h = h[Ti = -0.25u"s" .. 0.25u"s"]
            idxs = randperm(size(h, :trial))
            h = set(h, Dim{:trial} => lookup(h, :trial)[idxs])
            bac = classify_kfold(h; regcoef, k = folds, repeats)
        end

        D = @strdict bac_pred bac_pre bac_post bac_sur regcoef folds repeats
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
    push!(D, "W" => w)
    tagsave(datafile, D)
end

if isfile(datafile)
    D = jldopen(datafile, "r")
    @unpack bac_pre, bac_post, bac_sur, W = D
    close(D)
else
    error("Please produce or download $datafile")
end

begin # * Plot classification performance
    ax = Axis(gs[3]; xlabel = "Balanced accuracy", ylabel = "Density",
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

# -------------------------------- Rection-time prediction ------------------------------- #
begin # * K-fold validation with Huber (outlier-robust) regression
    reaction_times = [o[:trials].lick_latency for o in first(out)]
    regcoef = 1e5
    bac = pmap(eachindex(H)) do i
        regress_kfold(H[i][Ti = -0.0u"s" .. 0.25u"s"], reaction_times[i]; regcoef,
                      k = 5)
    end
    bac_sur = pmap(eachindex(H)) do i
        regress_kfold(H[i][Ti = -0.0u"s" .. 0.25u"s"],
                      reaction_times[i][randperm(length(reaction_times[i]))]; regcoef,
                      k = 5)
    end
end
begin # * Plot regression scores
    ax = Axis(gs[4]; xlabel = "Spearman correlation", ylabel = "Density",
              title = "Reaction time prediction")
    ziggurat!(ax, bac_sur; bins = 10, label = "Null", normalizaton = :pdf, color = crimson)
    # ziggurat!(ax, bac_post; bins = 10, label = "Post-stimulus", normalizaton = :pdf,
    #           color = cucumber)
    ziggurat!(ax, bac; bins = 10, label = "Pre-offset", normalizaton = :pdf,
              color = cornflowerblue)
    axislegend(ax, position = :lt, framevisible = true, labelsize = 12)
    tightlimits!(ax)
end

begin # * Save
    addlabels!(fig)
    wsave(plotdir("theta_waves_task", "theta_waves_task.pdf"), fig)
end
