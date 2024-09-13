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

vars = [:x, :k, :ω]
datafile = datadir("theta_waves_task.jld2")
session_table = load(datadir("posthoc_session_table.jld2"), "session_table")
oursessions = session_table.ecephys_session_id
path = datadir("calculations")

begin # * Set up main figure
    fig = FourPanel()
    gs = subdivide(fig, 2, 2)

    function orderparameter(Og)
        map(Og) do O
            o = dropdims.(nansafe(mean; dims = Trial).(O), dims = Trial)
            sessionids = [metadata(_o)[:sessionid] for _o in O]
            o = [_o[(end - minimum(length.(o)) + 1):end] for _o in o]
            stack(SessionID(sessionids), o)
        end
    end
end

# * Flash stimulus
stimulus = "flash_250ms"
Q = calcquality(path)[Structure = At(structures)]
quality = mean(Q[stimulus = At(stimulus)])
out = load_calculations(Q; stimulus, vars)

begin # * Calculate a global order parameter at each time point
    out = map(out) do o
        filter(o) do _o
            _o[:sessionid] ∈ oursessions
        end
    end
    Og = map(out) do o
        O = map(o) do p
            k = p[:k]
            ω = p[:ω]
            k[ustripall(ω) .< 0] .= NaN * unit(eltype(k))
            ret = dropdims(nansafe(mean; dims = Depth)(sign.(k)), dims = Depth) # Ignore negative frequencies
            set(ret, Dim{:changetime} => Trial)
        end
    end
end

begin # * Plot the mean order parameter across time
    Ō = orderparameter(Og)

    ax = Axis(gs[2]; xlabel = "Time (s)", ylabel = "Mean order parameter",
              title = "Order parameter during flashes",
              xautolimitmargin = (0, 0), xminorticksvisible = true,
              xminorticks = IntervalsBetween(5), yminorticksvisible = true,
              yminorticks = IntervalsBetween(5))
    hlines!(ax, [0]; color = (:black, 0.5), linestyle = :dash, linewidth = 2)
    vlines!(ax, [0, 0.25]; color = (:black, 0.5), linestyle = :dash, linewidth = 2)
    for (i, O) in reverse(collect(enumerate(Ō)))
        structure = metadata(O)[:structure]
        O = O[𝑡(SpatiotemporalMotifs.INTERVAL)]
        μ = dropdims(mean(O, dims = SessionID), dims = SessionID)
        σ = dropdims(std(O, dims = SessionID), dims = SessionID)
        σ = σ ./ 2
        bargs = [times(μ), μ .- σ, μ .+ σ] .|> ustripall .|> collect
        band!(ax, bargs..., color = (structurecolors[i], 0.3), label = structure)
        lines!(ax, times(μ) |> ustripall, μ |> ustripall |> collect,
               color = (structurecolors[i], 0.7),
               label = structure)
    end
    l = axislegend(ax, position = :lt, nbanks = 2, framevisible = true, labelsize = 12,
                   merge = true)
    reverselegend!(l)
end

# * Task stimulus
stimulus = r"Natural_Images"
Q = calcquality(path)[Structure = At(structures)]
quality = mean(Q[stimulus = At(stimulus)])
out = load_calculations(Q; stimulus, vars)

begin # * Calculate a global order parameter at each time point
    out = map(out) do o
        filter(o) do _o
            _o[:sessionid] ∈ oursessions
        end
    end
    Og = map(out) do o
        O = map(o) do p
            k = p[:k]
            ω = p[:ω]
            k[ustripall(ω) .< 0] .= NaN * unit(eltype(k))
            ret = dropdims(nansafe(mean, dims = Depth)(sign.(k)), dims = Depth)
            trials = p[:trials][1:size(k, :changetime), :]
            trialtimes = trials.change_time_with_display_delay
            @assert all(isapprox.(ustripall(lookup(k, :changetime)), trialtimes,
                                  atol = 0.05)) # These won't match exactly, because the data change times have been adjusted for rectification. This is OK.
            set(ret, Dim{:changetime} => Trial(trials.hit))
        end
    end

    Og_h = [[o[:, lookup(o, Trial) .== true] for o in O] for O in Og]
    Og_m = [[o[:, lookup(o, Trial) .== false] for o in O] for O in Og]
end

begin # * Plot the mean order parameter across time
    Ō = orderparameter(Og)
    Ō_h = orderparameter(Og_h)
    Ō_m = orderparameter(Og_m)

    ax = Axis(gs[1]; xlabel = "Time (s)", ylabel = "Mean order parameter",
              title = "Order parameter during task",
              xautolimitmargin = (0, 0), xminorticksvisible = true,
              xminorticks = IntervalsBetween(5), yminorticksvisible = true,
              yminorticks = IntervalsBetween(5))
    hlines!(ax, [0]; color = (:black, 0.5), linestyle = :dash, linewidth = 2)
    vlines!(ax, [0, 0.25]; color = (:black, 0.5), linestyle = :dash, linewidth = 2)
    for (i, O) in reverse(collect(enumerate(Ō)))
        structure = metadata(O)[:structure]
        O = O[𝑡(SpatiotemporalMotifs.INTERVAL)]
        μ = dropdims(mean(O, dims = SessionID), dims = SessionID)
        σ = dropdims(std(O, dims = SessionID), dims = SessionID)
        σ = σ ./ 2
        bargs = [times(μ), μ .- σ, μ .+ σ] .|> ustripall .|> collect
        band!(ax, bargs..., color = (structurecolors[i], 0.3), label = structure)
        lines!(ax, times(μ) |> ustripall, μ |> ustripall |> collect,
               color = (structurecolors[i], 0.7),
               label = structure)
    end
    l = axislegend(ax, position = :lt, nbanks = 2, framevisible = true, labelsize = 12,
                   merge = true)
    reverselegend!(l)

    # f = Figure()
    # ax = Axis(f[1, 1]; xlabel = "Time (s)", ylabel = "Mean order parameter (hit - miss)",
    #           xautolimitmargin = (0, 0), xminorticksvisible = true,
    #           xminorticks = IntervalsBetween(5), yminorticksvisible = true,
    #           yminorticks = IntervalsBetween(5))
    # hlines!(ax, [0]; color = (:black, 0.5), linestyle = :dash, linewidth = 2)
    # vlines!(ax, [0, 0.25]; color = (:black, 0.5), linestyle = :dash, linewidth = 2)
    # for (i, (O_h, O_m)) in enumerate(zip(Ō_h, Ō_m))
    #     structure = metadata(O_h)[:structure]
    #     x = O_h .- O_m
    #     μ = dropdims(mean(x, dims = SessionID), dims = SessionID)
    #     σ = dropdims(std(x, dims = SessionID), dims = SessionID)
    #     σ = σ ./ sqrt(size(O_h, 2)) # SEM
    #     bargs = [times(μ), μ .- σ, μ .+ σ] .|> ustripall .|> collect
    #     band!(ax, bargs..., color = (structurecolors[i], 0.3))
    #     lines!(ax, μ |> ustripall, color = (structurecolors[i], 0.7), label = structure)
    # end
    # axislegend(ax, position = :rb, framevisible = true, labelsize = 12)
    # display(f)
    # wsave(plotdir("theta_waves_task", "regional_orderparameter.pdf"), f)
    fig
end

# ---------------------------------- Hit/miss prediction --------------------------------- #
begin # * Filter out poor-performing mice and regenerate order parameters
    """
    f is a function applied to every element of the varaible var (e.g. `sign` for the
    direction of the phase velocity var=:k). One of these will take about 1 minute over 64 cores
    """
    # function filter_out(out, newsessions, var = :k, f = sign)
    #     out = [filter(x -> x.sessionid ∈ newsessions, o) for o in out]

    #     Og = progressmap(out; parallel = true) do o
    #         O = map(o) do p
    #             k = p[var]
    #             k = ZScore(k, dims = 1)(k)
    #             ret = dropdims(mean(f.(k), dims = Depth), dims = Depth)
    #             trials = p[:trials]
    #             trialtimes = trials.change_time_with_display_delay
    #             @assert all(isapprox.(ustripall(lookup(k, :changetime)), trialtimes,
    #                                   atol = 0.05)) # These won't match exactly, because the data change times have been adjusted for rectification. This is OK.
    #             ret = set(ret, Dim{:changetime} => Trial)
    #             set(ret,
    #                 DimensionalData.format(Trial(trials.hit), lookup(ret, Trial)))
    #         end
    #     end

    #     Og_h = [[o[:, lookup(o, Trial) .== true] for o in O] for O in Og]
    #     Og_m = [[o[:, lookup(o, Trial) .== false] for o in O] for O in Og]
    #     sessionids = [o[:sessionid] for o in out[1]]

    #     return Og, Og_h, Og_m, sessionids
    # end
    # Og, Og_h, Og_m, sessionids = filter_out(out, newsessions, :k, sign) # mean order parameter
    # GC.gc()
    # Ol, Ol_h, Ol_m, _ = filter_out(out, newsessions, :x, identity) # Mean theta LFP
    # GC.gc()

    Ol = map(out) do o
        O = map(o) do p
            k = p[:x]
            k = ZScore(k, dims = 1)(k)
            ret = dropdims(nansafe(mean; dims = Depth)(k), dims = Depth)
            trials = p[:trials][1:size(k, :changetime), :]
            trialtimes = trials.change_time_with_display_delay
            @assert all(isapprox.(ustripall(lookup(k, :changetime)), trialtimes,
                                  atol = 0.05)) # These won't match exactly, because the data change times have been adjusted for rectification. This is OK.
            set(ret, Dim{:changetime} => Trial(trials.hit))
        end
    end

    Ol_h = [[o[:, lookup(o, Trial) .== true] for o in O] for O in Ol]
    Ol_m = [[o[:, lookup(o, Trial) .== false] for o in O] for O in Ol]
end

# begin # Delete
#     f = Figure()
#     ax = Axis(f[1, 1])
#     a = 1
#     b = 22
#     m = mean(Og_h[a][b], dims = 2)[:]
#     sigma = std(Og_h[a][b], dims = 2) ./ 2
#     band!(1:length(m), m .- sigma[:], m .+ sigma[:]; color = (:cornflowerblue, 0.3))
#     lines!(m)
#     m = mean(Og_m[a][b], dims = 2)[:]
#     sigma = std(Og_m[a][b], dims = 2) ./ 2
#     band!(1:length(m), m .- sigma[:], m .+ sigma[:], color = (:crimson, 0.3))
#     lines!(m)
#     current_figure()
# end

begin # * LDA input data and downsampling
    H = [getindex.(Og, s) for s in eachindex(Og[1])] # Order parameter
    Hl = [getindex.(Ol, s) for s in eachindex(Ol[1])] # Mean LFP

    # * Make sure temporal dims same size
    ts = intersect([intersect(lookup.(h, 𝑡)...) for h in H]...)
    H = [[_h[𝑡 = At(ts)] for _h in h] for h in H]
    Hl = [[_h[𝑡 = At(ts)] for _h in h] for h in Hl]
    H = stack.([Structure(structures)], H, dims = 3)
    Hl = stack.([Structure(structures)], Hl, dims = 3)
    H = permutedims.(H, [(1, 3, 2)])
    Hl = permutedims.(Hl, [(1, 3, 2)])

    # H = [h[1:10:end, :, :] for h in H]
    # Hl = [h[1:10:end, :, :] for h in Hl]
    cg = h -> coarsegrain(h; dims = 1, newdim = 4)
    H = [dropdims(nansafe(mean; dims = 4)((cg ∘ cg ∘ cg)(h)); dims = 4) for h in H]
    @assert maximum(sum.([isnan.(x) for x in H])) .< 10
    map(H) do h
        for x in eachslice(h; dims = (Structure, Trial))
            while !isempty(findall(isnan.(x)))
                idxs = findall(isnan.(x))
                x[idxs] .= x[idxs .- 1]
            end
        end
    end
    @assert maximum(sum.([isnan.(x) for x in H])) .== 0
    Hl = [dropdims(nansafe(mean; dims = 4)((cg ∘ cg ∘ cg)(h)); dims = 4) for h in Hl]
    @assert maximum(sum.([isnan.(x) for x in Hl])) .< 10
    map(Hl) do h
        for x in eachslice(h; dims = (Structure, Trial))
            while !isempty(findall(isnan.(x)))
                idxs = findall(isnan.(x))
                x[idxs] .= x[idxs .- 1]
            end
        end
    end
    @assert maximum(sum.([isnan.(x) for x in Hl])) .== 0
    GC.gc()
end

if !isfile(datadir("hyperparameters", "theta_waves_task.jld2")) ||
   !isfile(datafile) # Run calculations; needs to be on a cluster
    if !isfile(datadir("hyperparameters", "theta_waves_task.jld2"))
        if haskey(ENV, "JULIA_DISTRIBUTED")
            using USydClusters
            procs = USydClusters.Physics.addprocs(32; mem = 15, ncpus = 2,
                                                  project = projectdir())
            @everywhere using SpatiotemporalMotifs
        else
            error("Calculations must be run on a cluster, set ENV[\"JULIA_DISTRIBUTED\"] to confirm this.")
        end
        begin # * Get a rough estimate for a good hyperparameter. Currently on pre-offset data. This gives us ~0.5 as a good estimate
            hfile = datadir("hyperparameters", "theta_waves_task.jld2")
            hyperr = pmap(SpatiotemporalMotifs.tuneclassifier, H)
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
    end

    begin # * Single-subject classifications, returning 5-fold balanced accuracy. Takes ages, about 1 hour
        regcoef = 0.5
        folds = 5
        repeats = 10

        bac_pred = pmap(H) do h
            h = h[𝑡 = -0.25u"s" .. 0.0u"s"]
            bac = classify_kfold(h; regcoef, k = folds, repeats)
        end

        bac_pre = pmap(H) do h
            h = h[𝑡 = -0.25u"s" .. 0.25u"s"]
            bac = classify_kfold(h; regcoef, k = folds, repeats)
        end

        bac_post = pmap(H) do h
            h = h[𝑡 = 0.25u"s" .. 0.75u"s"]
            bac = classify_kfold(h; regcoef, k = folds, repeats)
        end

        bac_sur = pmap(H) do h
            h = h[𝑡 = -0.25u"s" .. 0.25u"s"]
            idxs = randperm(size(h, Trial))
            h = set(h, Trial => lookup(h, Trial)[idxs])
            bac = classify_kfold(h; regcoef, k = folds, repeats)
        end

        bac_lfp = pmap(Hl) do h # * Mean LFP, pre-offset
            h = h[𝑡 = -0.25u"s" .. 0.25u"s"]
            bac = classify_kfold(h; regcoef, k = folds, repeats)
        end

        D = @strdict bac_pred bac_lfp bac_pre bac_post bac_sur regcoef folds repeats
    end

    begin # * Map of region-wise weightings
        W = map(H) do h
            h = h[𝑡 = -0.25u"s" .. 0.75u"s"] # !!!
            N, M = classifier(h; regcoef = 0.5) # !!!
            W = projection(M)
            W = reshape(W, size(h)[1:2])
            return ToolsArray(W, dims(h)[1:2])
        end
        W = stack(SessionID(sessionids), W, dims = 3)
        W = W ./ maximum(abs.(W))
    end
    push!(D, "W" => W)
    tagsave(datafile, D)
end

if isfile(datafile)
    D = jldopen(datafile, "r")
    @unpack bac_pre, bac_post, bac_sur, bac_lfp, W = D
    close(D)
else
    error("Please download $datafile or produce it")
end

begin # * Plot classification performance
    # ax = Axis(gs[3]; xlabel = "Balanced accuracy", ylabel = "Density",
    #   title = "Hit/miss classification")
    # ziggurat!(ax, bac_sur; bins = 10, label = "Null", normalizaton = :pdf, color = crimson)
    # ziggurat!(ax, bac_lfp; bins = 10, label = "Null", normalizaton = :pdf, color = crimson)
    # ziggurat!(ax, bac_post; bins = 10, label = "Post-offset", normalizaton = :pdf,
    #           color = cucumber)
    # ziggurat!(ax, bac_pre; bins = 10, label = "Pre-offset", normalizaton = :pdf,
    #           color = cornflowerblue)
    # axislegend(ax, position = :lt, framevisible = true, labelsize = 12)
    # tightlimits!(ax)
    ax = Axis(gs[3][1, 1],
              xticks = (1:4, ["Post-offset", "Pre-offset", "Mean LFP", "Null"]),
              ylabel = "Balanced accuracy", title = "Hit/miss classification",
              limits = (nothing, (0.35, 0.85)))
    boxargs = (; width = 0.75, strokewidth = 5, whiskerwidth = 0.2,
               strokecolor = (:gray, 0.0)) # !!!! Show outliers??
    boxplot!(ax, fill(1, length(bac_post)), bac_post; boxargs...)
    boxplot!(ax, fill(3, length(bac_lfp)), bac_lfp; boxargs...)
    boxplot!(ax, fill(2, length(bac_pre)), bac_pre; boxargs...)
    boxplot!(ax, fill(4, length(bac_sur)), bac_sur; color = :gray, boxargs...)
end

begin # * Plot region-wise weightings
    ax = Axis(gs[4]; xlabel = "Time (s)", ylabel = "Normalized LDA weight",
              title = "Regional classification weights")
    vlines!(ax, [0, 0.25], color = (:black, 0.2), linestyle = :dash)
    hlines!(ax, [0], color = (:black, 0.5), linewidth = 2)
    for structure in reverse(collect(lookup(W, Structure)))
        ws = W[Structure = At(structure)]
        ws = upsample(ws, 5, 1)
        ts = ustripall(lookup(ws, 𝑡))
        μ = mean(ws, dims = SessionID) |> vec
        σ = std(ws, dims = SessionID) ./ sqrt(size(ws, SessionID)) |> ustripall |> vec
        band!(ax, ts, μ - σ, μ + σ; color = (structurecolormap[structure], 0.3),
              label = structure)
        lines!(ax, ts, μ; color = structurecolormap[structure], label = structure)
    end
    tightlimits!(ax)
    l = axislegend(ax, position = :lt, nbanks = 2, framevisible = true, labelsize = 12,
                   merge = true)
    reverselegend!(l)
end

# -------------------------------- Rection-time prediction ------------------------------- #
# if false
#     begin # * K-fold validation with Huber (outlier-robust) regression
#         reaction_times = [o[:trials].lick_latency for o in first(out)]
#         regcoef = 1e5
#         bac = pmap(eachindex(H)) do i
#             regress_kfold(H[i][𝑡 = -0.0u"s" .. 0.25u"s"], reaction_times[i]; regcoef,
#                           k = 5)
#         end
#         bac_sur = pmap(eachindex(H)) do i
#             regress_kfold(H[i][𝑡 = -0.0u"s" .. 0.25u"s"],
#                           reaction_times[i][randperm(length(reaction_times[i]))]; regcoef,
#                           k = 5)
#         end
#     end
#     begin # * Plot regression scores
#         ax = Axis(gs[4]; xlabel = "Spearman correlation", ylabel = "Density",
#                   title = "Reaction time prediction")
#         ziggurat!(ax, bac_sur; bins = 10, label = "Null", normalizaton = :pdf,
#                   color = crimson)
#         # ziggurat!(ax, bac_post; bins = 10, label = "Post-stimulus", normalizaton = :pdf,
#         #           color = cucumber)
#         ziggurat!(ax, bac; bins = 10, label = "Pre-offset", normalizaton = :pdf,
#                   color = cornflowerblue)
#         axislegend(ax, position = :lt, framevisible = true, labelsize = 12)
#         tightlimits!(ax)
#     end
# end

begin # * Save
    addlabels!(fig, ["(d)", "(e)", "(f)", "(g)"])
    wsave(plotdir("theta_waves_task", "theta_waves_task.pdf"), fig)
    display(fig)
end
