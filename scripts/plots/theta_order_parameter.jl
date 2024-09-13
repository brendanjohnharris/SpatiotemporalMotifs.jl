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
    fig
end

# ---------------------------------- Hit/miss prediction --------------------------------- #
begin # * Generate mean LFP responses
    layergroups = [[1, 2], [3], [4, 5]] # Superficial, middle (L4), deep
    Olfp = map(layergroups) do ls
        map(out) do o
            map(o) do p
                k = p[:x]
                idxs = parselayernum.(metadata(k)[:layernames]) .∈ [ls]
                k = k[:, parent(idxs), :]
                k = ZScore(k, dims = 1)(k)
                ret = dropdims(nansafe(mean; dims = Depth)(k), dims = Depth)
                trials = p[:trials][1:size(k, :changetime), :]
                trialtimes = trials.change_time_with_display_delay
                @assert all(isapprox.(ustripall(lookup(k, :changetime)), trialtimes,
                                      atol = 0.05)) # These won't match exactly, because the data change times have been adjusted for rectification. This is OK.
                set(ret, Dim{:changetime} => Trial(trials.hit))
            end
        end
    end
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
    Hlfp = map(Olfp) do O
        [getindex.(O, s) for s in eachindex(O[1])]
    end

    # * Make sure temporal dims same size
    ts = intersect([intersect(lookup.(h, 𝑡)...) for h in H]...)
    H = [[_h[𝑡 = At(ts)] for _h in h] for h in H]
    H = stack.([Structure(structures)], H, dims = 3)
    H = permutedims.(H, [(1, 3, 2)])
    Hlfp = map(Hlfp) do H
        H = [[_h[𝑡 = At(ts)] for _h in h] for h in H]
        H = stack.([Structure(structures)], H, dims = 3)
        H = permutedims.(H, [(1, 3, 2)])
    end

    # H = [h[1:10:end, :, :] for h in H]
    # Hl = [h[1:10:end, :, :] for h in Hl]
    cg = h -> coarsegrain(h; dims = 1, newdim = 4)
    function consist(H)
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
        return H
    end
    H = consist(H)
    Hlfp = consist.(Hlfp)
    GC.gc()
end

if !isfile(datadir("hyperparameters", "theta_waves_task.jld2")) ||
   !isfile(datafile) # Run calculations; needs to be on a cluster
    if !isfile(datadir("hyperparameters", "theta_waves_task.jld2"))
        if haskey(ENV, "JULIA_DISTRIBUTED")
            using USydClusters
            procs = USydClusters.Physics.addprocs(26; mem = 22, ncpus = 4,
                                                  project = projectdir()) # ? Can reuse these for the following bac calculations
            @everywhere using SpatiotemporalMotifs
            @everywhere SpatiotemporalMotifs.@preamble
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
        regcoefs = first.(hyperr)
        folds = 5
        repeats = 10

        # bac_pred = pmap(H) do h
        #     h = h[𝑡 = -0.25u"s" .. 0.0u"s"]
        #     bac = classify_kfold(h; regcoef, k = folds, repeats)
        # end

        bac_pre = pmap(H, regcoefs) do h, regcoef
            h = h[𝑡 = -0.25u"s" .. 0.25u"s"]
            bac = classify_kfold(h; regcoef, k = folds, repeats)
        end

        bac_post = pmap(H, regcoefs) do h, regcoef
            h = h[𝑡 = 0.25u"s" .. 0.75u"s"]
            bac = classify_kfold(h; regcoef, k = folds, repeats)
        end

        bac_sur = pmap(H, regcoefs) do h, regcoef
            h = h[𝑡 = -0.25u"s" .. 0.25u"s"]
            idxs = randperm(size(h, Trial))
            h = set(h, Trial => lookup(h, Trial)[idxs])
            bac = classify_kfold(h; regcoef, k = folds, repeats)
        end

        bac_lfp = map(Hlfp) do H
            pmap(H, regcoefs) do h, regcoef # * Mean LFP, pre-offset
                h = h[𝑡 = -0.25u"s" .. 0.25u"s"]
                bac = classify_kfold(h; regcoef, k = folds, repeats)
            end
        end

        D = @strdict bac_lfp bac_pre bac_post bac_sur regcoefs folds repeats
    end

    begin # * Map of region-wise weightings
        W = map(H, regcoefs) do h, regcoef
            h = h[𝑡 = -0.25u"s" .. 0.75u"s"] # !!!
            N, M = classifier(h; regcoef) # !!!
            W = projection(M)
            W = reshape(W, size(h)[1:2])
            return ToolsArray(W, dims(h)[1:2])
        end
        W = stack(SessionID(oursessions), W, dims = 3)
        W = W ./ maximum(abs.(W)) # ! Normalized
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
    ax = Axis(gs[3][1, 1],
              xticks = (1:4, ["Post-offset", "Pre-offset", "Mean LFP", "Null"]),
              ylabel = "Balanced accuracy", title = "Hit/miss classification",
              limits = (nothing, (0.35, 0.85)))
    boxargs = (; width = 0.75, strokewidth = 5, whiskerwidth = 0.2,
               strokecolor = (:gray, 0.0)) # !!!! Show outliers??
    boxplot!(ax, fill(1, length(bac_post)), bac_post; boxargs...)
    boxplot!(ax, fill(2, length(bac_pre)), bac_pre; boxargs...)
    boxplot!(ax, vcat([fill(3 + i, length(bac_lfp[1])) for i in [-0.3, 0, 0.3]]...),
             vcat(bac_lfp...); boxargs..., width = 0.3)
    text!(ax, 3 .+ [-0.3, 0, 0.3], [0.8, 0.8, 0.8]; text = ["S", "M", "D"],
          align = (:center, :center))
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

begin # * Save
    addlabels!(fig, ["(d)", "(e)", "(f)", "(g)"])
    wsave(plotdir("theta_waves_task", "theta_waves_task.pdf"), fig)
    display(fig)
end
