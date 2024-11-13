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
datafile = datadir("theta_order_parameter.jld2")
session_table = load(datadir("posthoc_session_table.jld2"), "session_table")
oursessions = session_table.ecephys_session_id
path = datadir("calculations")

begin # * Use extra workers if we can
    if haskey(ENV, "JULIA_DISTRIBUTED") && length(procs()) == 1
        using USydClusters
        ourprocs = USydClusters.Physics.addprocs(10; mem = 22, ncpus = 4,
                                                 project = projectdir()) # ? Lower this for a smaller cluster
        @everywhere using SpatiotemporalMotifs
        @everywhere SpatiotemporalMotifs.@preamble
    end
end

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

begin
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

        ax = Axis(gs[2]; xlabel = "Time (s)",
                  ylabel = rich("Mean order parameter ",
                                rich("R", subscript("θ"), font = "Times Italic")),
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

    begin # * Generate mean LFP responses
        layergroups = [1:5]#, [1, 2], [3], [4, 5]] # Superficial, middle (L4), deep
        Olfp = map(layergroups) do ls
            map(out) do o
                map(o) do p
                    k = deepcopy(p[:x])
                    idxs = parselayernum.(metadata(k)[:layernames]) .∈ [ls]
                    k = k[:, parent(idxs), :]
                    N = ZScore(k, dims = 1)
                    normalize!(k, N)
                    ret = dropdims(mean(k; dims = Depth), dims = Depth)
                end
            end
        end
    end
    begin # * Plot mean LFP
        Ō = mean.(Olfp[1])
        f = Figure()
        ax = Axis(f[1, 1]; xlabel = "Time (s)",
                  ylabel = rich("Mean order parameter ",
                                rich("R", subscript("θ"), font = "Times Italic")),
                  title = "Order parameter during flashes",
                  xautolimitmargin = (0, 0), xminorticksvisible = true,
                  xminorticks = IntervalsBetween(5), yminorticksvisible = true,
                  yminorticks = IntervalsBetween(5))
        hlines!(ax, [0]; color = (:black, 0.5), linestyle = :dash, linewidth = 2)
        vlines!(ax, [0, 0.25]; color = (:black, 0.5), linestyle = :dash, linewidth = 2)
        for (i, O) in reverse(collect(enumerate(Ō)))
            structure = metadata(O)[:structure]
            O = O[𝑡(SpatiotemporalMotifs.INTERVAL)]
            μ = dropdims(mean(O, dims = 2), dims = 2)
            σ = dropdims(std(O, dims = 2), dims = 2)
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
        f
    end
end

# * Task stimulus
stimulus = r"Natural_Images"
Q = calcquality(path)[Structure = At(structures)]
quality = mean(Q[stimulus = At(stimulus)])
out = load_calculations(Q; stimulus, vars)

begin # * Calculate a global order parameter at each time point
    map(out) do o
        filter!(o) do _o
            begin # * Remove trials with a reaciton time < 0.25s (so we can reliably divide trial periods into 'pre-reaction' and 'post-reaction' based on the timestamps. There are typically only a few of these per session.
                reaction_times = _o[:trials].lick_latency
                idxs = (reaction_times .> 0.25) .| isnan.(reaction_times) # NaN means a miss trial

                @assert issorted(lookup(_o[:k], :changetime))
                _o[:x] = _o[:x][:, :, idxs] # Remove trials with low reaction time
                _o[:k] = _o[:k][:, :, idxs]
                _o[:ω] = _o[:ω][:, :, idxs]
                _o[:trials] = _o[:trials][idxs, :]
            end
            _o[:sessionid] ∈ oursessions
        end
    end
    Og = map(out) do o
        O = map(o) do p
            k = p[:k]
            ω = p[:ω]
            k[ustripall(ω) .< 0] .= NaN * unit(eltype(k))
            # idxs = parselayernum.(metadata(k)[:layernames]) .∈ [[2, 3, 4]] # Just select middle layers
            # k = k[:, parent(idxs), :]
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

    ax = Axis(gs[1]; xlabel = "Time (s)",
              ylabel = rich("Mean order parameter ",
                            rich("R", subscript("θ"), font = "Times Italic")),
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

begin # * Order parameter correlation to hierarchy
    N = 1e6
    statsfile = plotdir("theta_order_parameter", "theta_order_parameter.txt")
    close(open(statsfile, "w")) # Create the file or clear it

    x = getindex.([SpatiotemporalMotifs.hierarchy_scores], structures)
    _y = stack(Structure(structures), Ō)
    _y = stack(Depth([0]), [_y])

    for T in [0.0u"s" .. 0.15u"s", 0.25u"s" .. 0.4u"s"]
        y = _y[𝑡 = T]
        y = dropdims(mean(y; dims = 𝑡); dims = 𝑡)
        open(statsfile, "a+") do file
            write(file, "\n# $T\n")
        end

        # * Group level
        μ, σ, 𝑝 = hierarchicalkendall(x, y, :group; N) .|> first
        open(statsfile, "a+") do file
            write(file, "\n## Group level")
            write(file, "\nmedian τ = $μ")
            write(file, "\n95\\% conf. = $σ")
            write(file, "\n𝑝 = $𝑝")
            write(file, "\nN = $N")
            write(file, "\n")
        end

        # * Individual level
        μ, σ, 𝑝 = hierarchicalkendall(x, y, :individual; N) .|> first
        open(statsfile, "a+") do file
            write(file, "\n## Individual level")
            write(file, "\nmedian τ = $μ")
            write(file, "\nIQR = $(σ[2] - σ[1])")
            write(file, "\n𝑝 = $𝑝")
            write(file, "\n")
        end
    end
end

# ---------------------------------- Hit/miss prediction --------------------------------- #
begin # * Generate mean LFP responses
    layergroups = [1:5, [1, 2], [3], [4, 5]] # Superficial, middle (L4), deep
    Olfp = map(layergroups) do ls
        map(out) do o
            map(o) do p
                k = deepcopy(p[:x])
                idxs = parselayernum.(metadata(k)[:layernames]) .∈ [ls]
                k = k[:, parent(idxs), :]
                N = ZScore(k, dims = 1)
                normalize!(k, N)
                ret = dropdims(mean(k; dims = Depth), dims = Depth)
                trials = p[:trials][1:size(k, :changetime), :]
                trialtimes = trials.change_time_with_display_delay
                @assert all(isapprox.(ustripall(lookup(k, :changetime)), trialtimes,
                                      atol = 0.05)) # These won't match exactly, because the data change times have been adjusted for rectification. This is OK.
                set(ret, Dim{:changetime} => Trial(trials.hit))
            end
        end
    end
end

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

    cg = h -> coarsegrain(h; dims = 1, newdim = 4)
    function consist(H)
        HH = pmap(H) do h
            dropdims(nansafe(mean; dims = 4)((cg ∘ cg ∘ cg)(h)); dims = 4)
        end
        nanprop = maximum(sum.([isnan.(x) for x in HH]) ./ length.([isnan.(x) for x in HH]))
        @info "$(round(nanprop*100, sigdigits=2))% of order parameters are NaN"
        @assert nanprop < 0.05

        map(HH) do h
            h[isnan.(h)] .= 0.0 # Ok. Could be better.
        end
        @assert maximum(sum.([isnan.(x) for x in HH])) .== 0
        return HH
    end
    H = consist(H)
    Hlfp = consist.(Hlfp)
    GC.gc()
end

if !isfile(datadir("hyperparameters", "theta_waves_task.jld2")) ||
   !isfile(datafile) # Run calculations; needs to be on a cluster
    if !isfile(datadir("hyperparameters", "theta_waves_task.jld2"))
        if !haskey(ENV, "JULIA_DISTRIBUTED")
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
        # regcoefs = first.(hyperr)
        regcoef = 0.5
        folds = 5
        repeats = 20

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

        bac_lfp_pre = map(Hlfp) do H
            pmap(H) do h # * Mean LFP, pre-offset
                h = h[𝑡 = -0.25u"s" .. 0.25u"s"]
                bac = classify_kfold(h; regcoef, k = folds, repeats)
            end
        end
        bac_lfp_post = map(Hlfp) do H
            pmap(H) do h # * Mean LFP, pre-offset
                h = h[𝑡 = 0.25u"s" .. 0.75u"s"]
                bac = classify_kfold(h; regcoef, k = folds, repeats)
            end
        end

        D = @strdict bac_pre bac_post bac_sur bac_lfp_pre bac_lfp_post regcoef folds repeats
    end

    begin # * Map of region-wise weightings
        W = pmap(H) do h
            h = h[𝑡 = -0.25u"s" .. 0.75u"s"] # !!!
            N, M = classifier(h; regcoef) # !!!
            W = projection(M)
            W = reshape(W, size(h)[1:2])
            return ToolsArray(W, dims(h)[1:2])
        end
        W = stack(SessionID(oursessions), W, dims = 3)
        W = W ./ maximum(abs.(W)) # ! Normalized

        Wlfp = map(Hlfp) do H
            W = pmap(H) do h
                h = h[𝑡 = -0.25u"s" .. 0.75u"s"] # !!!
                N, M = classifier(h; regcoef) # !!!
                W = projection(M)
                W = reshape(W, size(h)[1:2])
                return ToolsArray(W, dims(h)[1:2])
            end
            W = stack(SessionID(oursessions), W, dims = 3)
            W = W ./ maximum(abs.(W)) # ! Normalized
        end
    end
    push!(D, "W" => W)
    push!(D, "Wlfp" => Wlfp)
    tagsave(datafile, D)
end

if isfile(datafile)
    D = jldopen(datafile, "r")
    @unpack bac_pre, bac_post, bac_sur, bac_lfp_pre, bac_lfp_post, regcoef, folds, repeats, W, Wlfp = D
    close(D)
else
    error("Please download $datafile or produce it")
end

begin # * Plot classification performance
    ax = Axis(fig[2, 1];
              xticks = ([2, 5, 7], ["Post-offset", "Pre-offset", "Null"]),
              ylabel = "Balanced accuracy", title = "Hit/miss classification",
              limits = ((0.5, 7.5), (0.35, 0.95)), xminorticksvisible = false)
    boxargs = (; width = 0.75, strokewidth = 5, whiskerwidth = 0.2,
               strokecolor = (:gray, 0.0)) # !!!! Show outliers??

    vspan!(ax, 0.5, 3.5; color = (california, 0.2))
    vspan!(ax, 3.5, 6.5; color = (cucumber, 0.2))
    vspan!(ax, 6.5, 7.5; color = (:gray, 0.2))
    hlines!(ax, [0.5]; color = (:black, 0.5), linestyle = :dash, linewidth = 2)

    begin # * Write statistics
        open(statsfile, "a+") do file
            write(file, "\n# Classification performance\n")
            write(file, "\n## Order parameter pre-offset\n")
            write(file, "\nMedian BAc = $(median(bac_pre))")
            write(file, "\nIQR = $(iqr(bac_pre))")
            𝑝 = HypothesisTests.pvalue(HypothesisTests.MannWhitneyUTest(bac_pre, bac_sur);
                                       tail = :right)
            write(file, "\nU-test, right-sided 𝑝 to sur= $𝑝")
            write(file, "\n")

            write(file, "\n## Order parameter post-offset\n")
            write(file, "\nMedian BAc = $(median(bac_post))")
            write(file, "\nIQR = $(iqr(bac_post))")
            𝑝 = HypothesisTests.pvalue(HypothesisTests.MannWhitneyUTest(bac_pre, bac_post);
                                       tail = :both)
            write(file, "\nU-test, two-sided 𝑝 to pre. = $𝑝")
            write(file, "\n")

            write(file, "\n## Mean LFP pre-offset \n")
            write(file, "\nMedian BAc = $(median(bac_lfp_pre[1]))")
            write(file, "\nIQR = $(iqr(bac_lfp_pre[1]))")
            𝑝 = HypothesisTests.pvalue(HypothesisTests.MannWhitneyUTest(bac_pre,
                                                                        bac_lfp_pre[1]);
                                       tail = :right)
            write(file, "\nU-test, two-sided 𝑝 to pre. = $𝑝")
            write(file, "\n")

            write(file, "\n## Order parameter pre-offset null\n")
            write(file, "\nMedian BAc = $(median(bac_sur))")
            write(file, "\nIQR = $(iqr(bac_sur))")
            write(file, "\n")
        end
    end

    boxplot!(ax, fill(1, length(bac_post)), bac_post; boxargs..., color = cornflowerblue,
             label = "Order parameter")
    boxplot!(ax, fill(2, length(bac_lfp_post[1])), bac_lfp_post[1]; boxargs...,
             color = juliapurple, label = "Mean LFP")
    boxplot!(ax, vcat([fill(3 + i, length(bac_lfp_post[1])) for i in [-0.3, 0, 0.3]]...),
             reverse(vcat(bac_lfp_post[2:end]...)); boxargs..., width = 0.3,
             color = crimson, label = "Compartmental LFP")
    text!(ax, 3 .+ [-0.3, 0, 0.3], [0.9, 0.9, 0.9]; text = reverse(["S", "M", "D"]),
          align = (:center, :center))

    boxplot!(ax, fill(4, length(bac_pre)), bac_pre; boxargs..., color = cornflowerblue,
             label = "Order parameter")
    boxplot!(ax, fill(5, length(bac_lfp_pre[1])), bac_lfp_pre[1]; boxargs...,
             color = juliapurple, label = "Mean LFP")
    boxplot!(ax, vcat([fill(6 + i, length(bac_lfp_pre[1])) for i in [-0.3, 0, 0.3]]...),
             reverse(vcat(bac_lfp_pre[2:end]...)); boxargs..., width = 0.3, color = crimson,
             label = "Compartmental LFP")
    text!(ax, 6 .+ [-0.3, 0, 0.3], [0.8, 0.8, 0.8]; text = reverse(["S", "M", "D"]),
          align = (:center, :center))

    boxplot!(ax, fill(7, length(bac_sur)), bac_sur; color = :gray, boxargs...)
    axislegend(ax; merge = true)
end

begin # * Plot region-wise weightings
    ax = Axis(fig[2, 2]; xlabel = "Time (s)", ylabel = "Normalized LDA weight",
              title = "Regional classification weights")
    vlines!(ax, [0, 0.25], color = (:black, 0.2), linestyle = :dash)
    hlines!(ax, [0], color = (:black, 0.5), linewidth = 2)
    for structure in (collect(lookup(W, Structure)))
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
    # reverselegend!(l)
end

begin # * Save
    addlabels!(fig, ["(d)", "(e)", "(f)", "(g)"])
    wsave(plotdir("theta_order_parameter", "theta_order_parameter.pdf"), fig)
    display(fig)
end
