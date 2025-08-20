#! /bin/bash
#=
exec julia +1.10.10 -t auto "${BASH_SOURCE[0]}" "$@"
=#
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

plot_data, data_file = produce_or_load(Dict("thr" => 0.5), calcdir("plots");
                                       filename = savepath("figS3")) do config
    thr = config["thr"]
    begin # * Load active vs passive order parameters
        path = calcdir("calculations")
        vars = [:k, :ω]

        session_table = load(calcdir("plots", "posthoc_session_table.jld2"),
                             "session_table")
        oursessions = session_table.ecephys_session_id
        begin # * Active natural images
            stimulus = r"Natural_Images"
            Q = calcquality(path; require = ["k", "ω"])[Structure = At(structures)]
            quality = mean(Q[stimulus = At(stimulus)])
            out = load_calculations(Q; stimulus, vars, path)

            begin # * Calculate a global order parameter and mean LFP at each time point
                out = map(out) do o
                    filter(o) do _o
                        _o[:sessionid] ∈ oursessions
                    end
                end
                sessionids = map(first(out)) do o
                    metadata(o[:k])[:sessionid]
                end
                Og = map(out) do o
                    O = map(o) do p
                        k = p[:k]
                        ω = p[:ω]
                        k[ustripall(ω) .< 0] .= NaN * unit(eltype(k))
                        ret = dropdims(nansafe(mean; dims = Depth)(sign.(k)), dims = Depth) # Ignore negative frequencies
                        set(ret, Dim{:changetime} => Trial) .|> Float32
                    end
                    O = ToolsArray(O, (SessionID(sessionids),))
                end
                Og = ToolsArray(Og, (Structure(structures),))
            end
            active = Og
            out = []
            GC.gc()
        end

        begin # * Passive natural image changes
            stimulus = "Natural_Images_passive"
            Q = calcquality(path; require = ["k", "ω"])[Structure = At(structures)]
            quality = mean(Q[stimulus = At(stimulus)])
            out = load_calculations(Q; stimulus, vars, path)

            begin # * Calculate a global order parameter and mean LFP at each time point
                out = map(out) do o
                    filter(o) do _o
                        _o[:sessionid] ∈ oursessions
                    end
                end
                sessionids = map(first(out)) do o
                    metadata(o[:k])[:sessionid]
                end
                Og = map(out) do o
                    O = map(o) do p
                        k = p[:k]
                        ω = p[:ω]
                        k[ustripall(ω) .< 0] .= NaN * unit(eltype(k))
                        ret = dropdims(nansafe(mean; dims = Depth)(sign.(k)), dims = Depth) # Ignore negative frequencies
                        set(ret, Dim{:changetime} => Trial) .|> Float32
                    end
                    O = ToolsArray(O, (SessionID(sessionids),))
                end
                Og = ToolsArray(Og, (Structure(structures),))
            end
            passive = Og
            out = []
            GC.gc()
        end
    end

    begin # * Calculate contrast between active and passive mean order parameters
        ma = SpatiotemporalMotifs.ramap(AbstractArray{<:Number}, active) do x
                 x = nansafe(mean, dims = 2)(x)
                 dropdims(x, dims = Trial)[1:1562]
             end .|> stack |> stack
        mp = SpatiotemporalMotifs.ramap(AbstractArray{<:Number}, passive) do x
                 x = nansafe(mean, dims = 2)(x)
                 dropdims(x, dims = Trial)[1:1562]
             end .|> stack |> stack
        m = ma .- mp
    end

    begin # * Load spontaneous order parameter
        pQ = calcquality(calcdir("power_spectra"))
        stimuli = [r"Natural_Images", "flash_250ms", "spontaneous"]
        session_table = load(calcdir("plots", "posthoc_session_table.jld2"),
                             "session_table")
        oursessions = session_table.ecephys_session_id

        data = map(stimuli) do stimulus
            _Q = pQ[stimulus = At(stimulus), Structure = At(structures)]
            subsessions = intersect(oursessions, lookup(_Q, SessionID))
            if length(subsessions) < length(oursessions)
                @warn "Power spectra calculations are incomplete, proceeding regardless"
            end
            _Q = _Q[SessionID(At(subsessions))]
            filebase = stimulus == "spontaneous" ? "" : "_$stimulus"

            R = map(lookup(_Q, Structure)) do structure # * Load data
                out = map(lookup(_Q, SessionID)) do sessionid
                    if _Q[SessionID = At(sessionid), Structure = At(structure)] == 0
                        return nothing
                    end
                    filename = savepath((@strdict sessionid structure stimulus), "jld2",
                                        calcdir("power_spectra"))
                    R = load(filename, "R")[1:5:end] # Downsample to keep things manageable
                    sR = load(filename, "sR")[1:5:end] # Downsample to keep things manageable
                    # S = load(filename, "sC")
                    # return (C .- S) ./ median(S)
                    return R, sR
                end
                idxs = .!isnothing.(out)
                out = out[idxs]

                # out = Iterators.flatten(parent.(out))
                # out = collect(out)
                return out
            end
            return R
        end
    end

    begin # * Extract the probability mass of 'coherent' events; |OP| > 0.5
        coherent_events = map(data) do R
            R = SpatiotemporalMotifs.ramap(Base.Fix2(getindex, 1), R)
            out = map(R) do y
                out = map(y) do x
                    l = sum(x -> x .< -thr, x) / length(x)
                    u = sum(x -> x .> thr, x) / length(x)
                    return [l, u]
                end |> stack
            end |> stack
        end |> stack
        coherent_events_sur = map(data) do R
            R = SpatiotemporalMotifs.ramap(Base.Fix2(getindex, 2), R)
            out = map(R) do y
                out = map(y) do x
                    l = sum(x -> x .< -thr, x) / length(x)
                    u = sum(x -> x .> thr, x) / length(x)
                    return [l, u]
                end |> stack
            end |> stack
        end |> stack
        ds = (Dim{:side}([:l, :u]), SessionID(oursessions), Structure(structures),
              Dim{:stimulus}(stimuli))
        coherent_events = ToolsArray(coherent_events, ds)
        coherent_events_sur = ToolsArray(coherent_events_sur, ds)
    end

    return Dict("order_parameters" => Dict("active" => ma,
                                           "passive" => mp,
                                           "contrast" => m),
                "coherent_events" => coherent_events,
                "coherent_events_surrogate" => coherent_events_sur)
end

begin
    f = SixPanel()
    gs = subdivide(f, 3, 2)
end

begin # * Order parameter for passive natural images
    mp = plot_data["order_parameters"]["passive"]
    ax = Axis(gs[1]; xlabel = "Time (s)",
              ylabel = rich("Mean order parameter ",
                            rich("R", subscript("θ"), font = "Times Italic")),
              title = "Passive order parameter",
              xautolimitmargin = (0, 0), xminorticksvisible = true,
              xminorticks = IntervalsBetween(5), yminorticksvisible = true,
              yminorticks = IntervalsBetween(5),
              xtickformat = terseticks,
              ytickformat = terseticks)
    hlines!(ax, [0]; color = (:black, 0.5), linestyle = :dash, linewidth = 2)
    vlines!(ax, [0, 0.25]; color = (:black, 0.5), linestyle = :dash, linewidth = 2)
    for structure in reverse(structures)
        O = mp[Structure = At(structure)]
        O = O[𝑡(SpatiotemporalMotifs.INTERVAL)]
        μ = dropdims(mean(O, dims = SessionID), dims = SessionID)
        σ = dropdims(std(O, dims = SessionID), dims = SessionID)
        σ = σ ./ 2
        bargs = [times(μ), μ .- σ, μ .+ σ] .|> ustripall .|> collect
        band!(ax, bargs..., color = (structurecolormap[structure], 0.3), label = structure)
        lines!(ax, times(μ) |> ustripall, μ |> ustripall |> collect,
               color = (structurecolormap[structure], 0.7),
               label = structure)
    end
    l = axislegend(ax, position = :lt, nbanks = 2, framevisible = true,
                   labelsize = 12,
                   merge = true)
    reverselegend!(l)
end

begin # * Plot contrast of active and passive
    m = plot_data["order_parameters"]["contrast"]
    ax = Axis(gs[2]; xlabel = "Time (s)",
              ylabel = rich("Mean active - passive ",
                            rich("R", subscript("θ"), font = "Times Italic")),
              title = "Order parameter contrast", xautolimitmargin = (0, 0),
              xminorticksvisible = true,
              xminorticks = IntervalsBetween(5), yminorticksvisible = true,
              yminorticks = IntervalsBetween(5),
              xtickformat = terseticks,
              ytickformat = terseticks)

    hlines!(ax, [0]; color = (:black, 0.5), linestyle = :dash, linewidth = 2)
    vlines!(ax, [0, 0.25]; color = (:black, 0.5), linestyle = :dash, linewidth = 2)

    map(reverse(lookup(m, Structure))) do structure
        mstructure = m[Structure = At(structure), 𝑡(SpatiotemporalMotifs.INTERVAL)] |>
                     ustripall

        μ = nansafe(mean, dims = 2)(mstructure)
        μ = dropdims(μ, dims = SessionID)

        # σ = nansafe(std, dims = 2)(mstructure)
        # σ = dropdims(σ, dims = SessionID)
        # σl = μ .- σ / 2
        # σh = μ .+ σ / 2

        # μ, (σl, σh) = bootstrapaverage(mean, mstructure[1:2:end, :]; dims = 2)

        ts = collect(lookup(μ, 𝑡))
        # band!(ax, ts, parent(σl), parent(σh), label = structure,
        #   color = structurecolormap[structure], alpha = 0.2)
        lines!(ax, ts, parent(μ), label = structure,
               color = structurecolormap[structure])
    end
    # l = axislegend(ax, position = :lt, nbanks = 2, framevisible = true,
    #                labelsize = 12,
    #                merge = true)
    # reverselegend!(l)
    # display(f)
end

begin # * Plot boxplots of probability of coherent events
    coherent_events = plot_data["coherent_events"]
    coherent_events_sur = plot_data["coherent_events_surrogate"]
    stimuli = [r"Natural_Images", "flash_250ms", "spontaneous"]
    slabels = ["Natural_Images", "flash_250ms", "spontaneous"]
    for (i, side) in enumerate([:u, :l])
        for (j, surr) in enumerate([false, true])
            title = "Coherent events ("

            if surr
                title = title * "null, "
            end

            if side === :u
                title = title * "downward)"
            else
                title = title * "upward)"
            end

            ax = Axis(f[i + 1, j]; ylabel = "Probability of coherent events",
                      xticks = (1.5:1:3.5, string.(slabels)), title, xticklabelsize = 15,
                      ytickformat = terseticks)

            vlines!(ax, [2, 3]; color = (:black, 0.5), linewidth = 2)

            if surr
                levents = coherent_events_sur[side = At(side)]
            else
                levents = coherent_events[side = At(side)]
            end

            for (i, stimulus) in enumerate(lookup(levents, :stimulus))
                y = levents[stimulus = At(stimulus)]
                structures = lookup(y, Structure)[:]
                x = repeat(eachindex(structures)', size(y, SessionID), 1)[:]

                colors = getindex.([structurecolors], x)[:]

                alphacolors = map(colors) do c
                    Makie.coloralpha(convert(Makie.RGB, c), 0.4)
                end

                map(unique(x)) do _x
                    idxs = x .== _x
                    c = colors[idxs]
                    a = alphacolors[idxs]
                    boxplot!(ax, i .+ x[idxs] ./ (length(structures) + 1), y[idxs],
                             width = 0.075,
                             color = a,
                             show_outliers = false,
                             whiskerlinewidth = 1,
                             whiskercolor = (:black, 0.8),
                             mediancolor = c,
                             medianlinewidth = 7, label = structures[_x])
                end
                # ms = median(levents[stimulus = At(stimulus)], dims = SessionID)
                # ms = dropdims(ms, dims = SessionID)
                # scatter!(ax, i .+ range(0, 1, length = length(structures)), collect(ms))
            end
        end
    end
    # for (i, structure) in enumerate(lookup(levents, Structure))
    #     y = levents[Structure = At(structure)]
    #     stimuli = lookup(y, :stimulus)
    #     x = repeat(eachindex(stimuli)', size(y, SessionID), 1)

    #     boxplot!(ax, i .+ x[:] ./ (length(stimuli) + 1), y[:], width = 0.1,
    #              color = getindex.([structurecolors], x)[:], show_outliers = false)
    # end
    addlabels!(f, labelformat)
    display(f)
    wsave(plotdir("figS3", "figS3.pdf"), f)
end
