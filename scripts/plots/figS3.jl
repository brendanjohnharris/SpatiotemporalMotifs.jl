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

begin # * Load active vs passive order parameters
    path = calcdir("calculations")
    vars = [:k, :ω]

    session_table = load(calcdir("posthoc_session_table.jld2"), "session_table")
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
    ma = ramap(AbstractArray{<:Number}, active) do x
             x = mean(x, dims = Trial)
             dropdims(x, dims = Trial)[1:1562]
         end .|> stack |> stack
    mp = ramap(AbstractArray{<:Number}, passive) do x
             x = mean(x, dims = Trial)
             dropdims(x, dims = Trial)[1:1562]
         end .|> stack |> stack
    m = ma .- mp
end

begin # * Plot contrast of active and passive
    f = Figure()
    ax = Axis(f[1, 1]; xlabel = "Time (s)", ylabel = "Order parameter contrast",
              title = "Active vs passive mean order parameter")

    map(lookup(m, Structure)) do structure
        mstructure = m[Structure = At(structure)] |> ustripall

        μ = nansafe(mean, dims = 2)(mstructure)
        μ = dropdims(μ, dims = SessionID)

        σ = nansafe(std, dims = 2)(mstructure)
        σ = dropdims(σ, dims = SessionID)
        σl = μ .- σ / 2
        σh = μ .+ σ / 2

        ts = collect(lookup(μ, 𝑡))
        # band!(ax, ts, parent(σl), parent(σh), label = structure,
        #       color = structurecolormap[structure], alpha = 0.2)
        lines!(ax, ts, parent(μ), label = structure,
               color = structurecolormap[structure])
    end
    display(f)
end

begin # * Load spontaneous order parameter
    pQ = calcquality(calcdir("power_spectra"))
    stimuli = [r"Natural_Images", "flash_250ms", "spontaneous"]
    session_table = load(calcdir("posthoc_session_table.jld2"), "session_table")
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

begin # * Extract the probability mass of 'coherent' events; |OP| > 0.75
    coherent_events = map(data) do R
        R = SpatiotemporalMotifs.ramap(Base.Fix2(getindex, 1), R)
        out = map(R) do y
            out = map(y) do x
                l = sum(x -> x .< -0.5, x) / length(x)
                u = sum(x -> x .> 0.5, x) / length(x)
                return [l, u]
            end |> stack
        end |> stack
    end |> stack
    coherent_events_sur = map(data) do R
        R = SpatiotemporalMotifs.ramap(Base.Fix2(getindex, 2), R)
        out = map(R) do y
            out = map(y) do x
                l = sum(x -> x .< -0.5, x) / length(x)
                u = sum(x -> x .> 0.5, x) / length(x)
                return [l, u]
            end |> stack
        end |> stack
    end |> stack
    ds = (Dim{:side}([:l, :u]), SessionID(oursessions), Structure(structures),
          Dim{:stimulus}(stimuli))
    coherent_events = ToolsArray(coherent_events, ds)
    coherent_events_sur = ToolsArray(coherent_events_sur, ds)
end

begin # * Plot boxplots of probability of coherent events
    f = SixPanel()

    gs = subdivide(f, 3, 2)

    for (i, side) in enumerate([:u, :l])
        for (j, surr) in enumerate([false, true])
            title = "Coherent events ("

            if surr
                title = title * "null, "
            end

            if side === :u
                title = title * "upward)"
            else
                title = title * "downward)"
            end

            ax = Axis(f[i + 1, j]; ylabel = "Probability of coherent events",
                      xticks = (1.5:1:3.5, string.(stimuli)), title)

            if surr
                levents = coherent_events_sur[side = At(side)]
            else
                levents = coherent_events[side = At(side)]
            end

            for (i, stimulus) in enumerate(lookup(levents, :stimulus))
                y = levents[stimulus = At(stimulus)]
                structures = lookup(y, Structure)
                x = repeat(eachindex(structures)', size(y, SessionID), 1)

                boxplot!(ax, i .+ x[:] ./ (length(structures) + 1), y[:], width = 0.1,
                         color = getindex.([structurecolors], x)[:], show_outliers = false,
                         whiskerlinewidth = 2)
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
    display(f)
end

# begin # * Plot histograms
#     f = OnePanel()
#     ax = Axis(f[1, 1])
#     structureidx = findfirst(lookup(_Q, Structure) .== ["VISl"])
#     map(data, stimuli, eachindex(stimuli)) do d, stimulus, i
#         sR = d["sR"][structureidx]

#         hill!(ax, sR, bandwidth = 0.05, facealpha = 0.1, boundary = (-1, 1),
#               label = string(stimulus))
#     end
#     axislegend(ax, position = :lt, title = "Stimulus")
#     display(f)
# end

begin # * Supplementary figure: contrast between active and passive stimulus
    begin # * Active rder parameter
        @unpack Og = plot_data["order_parameters"]["Natural_Images"]
        Ō_active = orderparameter(Og)

        @unpack Og = plot_data["order_parameters"]["passive"]
        Ō_passive = orderparameter(Og)

        begin # * Plot the mean order parameter across time
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
            l = axislegend(ax, position = :lt, nbanks = 2, framevisible = true,
                           labelsize = 12,
                           merge = true)
            reverselegend!(l)
        end
    end
end
