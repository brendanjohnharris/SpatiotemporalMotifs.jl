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

begin # * Load spontaneous order parameter
    pQ = calcquality(datadir("power_spectra"))
    stimuli = [r"Natural_Images", "flash_250ms", "spontaneous"]
    session_table = load(datadir("posthoc_session_table.jld2"), "session_table")
    oursessions = session_table.ecephys_session_id

    data = map(stimuli) do stimulus
        _Q = pQ[stimulus = At(stimulus), Structure = At(structures)]
        subsessions = intersect(oursessions, lookup(_Q, SessionID))
        if length(subsessions) < length(oursessions)
            @warn "Power spectra calculations are incomplete, proceeding regardless"
        end
        _Q = _Q[SessionID(At(subsessions))]
        filebase = stimulus == "spontaneous" ? "" : "_$stimulus"

        sR = map(lookup(_Q, Structure)) do structure # * Load data
            out = map(lookup(_Q, SessionID)) do sessionid
                if _Q[SessionID = At(sessionid), Structure = At(structure)] == 0
                    return nothing
                end
                filename = savepath((@strdict sessionid structure stimulus), "jld2",
                                    datadir("power_spectra"))
                sR = load(filename, "sR")[1:5:end] # Downsample to keep things manageable
                # S = load(filename, "sC")
                # return (C .- S) ./ median(S)
            end
            idxs = .!isnothing.(out)
            out = out[idxs]

            # out = Iterators.flatten(parent.(out))
            # out = collect(out)
            return out
        end
        return sR
    end
end

begin # * Extract the probability mass of 'coherent' events; |OP| > 0.75
    coherent_events = map(data) do sR
        out = map(sR) do y
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
end

begin # * Plot boxplots of probability of coherent events
    f = Figure()
    ax = Axis(f[1, 1], ylabel = "Probability of coherent events",
              xticks = (1.5:1:3.5, string.(stimuli)))
    levents = coherent_events[side = At(:l)]
    for (i, stimulus) in enumerate(lookup(levents, :stimulus))
        y = levents[stimulus = At(stimulus)]
        structures = lookup(y, Structure)
        x = repeat(eachindex(structures)', size(y, SessionID), 1)

        boxplot!(ax, i .+ x[:] ./ (length(structures) + 1), y[:], width = 0.1,
                 color = getindex.([structurecolors], x)[:], show_outliers = false)
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

begin # * Omission trials or passive trials?
    stimulus = "Natural_Images_passive_nochange"
    path = datadir("calculations")
    Q = calcquality(path)[Structure = At(structures)]
    quality = mean(Q[stimulus = At(stimulus)])
    vars = [:csd, :k, :ω]
    out = load_calculations(Q; stimulus, vars)

    session_table = load(datadir("posthoc_session_table.jld2"), "session_table")
    oursessions = session_table.ecephys_session_id

    begin # * Calculate a global order parameter and mean LFP at each time point
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
    omission = @strdict Og
end

begin # * Order parameter during omission/passive trials
    f = Figure()
    ax = Axis(f[1, 1], xlabel = "Time (s)", ylabel = "Order parameter")
    t = times(first(Og)[1])[1:1562]
    for i in 1:length(Og)
        x = parent.(Og[i])
        x = getindex.(x, (1:1562,), (:,))
        x = hcat(x...)
        x = mean(x, dims = 2)
        x = dropdims(x, dims = 2)
        lines!(ustripall(t), x, label = structures[i], color = structurecolors[i])
    end
    axislegend(ax, position = :rt, title = "Structure")
    ax.limits = ((-0.25, 0.75), (-0.4, 0.8))
    display(f)
end
# * Now plot this as a third panel in main text.....
