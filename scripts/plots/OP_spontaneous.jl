#! /bin/bash
#=
exec julia +1.10.9 -t auto "${BASH_SOURCE[0]}" "$@"
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

            out = Iterators.flatten(parent.(out))
            out = collect(out)
            return out
        end
        return (@strdict sR)
    end
end

begin # * Plot histograms
    f = OnePanel()
    ax = Axis(f[1, 1])
    structureidx = findfirst(lookup(_Q, Structure) .== ["VISl"])
    map(data, stimuli, eachindex(stimuli)) do d, stimulus, i
        sR = d["sR"][structureidx]

        hill!(ax, sR, bandwidth = 0.05, facealpha = 0.1, boundary = (-1, 1),
              label = string(stimulus))
    end
    axislegend(ax, position = :lt, title = "Stimulus")
    display(f)
end
