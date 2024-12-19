#! /bin/bash
# -*- mode: julia -*-
#=
exec julia -t auto --color=yes "${BASH_SOURCE[0]}" "$@"
=#
using DrWatson
@quickactivate "SpatiotemporalMotifs"

using SpatiotemporalMotifs
import TimeseriesTools: freqs
using Peaks
using FileIO
using SpatiotemporalMotifs
import SpatiotemporalMotifs.layers
@preamble
set_theme!(foresight(:physics))

stimuli = ["r\"Natural_Images\"", "spontaneous", "flash_250ms"]
xtickformat = x -> string.(round.(Int, x))

session_table = load(datadir("posthoc_session_table.jld2"), "session_table")
oursessions = session_table.ecephys_session_id
path = datadir("power_spectra")
Q = calcquality(path)

for stimulus in stimuli
    _Q = Q[stimulus = At(stimulus), Structure = At(structures)]
    subsessions = intersect(oursessions, lookup(_Q, SessionID))
    if length(subsessions) < length(oursessions)
        @warn "Power spectra calculations are incomplete, proceeding regardless"
    end
    _Q = _Q[SessionID(At(subsessions))]
    filebase = stimulus == "spontaneous" ? "" : "_$stimulus"
    f = Figure(size = (900, 1080))

    begin # * Load data
        S = map(lookup(_Q, Structure)) do structure
            out = map(lookup(_Q, SessionID)) do sessionid
                if _Q[SessionID = At(sessionid), Structure = At(structure)] == 0
                    return nothing
                end
                filename = savepath((@strdict sessionid structure stimulus), "jld2", path)
                C = load(filename, "C")
                # S = load(filename, "sC")
                # return (C .- S) ./ median(S)
            end
            idxs = .!isnothing.(out)
            out = out[idxs]

            m = DimensionalData.metadata.(out)
            out = map(out) do o
                dropdims(mean(o, dims = Chan), dims = Chan)
            end
            out = stack(SessionID(lookup(_Q, SessionID)[idxs]), out, dims = 3)
            return out
        end
    end

    begin # * Supplemental average comodulograms
        f = SixPanel()
        gs = subdivide(f, 3, 2)
        map(gs, structures, S) do g, structure, s
            ax = Axis(g[1, 1]; title = structure, xlabel = "Phase frequency (Hz)",
                      ylabel = "Amplitude frequency (Hz)")
            s = dropdims(mean(s, dims = SessionID); dims = SessionID)
            s = upsample(s, 5, 1)
            s = upsample(s, 5, 2)
            p = heatmap!(ax, s; colormap = seethrough(reverse(sunrise)), rasterize = 5)
            Colorbar(g[1, 2], p; label = "Modulation index")
        end
        addlabels!(f)
        display(f)
        wsave(plotdir("comodulograms", "comodulograms_$stimulus.pdf"), f)
    end
end
