#! /bin/bash
# -*- mode: julia -*-
#=
exec julia +1.10.9 -t auto --color=yes "${BASH_SOURCE[0]}" "$@"
=#
using DrWatson
@quickactivate "SpatiotemporalMotifs"

using SpatiotemporalMotifs
import TimeseriesTools: freqs
using Peaks
using FileIO
using SpatiotemporalMotifs
import SpatiotemporalMotifs.layers
using Random
@preamble
set_theme!(foresight(:physics))

# begin # * So we want to load in the i) theta phases for each ii) spike time at iii) each depth, then do some comparison of these phases across layers for hit vs miss trials...
#     stimulus = r"Natural_Images"
#     vars = [:ϕ, :r]

#     session_table = load(datadir("posthoc_session_table.jld2"), "session_table")
#     oursessions = session_table.ecephys_session_id

#     path = datadir("calculations")
#     Q = calcquality(path)[Structure = At(structures)][SessionID = At(oursessions)]
#     out = load_calculations(Q; stimulus, vars)
#     Qs = calcquality(datadir("power_spectra"))[Structure = At(structures)][SessionID = At(oursessions)]
#     unitdepths = load_unitdepths(Qs)

#     begin # * Format spikes to depth dataframe. These have already been rectified
#         spikes = map(out) do o
#             regionalspikes = map(o) do _o
#                 DataFrame(collect.([keys(_o[:spiketimes]), values(_o[:spiketimes])]),
#                           [:ecephys_unit_id, :spiketimes])
#             end
#             vcat(regionalspikes...)
#         end
#         D = vcat(spikes...)
#         spikes = innerjoin(D, unitdepths, on = :ecephys_unit_id)
#         filter!(:streamlinedepth => ∈(0 .. 1), spikes) # Cortex only
#     end
# end

begin # * Load trial-by-trial spc
    file = datadir("spike_lfp.jld2")
    pspikes = load(file, "pspikes")
end

begin # * Compare distribution of PPC for hit/miss trials
    sessionids = unique(pspikes.ecephys_session_id)
    outfile = datadir("out&stimulus=Natural_Images.jld2")
    trials = jldopen(outfile, "r") do f
        map(sessionids) do sessionid
            trials = f["VISp/$(sessionid)/trials"]
        end
    end
    trials = ToolsArray(trials, (SessionID(sessionids),))
    pspikes.hitmiss = map(eachrow(pspikes)) do row
        isempty(row[:trial_pairwise_phase_consistency]) && return
        _trials = trials[SessionID = At(row[:ecephys_session_id])]
        return _trials.hit
    end
end

begin # * Plot mean layer-wise PPC for hit/miss trials
    sessionid = sessionids[1]
    x = pspikes[pspikes.ecephys_session_id .== sessionid, :]
    structure = "VISl"
    x = x[x.structure_acronym .== structure, :]
    x = x[0.0 .< x.streamlinedepth .< 0.4, :]

    ppcs = x.trial_pairwise_phase_consistency
    hitmiss = x.hitmiss
    hits = map(ppcs, hitmiss) do p, h
        filter(!isnan, p[h])
    end
    misses = map(ppcs, hitmiss) do p, h
        filter(!isnan, p[.!h])
    end
end
begin
    density(filter(!isnan, vcat(hits...)))
    vlines!(mean(filter(!isnan, vcat(hits...))), color = cornflowerblue, label = "mean hit")
    density!(filter(!isnan, vcat(misses...)))
    vlines!(mean(filter(!isnan, vcat(misses...))), color = crimson, label = "mean miss")
    axislegend(; position = :lt)
    current_axis().xlabel = "PPC"
    current_axis().ylabel = "Density"
    current_figure() |> display
end
