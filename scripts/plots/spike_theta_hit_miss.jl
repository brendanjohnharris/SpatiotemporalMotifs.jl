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

begin # * So we want to load in the i) theta phases for each ii) spike time at iii) each depth, then do some comparison of these phases across layers for hit vs miss trials...
    stimulus = r"Natural_Images"
    vars = [:ϕ, :r]

    session_table = load(datadir("posthoc_session_table.jld2"), "session_table")
    oursessions = session_table.ecephys_session_id

    path = datadir("calculations")
    Q = calcquality(path)[Structure = At(structures)][SessionID = At(oursessions)]
    out = load_calculations(Q; stimulus, vars)
    Qs = calcquality(datadir("power_spectra"))[Structure = At(structures)][SessionID = At(oursessions)]
    unitdepths = load_unitdepths(Qs)

    begin # * Format spikes to depth dataframe. These have already been rectified
        spikes = map(out) do o
            regionalspikes = map(o) do _o
                DataFrame(collect.([keys(_o[:spiketimes]), values(_o[:spiketimes])]),
                          [:ecephys_unit_id, :spiketimes])
            end
            vcat(regionalspikes...)
        end
        D = vcat(spikes...)
        spikes = innerjoin(D, unitdepths, on = :ecephys_unit_id)
        filter!(:streamlinedepth => ∈(0 .. 1), spikes) # Cortex only
    end
end
