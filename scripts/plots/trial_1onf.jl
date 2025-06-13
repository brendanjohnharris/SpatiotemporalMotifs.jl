#! /bin/bash
# -*- mode: julia -*-
#=
exec julia +1.10.9 -t auto --color=yes "${BASH_SOURCE[0]}" "$@"
=#
using DrWatson
@quickactivate "SpatiotemporalMotifs"

using SpatiotemporalMotifs
import TimeseriesTools: freqs
using SpatiotemporalMotifs

using Random
@preamble
set_theme!(foresight(:physics))

begin # * Load the trial LFP for natural images
    session_table = load(datadir("posthoc_session_table.jld2"), "session_table")
    oursessions = session_table.ecephys_session_id
    path = datadir("calculations")
    Q = calcquality(path)[Structure = At(structures)]
    Q = Q[SessionID(At(oursessions))]
    @assert mean(Q[stimulus = At(r"Natural_Images")]) == 1
    out = load_calculations(Q; stimulus = r"Natural_Images", vars = [:V])
end

begin # * Calculate 1/f exponent for each trial
    function trial_1onf(x::UnivariateRegular)
    end
end

begin
    x = out[1][1][:V][:, 10, 10]
end

begin # * 1/f across all sessions
    oneonf = map(out) do out_structure
        map(out_structure) do out_session
            V = out_session[:V]
            V = V[𝑡 = SpatiotemporalMotif.INTERVAL]
        end
    end
end
