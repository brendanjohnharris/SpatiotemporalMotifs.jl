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
using Distributed

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
    addprocs(50)
    @everywhere using DrWatson
    @everywhere DrWatson.@quickactivate "SpatiotemporalMotifs"
    @everywhere using Unitful
    @everywhere using TimeseriesTools
    @everywhere using SpatiotemporalMotifs

    x = out[1][1][:V][:, 10, 10]
    @spawnat 5 SpatiotemporalMotifs.trial_1onf(x)
end

begin # * 1/f across all sessions. Should take about 15 minutes over 50 workers
    oneoneff = map(out) do out_structure
        @info "Calculating 1/f exponent for $(metadata(out_structure[1][:V])[:structure])"
        pmap(out_structure) do out_session
            V = out_session[:V]
            V = V[𝑡 = SpatiotemporalMotifs.INTERVAL]
            map(SpatiotemporalMotifs.trial_1onf, eachslice(V, dims = (:Depth, :changetime)))
        end
    end
end
begin
    χ = map(oneoneff) do oneoneff_structure
        map(oneoneff_structure) do oneoneff_session
            getindex.(oneoneff_session, :χ)
        end
    end
    b = map(oneoneff) do oneoneff_structure
        map(oneoneff_structure) do oneoneff_session
            getindex.(oneoneff_session, :b)
        end
    end
    D = @strdict χ b
    tagsave(datadir("trial_1onf.jld2"), D)
end
