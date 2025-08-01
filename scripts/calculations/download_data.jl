#! /bin/bash
#=
exec julia +1.10.10 -t 1 "${BASH_SOURCE[0]}" "$@"
=#
using DrWatson
@quickactivate "SpatiotemporalMotifs"
import SpatiotemporalMotifs as SM

session_table = load(calcdir("session_table.jld2"), "session_table")
oursessions = session_table.ecephys_session_id

for s in oursessions # * Will take about 10 minutes if files are already downloaded
    @info "Loading session $s"
    session = AN.Session(s)
    probes = AN.getprobes(session)
    for probeid in probes.id
        @info "Loading probe $probeid" # Load a piece of data to ensure the files have downloaded
        LFP = AN.getlfp(session, probeid; times = 100 .. 101)
        LFP = []
        GC.gc()
    end
end
