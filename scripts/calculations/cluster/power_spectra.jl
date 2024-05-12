#! /bin/bash
#=
exec julia -t auto "${BASH_SOURCE[0]}" "$@"
=#
using DrWatson
@quickactivate "SpatiotemporalMotifs"
import SpatiotemporalMotifs as SM
using USydClusters
SM.@preamble
set_theme!(foresight(:physics))
ENV["JULIA_DEBUG"] = "AllenNeuropixelsBase"

session_table = load(datadir("session_table.jld2"), "session_table")
oursessions = session_table.ecephys_session_id

if haskey(ENV, "JULIA_DISTRIBUTED")
    proclist = USydClusters.Physics.addprocs(4; ncpus = 2, mem = 125,
                                             walltime = 96, project = projectdir())
    @everywhere import SpatiotemporalMotifs as SM
    O = [remotecall(SM.send_powerspectra, p, s; rewrite = false)
         for (p, s) in zip(proclist, oursessions)]
    wait.(O) # Wait for workers to finish
    display("All workers completed")
    SM.send_powerspectra.(oursessions; rewrite = false) #  A final local check, in case a few workers failed.
# @sync USydClusters.Physics.selfdestruct()
else
    SM.send_powerspectra.(oursessions; rewrite = false)
end
