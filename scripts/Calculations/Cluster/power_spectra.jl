#! /bin/bash
#=
exec julia -t 1 "${BASH_SOURCE[0]}" "$@"
=#
using DrWatson
@quickactivate "SpatiotemporalMotifs"

import SpatiotemporalMotifs as SM
using FileIO
SM.@preamble
set_theme!(foresight(:physics))

session_table = load(datadir("session_table.jld2"), "session_table")
oursessions = session_table.ecephys_session_id

preamble = "import SpatiotemporalMotifs as SM"
O = USydClusters.Physics.addprocs(oursessions; preamble, ncpus = 2, mem = 15,
                                  walltime = 96) do s
    SM.send_powerspectra(s; rewrite = false)
end
wait.(O)
fetch.(O) # Wait for workers to finish
display("All workers completed")
USydClusters.Physics.selfdestruct()
