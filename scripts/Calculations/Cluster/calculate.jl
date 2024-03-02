@info @__FILE__
#! /bin/bash
#=
exec julia -t auto "${BASH_SOURCE[0]}" "$@"
=#
using DrWatson
@quickactivate "SpatiotemporalMotifs"
using Unitful
import SpatiotemporalMotifs as SM
using USydClusters
SM.@preamble
set_theme!(foresight(:dark, :serif, :physics))

session_table = load(datadir("session_table.jld2"), "session_table")
oursessions = session_table.ecephys_session_id

preamble = "import SpatiotemporalMotifs as SM"
O = USydClusters.Physics.addprocs(oursessions; preamble, ncpus = 4, mem = 31,
                                  walltime = 96) do s
    SM.send_calculations(s; outpath = datadir("Calculations"), rewrite = false)
end
wait.(O)
fetch.(O) # Wait for workers to finish
display("All workers completed")
exit()
