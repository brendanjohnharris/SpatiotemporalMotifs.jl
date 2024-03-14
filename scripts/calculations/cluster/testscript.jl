#! /bin/bash
#=
exec julia -t auto "${BASH_SOURCE[0]}" "$@"
=#
using DrWatson
@quickactivate "SpatiotemporalMotifs"
import SpatiotemporalMotifs as SM
using Distributed
using USydClusters

project = "/headnode2/bhar9988/code/DDC/SpatiotemporalMotifs/"
procs = USydClusters.Physics.addprocs(5; ncpus = 2, mem = 15,
                                      walltime = 96, project)
@everywhere import SpatiotemporalMotifs as SM
function send()
    [remotecall(SM.test_calculations, p, nothing) for p in procs]
end
O = send()
wait.(O) # Wait for workers to finish
fetch.(O)
display("All workers completed")
@sync USydClusters.Physics.selfdestruct()
