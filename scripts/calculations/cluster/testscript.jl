#! /bin/bash
#=
exec julia -t auto "${BASH_SOURCE[0]}" "$@"
=#
using DrWatson
@quickactivate "SpatiotemporalMotifs"
using Distributed
import USydClusters.Physics: addprocs

project = "/headnode2/bhar9988/code/DDC/SpatiotemporalMotifs/"
procs = addprocs(5; ncpus = 2, mem = 3,
                 walltime = 96, project)
@everywhere import SpatiotemporalMotifs: test_calculations
pmap(test_calculations, p)
display("All workers completed")
@sync USydClusters.Physics.selfdestruct()
