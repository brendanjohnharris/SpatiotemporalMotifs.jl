#! /bin/bash
#=
exec julia +1.10.10 -t auto "${BASH_SOURCE[0]}" "$@"
=#
using DrWatson
@quickactivate "SpatiotemporalMotifs"
using Distributed
import USydClusters.Physics: addprocs, selfdestruct

project = "/taiji1/bhar9988/code/DDC/SpatiotemporalMotifs/"
procs = addprocs(5; ncpus = 2, mem = 4,
                 walltime = 96, project)
@everywhere using Pkg;
@everywhere Pkg.instantiate();
@everywhere import SpatiotemporalMotifs: test_calculations
pmap(test_calculations, procs)
display("All workers completed")
@sync selfdestruct()
