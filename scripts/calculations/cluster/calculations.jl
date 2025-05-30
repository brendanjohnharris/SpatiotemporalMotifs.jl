#! /bin/bash
#=
exec julia +1.10.9 -t auto "${BASH_SOURCE[0]}" "$@"
=#
using DrWatson
@quickactivate "SpatiotemporalMotifs"
project = Base.active_project()
import AllenNeuropixels as AN
import SpatiotemporalMotifs as SM
import USydClusters.Physics: addprocs, selfdestruct
using USydClusters
SM.@preamble

session_table = load(datadir("session_table.jld2"), "session_table")
oursessions = session_table.ecephys_session_id

outpath = datadir("calculations")
rewrite = false

if haskey(ENV, "JULIA_DISTRIBUTED") # ? Should take a night or so
    exprs = map(oursessions) do sessionid
        expr = quote
            using Pkg
            Pkg.instantiate()
            import SpatiotemporalMotifs: send_calculations
            send_calculations($sessionid; outpath = $outpath, rewrite = $rewrite)
        end
    end
    USydClusters.Physics.runscripts(exprs; ncpus = 8, mem = 60, walltime = 8,
                                    project = projectdir(), exeflags = `+1.10.9`,
                                    queue = "taiji")

    display("All workers submitted")
else
    SM.send_calculations.(reverse(oursessions); outpath, rewrite) # ? This version will take a few days if the above calculations errored, otherwise a few minutes (checks all calculations are correct)
    Q = SM.calcquality(datadir("calculations"))
    @assert all(oursessions .∈ [lookup(Q, 3)])
    @assert mean(Q) == 1
end
# @sync selfdestruct()
