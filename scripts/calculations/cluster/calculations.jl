#! /bin/bash
#=
exec julia -t auto "${BASH_SOURCE[0]}" "$@"
=#
using DrWatson
@quickactivate "SpatiotemporalMotifs"
project = Base.active_project()
import SpatiotemporalMotifs as SM
import USydClusters.Physics: addprocs, selfdestruct
SM.@preamble
set_theme!(foresight(:dark, :serif, :physics))

session_table = load(datadir("session_table.jld2"), "session_table")
oursessions = session_table.ecephys_session_id

outpath = datadir("calculations")
rewrite = false
distribute = true

if distribute
    procs = addprocs(10; ncpus = 8, mem = 31,
                     walltime = 96, project) # ! If you have workers dying unexpectedly, try increasing the memory for each job
    @everywhere begin
        import SpatiotemporalMotifs: send_calculations, on_error
        outpath = $outpath
        rewrite = $rewrite
    end
    O = pmap(x -> send_calculations.(x; outpath, rewrite), oursessions; on_error)
    fetch.(O)
    display("All workers completed")
    @sync selfdestruct()
else
    send_calculations.(oursessions; outpath, rewrite)
end
