#! /bin/bash
#=
exec julia +1.10.9 -t auto "${BASH_SOURCE[0]}" "$@"
=#
using DrWatson
@quickactivate "SpatiotemporalMotifs"
project = Base.active_project()
import AllenNeuropixels as AN
using SpatiotemporalMotifs
import USydClusters.Physics: addprocs, selfdestruct
@preamble

session_table = load(datadir("session_table.jld2"), "session_table")
thalamic_structures = ["LGd-co", "LGd-sh"]
session_table = subset(session_table,
                       :ecephys_structure_acronyms => ByRow(x -> any(thalamic_structures .∈
                                                                     [x])))
oursessions = session_table.ecephys_session_id

outpath = datadir("thalamus_calculations")
rewrite = false

if haskey(ENV, "JULIA_DISTRIBUTED")
    procs = addprocs(10; ncpus = 8, mem = 45,
                     walltime = 96, project) # ! If you have workers dying unexpectedly, try increasing the memory for each job
    @everywhere begin
        using SpatiotemporalMotifs
        import SpatiotemporalMotifs.on_error
        thalamic_structures = $thalamic_structures
        outpath = $outpath
        rewrite = $rewrite
    end
    O = pmap(x -> send_thalamus_calculations.(x; structures = ["LGd"],
                                              rewrite, outpath),
             oursessions; on_error)
    fetch.(O)
    display("All workers completed")
end
send_thalamus_calculations.(oursessions; rewrite, outpath,
                            structures = ["LGd"])
