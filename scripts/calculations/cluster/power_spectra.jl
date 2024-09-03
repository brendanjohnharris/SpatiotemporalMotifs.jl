#! /bin/bash
#=
exec julia -t auto --heap-size-hint=`grep MemFree /proc/meminfo | awk '{print int($2 * 0.4) "k"}'` "${BASH_SOURCE[0]}" "$@"
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
    exprs = map(oursessions) do o
        expr = quote
            using Pkg
            Pkg.instantiate()
            import SpatiotemporalMotifs as SM
            SM.send_powerspectra($o; rewrite = false, retry_errors = true)
        end
    end
    USydClusters.Physics.runscript.(exprs; ncpus = 16, mem = 90, walltime = 4,
                                    project = projectdir())
else
    for o in reverse(oursessions)
        SM.send_powerspectra(o; rewrite = false, retry_errors = true)
        GC.gc()
    end
end
