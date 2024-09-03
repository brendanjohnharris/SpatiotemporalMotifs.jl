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

stimuli = ["spontaneous", "flash_250ms", r"Natural_Images"]
session_table = load(datadir("session_table.jld2"), "session_table")
oursessions = session_table.ecephys_session_id

if haskey(ENV, "JULIA_DISTRIBUTED")
    params = Iterators.product(oursessions, stimuli) |> collect
    idxs = map(xy -> SM.powerspectra_quality(xy...; rewrite = false,
                                             retry_errors = true), params)
    params = params[.!idxs]
    exprs = map(params) do (o, stimulus)
        expr = quote
            using Pkg
            Pkg.instantiate()
            import SpatiotemporalMotifs as SM
            SM.send_powerspectra($o, $stimulus; rewrite = false, retry_errors = true)
        end
    end
    USydClusters.Physics.runscripts(exprs; ncpus = 14, mem = 60, walltime = 2,
                                    project = projectdir())
else
    for o in reverse(oursessions)
        for stimulus in stimuli
            SM.send_powerspectra(o, stimulus; rewrite = false, retry_errors = true)
            GC.gc()
        end
    end
end
