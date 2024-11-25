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

isdir(datadir("power_spectra")) || mkpath(datadir("power_spectra"))

rewrite = false
retry_errors = true

structures = deepcopy(SM.structures)
# Add thalamic regions
push!(structures, "LGd")
push!(structures, "LGd-sh")
push!(structures, "LGd-co")
params = Iterators.product(oursessions, stimuli, unique(structures)) |> collect
idxs = map(xy -> SM.powerspectra_quality(xy...; rewrite,
                                         retry_errors), params)
params = params[.!idxs]

if haskey(ENV, "JULIA_DISTRIBUTED")
    exprs = map(params) do (o, stimulus, structure)
        expr = quote
            using Pkg
            Pkg.instantiate()
            import SpatiotemporalMotifs as SM
            SM.send_powerspectra($o, $stimulus, $structure; rewrite = $rewrite,
                                 retry_errors = $retry_errors)
        end
    end
    nn = length(exprs) ÷ 2
    USydClusters.Physics.runscripts(exprs[1:nn]; ncpus = 20, mem = 120, walltime = 1,
                                    project = projectdir(), qsub_flags = "-q yossarian")
    USydClusters.Physics.runscripts(exprs[(nn + 1):end]; ncpus = 20, mem = 120,
                                    walltime = 1,
                                    project = projectdir(), qsub_flags = "-q yossarian")
else
    for param in params
        SM.send_powerspectra(param...; rewrite, retry_errors)
    end
end
