#! /bin/bash
#=
exec julia +1.10.10 -t auto --heap-size-hint=`grep MemFree /proc/meminfo | awk '{print int($2 * 0.4) "k"}'` "${BASH_SOURCE[0]}" "$@"
=#
using DrWatson
@quickactivate "SpatiotemporalMotifs"
import SpatiotemporalMotifs as SM
using USydClusters
SM.@preamble
set_theme!(foresight(:physics))
ENV["JULIA_DEBUG"] = "AllenNeuropixelsBase"

stimuli = ["spontaneous", "flash_250ms", r"Natural_Images"]
session_table = load(calcdir("session_table.jld2"), "session_table")
oursessions = session_table.ecephys_session_id

isdir(calcdir("power_spectra")) || mkpath(calcdir("power_spectra"))

rewrite = false
retry_errors = true

pstructures = deepcopy(SM.structures)
# Add thalamic regions
push!(pstructures, "LGd")
push!(pstructures, "LGd-sh")
push!(pstructures, "LGd-co")
_params = Iterators.product(oursessions, stimuli, unique(pstructures)) |> collect
idxs = map(xy -> SM.powerspectra_quality(xy...; rewrite,
                                         retry_errors), _params)
params = _params[.!idxs]

if haskey(ENV, "JULIA_DISTRIBUTED") && !isempty(params)
    exprs = map(params) do (o, stimulus, structure)
        expr = quote
            using Pkg
            Pkg.instantiate()
            import SpatiotemporalMotifs as SM
            SM.send_powerspectra($o, $stimulus, $structure; rewrite = $rewrite,
                                 retry_errors = $retry_errors)
        end
    end
    SM.submit_calculations(exprs, mem = 100)
else
    for param in _params
        SM.send_powerspectra(param...; rewrite, retry_errors)
    end
end
