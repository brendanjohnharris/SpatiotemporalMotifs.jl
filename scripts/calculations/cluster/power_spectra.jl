#! /bin/bash
#=
exec julia +1.10.10 -t auto "${BASH_SOURCE[0]}" "$@"
=#
using DrWatson
@quickactivate "SpatiotemporalMotifs"
import SpatiotemporalMotifs as SM
using USydClusters
SM.@preamble
set_theme!(foresight(:physics))
ENV["JULIA_DEBUG"] = "AllenNeuropixelsBase"

path = calcdir("power_spectra")
mkpath(path)
stimuli = ["spontaneous", "flash_250ms", r"Natural_Images"]
session_table = load(calcdir("session_table.jld2"), "session_table")
oursessions = session_table.ecephys_session_id

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

Q = SM.calcquality(path)

if haskey(ENV, "SM_CLUSTER") && !isempty(params)
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
elseif ENV["HOSTNAME"] == "cartman.physics.usyd.edu.au"
    addprocs(3)
    @everywhere import SpatiotemporalMotifs as SM
    pmap(_params) do param
        SM.send_powerspectra(param...; rewrite, retry_errors)
        GC.gc()
    end
else
    @progress for param in _params
        SM.send_powerspectra(param...; rewrite, retry_errors)
        GC.gc()
    end
end

begin
    dir = "/import/taiji1/bhar9988/code/DDC/SpatiotemporalMotifs/data/power_spectra"
    files = readdir(dir)

    map(files) do file
        _, params, _ = parse_savename(file; connector = SpatiotemporalMotifs.connector)
        plotfile = joinpath(calcdir("plots", "power_spectra"), "$(params["sessionid"])",
                            "$(params["stimulus"])_$(params["structure"]).pdf")
        plotfile = relpath(plotfile, projectdir())

        pac_plotfile = joinpath(calcdir("plots", "power_spectra"),
                                "$(params["sessionid"])",
                                "$(params["stimulus"])_$(params["structure"])_pac.pdf")
        pac_plotfile = relpath(pac_plotfile, projectdir())

        jldopen(joinpath(dir, file), "a+") do f
            f["plotfiles"] = [plotfile, pac_plotfile]
        end
    end
end
