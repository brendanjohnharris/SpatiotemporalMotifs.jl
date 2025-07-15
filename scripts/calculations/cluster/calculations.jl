#! /bin/bash
#=
exec julia +1.10.10 -t auto "${BASH_SOURCE[0]}" "$@"
=#
using DrWatson
@quickactivate "SpatiotemporalMotifs"
project = Base.active_project()
import AllenNeuropixels as AN
import SpatiotemporalMotifs as SM
import USydClusters.Physics: addprocs, selfdestruct
using USydClusters
SM.@preamble

session_table = load(calcdir("session_table.jld2"), "session_table")
oursessions = session_table.ecephys_session_id

outpath = calcdir("calculations")
rewrite = false
check_quality = true
mkpath(outpath)

if haskey(ENV, "JULIA_DISTRIBUTED") # ? Should take a night or so
    exprs = map(oursessions) do sessionid
        expr = quote
            using Pkg
            Pkg.instantiate()
            import SpatiotemporalMotifs: send_calculations
            send_calculations($sessionid; outpath = $outpath, rewrite = $rewrite)
        end
    end

    if check_quality && isfile(calcdir("posthoc_session_table.jld2")) &&
       !isempty(readdir(outpath))
        Q = SM.calcquality(outpath)[Structure = At(SM.structures)]
        # * Delete bad files
        filenames = collect(Iterators.product(lookup(Q)...))[.!Q]
        filenames = map(filenames) do f
            Dict{String, Any}("structure" => f[1], "stimulus" => f[2], "sessionid" => f[3])
        end
        filenames = calcdir.(["calculations"], map(SM.savepath, filenames) .* [".jld2"])
        rm.(filenames; force = true)

        Q = any(.!Q, dims = (SM.Structure, Dim{:stimulus})) # Sessions that have bad files
        Q = dropdims(Q, dims = (SM.Structure, Dim{:stimulus}))
        donesessions = lookup(Q, :SessionID)[findall(Q)]
        exprs = exprs[.!(oursessions .∈ [donesessions])]
    end

    N3 = length(exprs) ÷ 3
    USydClusters.Physics.runscripts(exprs[1:N3]; ncpus = 8, mem = 42, walltime = 8,
                                    project = projectdir(), exeflags = `+1.10.10`,
                                    queue = "l40s")
    USydClusters.Physics.runscripts(exprs[(N3 + 1):(2 * N3)]; ncpus = 8, mem = 42,
                                    walltime = 8,
                                    project = projectdir(), exeflags = `+1.10.10`,
                                    queue = "h100")
    USydClusters.Physics.runscripts(exprs[(2 * N3 + 1):end]; ncpus = 8, mem = 42,
                                    walltime = 8,
                                    project = projectdir(), exeflags = `+1.10.10`,
                                    qsub_flags = `-l node=cmt01`)

    display("All workers submitted")
else
    SM.send_calculations.(reverse(oursessions); outpath, rewrite) # ? This version will take a few days if the above calculations errored, otherwise a few minutes (checks all calculations are correct)
    Q = SM.calcquality(calcdir("calculations"))
    @assert all(oursessions .∈ [lookup(Q, 3)])
    @assert mean(Q) == 1
end
# @sync selfdestruct()
