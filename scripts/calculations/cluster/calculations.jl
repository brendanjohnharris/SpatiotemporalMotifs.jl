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

session_table = load(datadir("session_table.jld2"), "session_table")
oursessions = session_table.ecephys_session_id

outpath = datadir("calculations")
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

    if check_quality && !isempty(readdir(outpath))
        Q = SM.calcquality(outpath)[Structure = At(SM.structures)]
        # * Delete bad files
        filenames = collect(Iterators.product(lookup(Q)...))[.!Q]
        filenames = map(filenames) do f
            Dict{String, Any}("structure" => f[1], "stimulus" => f[2], "sessionid" => f[3])
        end
        filenames = calcdir.(["calculations"], map(SM.savepath, filenames) .* [".jld2"])
        # rm.(filenames; force = true)

        Q = any(.!Q, dims = (SM.Structure, Dim{:stimulus})) # Sessions that have bad files
        Q = dropdims(Q, dims = (SM.Structure, Dim{:stimulus}))
        notgood = lookup(Q, :SessionID)[findall(Q)] |> collect
        Base.Fix1(push!, notgood).(setdiff(oursessions, lookup(Q, 1)))
        exprs = exprs[(oursessions .∈ [notgood])]
    end

    SM.submit_calculations(exprs)

    display("All workers submitted")
else
    SM.send_calculations.(reverse(oursessions); outpath, rewrite) # ? This version will take a few days if the above calculations errored, otherwise a few minutes (checks all calculations are correct)
    Q = SM.calcquality(datadir("calculations"))
    @assert all(oursessions .∈ [lookup(Q, 3)])
    @assert mean(Q) == 1
end
# @sync selfdestruct()
