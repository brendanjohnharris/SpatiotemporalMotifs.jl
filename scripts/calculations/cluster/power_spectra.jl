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

path = calcdir("power_spectra")
mkpath(path)
stimuli = ["spontaneous", "flash_250ms", r"Natural_Images"]
session_table = load(calcdir("session_table.jld2"), "session_table")
oursessions = session_table.ecephys_session_id

pstructures = deepcopy(SM.structures)
# Add thalamic regions
push!(pstructures, "LGd")
push!(pstructures, "LGd-sh")
push!(pstructures, "LGd-co")
_params = Iterators.product(oursessions, stimuli, unique(pstructures)) |> collect

Q = SM.calcquality(path)
params = []
map(_params) do p
    sessionid, stimulus, structure = p
    try
        if !Q[stimulus = At(stimulus),
              Structure = At(structure),
              SessionID = At(sessionid)]
            push!(params, p)
        end
    catch
        push!(params, p)
    end
end

if !isempty(params)
    if SM.CLUSTER()
        exprs = map(params) do (o, stimulus, structure)
            expr = quote
                using Pkg
                Pkg.instantiate()
                import SpatiotemporalMotifs as SM
                SM.send_powerspectra($o, $stimulus, $structure)
            end
        end
        SM.submit_calculations(exprs, mem = 49, ncpus = 8, walltime = 8)
    elseif ENV["HOSTNAME"] ∈ ["cartman.physics.usyd.edu.au", "stan.physics.usyd.edu.au"]
        addprocs(7)
        @everywhere import SpatiotemporalMotifs as SM
        @everywhere using Suppressor

        channel = RemoteChannel(() -> Channel{Bool}(length(_params)), 1)
        threadlog = Threads.Atomic{Int}(0)
        threadmax = length(_params)
        lk = Threads.ReentrantLock()
        @withprogress name="Power spectra" begin
            @async while take!(channel)
                lock(lk) do
                    Threads.atomic_add!(threadlog, 1)
                    @logprogress threadlog[] / threadmax
                end
            end

            @sync begin
                pmap(params) do param
                    @info "Calculating power spectra for $(param)"

                    @suppress SM.send_powerspectra(param...)
                    GC.gc()

                    put!(channel, true)
                    yield()
                end
                put!(channel, false)
            end
        end
    else
        @progress for param in params
            SM.send_powerspectra(param...)
            GC.gc()
        end
    end
end
