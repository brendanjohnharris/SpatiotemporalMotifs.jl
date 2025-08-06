#! /bin/bash
#=
exec julia +1.10.10 -t auto "${BASH_SOURCE[0]}" "$@"
=#
using DrWatson
@quickactivate "SpatiotemporalMotifs"

import TimeseriesTools: freqs
using FileIO
using Unitful
import DimensionalData: metadata
using MultivariateStats
using SpatiotemporalMotifs
@preamble
set_theme!(foresight(:physics))
import SpatiotemporalMotifs: prony, ramap

begin
    config = Dict("p" => 100,
                  "band" => (1, 30),
                  "interval" => (0.25, 1))
end

if !isfile(calcdir("plots", savepath("figS2", config, "jld2")))
    if nprocs() == 1
        if SpatiotemporalMotifs.CLUSTER()
            using USydClusters
            ourprocs = USydClusters.Physics.addprocs(6; mem = 33, ncpus = 5,
                                                     project = projectdir(),
                                                     queue = "taiji")
        else
            addprocs(6) # Equal to the number of structures
        end
    end
    @everywhere using SpatiotemporalMotifs
    @everywhere SpatiotemporalMotifs.@preamble
end

plot_data, data_file = produce_or_load(config, calcdir("plots");
                                       filename = savepath("figS2")) do config
    outdict = Dict{String, Dict{String, Any}}()

    interval = Interval((config["interval"] .* [u"s"])...) # Offset onwards
    p = config["p"]
    pronyband = config["band"]

    function pronyextract(P, x)
        freqs = frequencies(P)
        amps = amplitudes(P)

        # * Find strongest positive frequency
        idxs = freqs .> 0
        _, idx = findmax(amps[idxs])
        topfreq = freqs[idxs][idx]

        # * Reconstruct the signal from this component and its reflection
        topidx = findfirst(isequal(topfreq), freqs)
        idxs = [topidx, findfirst(isequal(-topfreq), freqs)]

        y = reconstruct(P, x; idxs)
        rmse = sqrt(mean(abs2, y - x))

        return amps[topidx], freqs[topidx], rmse
    end

    begin
        stimulus = r"Natural_Images"
        begin # * Load the trial LFP for natural images
            session_table = load(calcdir("posthoc_session_table.jld2"), "session_table")
            oursessions = session_table.ecephys_session_id
            path = calcdir("calculations")
            Q = calcquality(path)[Structure = At(structures)]
            Q = Q[SessionID(At(oursessions))]
            @assert mean(Q[stimulus = At(stimulus)]) == 1
            out = load_calculations(Q; stimulus = stimulus, vars = [:V])
        end
        begin # * Calculate prony fits for each channel. Takes about 30 minutes.
            @info "Calculating prony fits for $stimulus"
            pronies = pmap(out) do out_structure
                sessionids = getindex.(out_structure, :sessionid)

                P = map(out_structure) do out_session
                    V = out_session[:V]
                    V = set(V, 𝑡 => Float32.(times(V)))
                    V = V[𝑡 = interval] |> ustripall
                    V = ZScore(V, dims = 1)(V)
                    V = bandpass(V, pronyband)
                    P = map(Base.Fix2(prony, p), eachslice(V, dims = (Depth, :changetime)))
                end
                pronies = ToolsArray(P, (SessionID(sessionids),))
            end
        end

        ampfreqfits = ramap(pronyextract, pronies)

        outdict[string(stimulus)] = Dict("prony_fit" => pronies)
    end
    uni = []
    GC.gc()

    begin # * Load the evoked oscillations for the flashes stimulus
        stimulus = "flash_250ms"
        begin # * Load the trial LFP for natural images
            session_table = load(calcdir("posthoc_session_table.jld2"), "session_table")
            oursessions = session_table.ecephys_session_id
            path = calcdir("calculations")
            Q = calcquality(path)[Structure = At(structures)]
            Q = Q[SessionID(At(oursessions))]
            @assert mean(Q[stimulus = At(stimulus)]) == 1
            out = load_calculations(Q; stimulus = stimulus, vars = [:V])
        end

        begin # * Calculate prony fits for each channel. Takes about 30 minutes.
            @info "Calculating prony fits for $stimulus"
            pronies = pmap(out) do out_structure
                sessionids = getindex.(out_structure, :sessionid)

                P = map(out_structure) do out_session
                    V = out_session[:V]
                    V = V[𝑡 = interval] |> ustripall
                    V = ZScore(V, dims = 1)(V)
                    V = bandpass(V, pronyband)
                    P = map(Base.Fix2(prony, p), eachslice(V, dims = (Depth, :changetime)))
                end
                pronies = ToolsArray(P, (SessionID(sessionids),))
            end
        end
        outdict[string(stimulus)] = Dict("prony_fit" => pronies)
    end
    out = []
    GC.gc()
    return outdict
end

begin # * Select the peak amplitude and frequency for each channel
    pronies = plot_data["flash_250ms"]["prony_fit"]
    af = ramap(pronies) do x
        frequencies = x.frequencies
        amplitudes = x.amplitudes

        idxs = frequencies .> 0
        frequencies = frequencies[idxs]
        amplitudes = amplitudes[idxs]

        idx = findmax(amplitudes) |> last

        return amplitudes[idx], frequencies[idx]
    end
    amplitudes = ramap(first, af)
    frequencies = ramap(last, af)
end
# begin # * Plot normalized amplitude distribution
# end

# begin
#     f = TwoPanel()
#     gs = subdivide(f, 1, 3)
#     layerints = load(calcdir("plots", "grand_unified_layers.jld2"), "layerints")

#     titles = ["Natural images post-onset", "Flashes post-onset"]

#     for (i, stimulus) in enumerate(["r\"Natural_Images\"", "flash_250ms"])
#         ax = Axis(gs[i], xlabel = "Average frequency (Hz)", ylabel = "Cortical depth (%)",
#                   title = "$(titles[i])",
#                   ytickformat = depthticks, yreversed = true,
#                   limits = ((4.7, 7), (0, 1)))
#         ω = plot_data[stimulus]["ω"]

#         plotlayerints!(ax, layerints; axis = :y, flipside = false, newticks = false,
#                        bgcolor = Makie.RGBA(0, 0, 0, 0))

#         map(lookup(ω, Structure), eachslice(ω, dims = Structure)) do structure, x
#             x = ustripall(x ./ (2π))

#             μ, (σl, σh) = bootstrapmedian(x, dims = 2)

#             μ, σl, σh = upsample.((μ, σl, σh), 5)

#             band!(ax, Point2f.(collect(σl), lookup(μ, 1)),
#                   Point2f.(collect(σh), lookup(μ, 1));
#                   color = (structurecolormap[structure], 0.32), label = structure)

#             lines!(ax, collect(μ), lookup(μ, 1);
#                    color = (structurecolormap[structure], 0.8),
#                    label = structure)
#         end
#     end
#     Legend(gs[3], contents(gs[2]) |> only; merge = true, title = "Structure")
#     addlabels!(f, labelformat)
#     display(f)
#     wsave(plotdir("figS2", "figS2.pdf"), f)
# end
