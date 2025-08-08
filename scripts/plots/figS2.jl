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

begin
    config = Dict("p" => 50,
                  "band" => (1, 30),
                  "interval" => (0.25, 1),
                  "lambda" => 1e-8)
end

if !isfile(calcdir("plots", savepath("figS2", config, "jld2")))
    if nprocs() == 1
        if SpatiotemporalMotifs.CLUSTER()
            using USydClusters
            ourprocs = USydClusters.Physics.addprocs(16; mem = 12, ncpus = 2,
                                                     project = projectdir(),
                                                     queue = "taiji")
        else
            addprocs(9)
        end
    end
    @everywhere using SpatiotemporalMotifs
    @everywhere SpatiotemporalMotifs.@preamble
end

plot_data, data_file = produce_or_load(config, calcdir("plots");
                                       filename = savepath("figS2")) do config
    outdict = Dict{String, Dict{String, Any}}()

    interval = Interval(config["interval"]...) # Offset onwards
    p = config["p"]
    λ = config["lambda"]
    pronyband = config["band"]

    function prony_fit(out)
        @withprogress name="Prony" begin
            threadmax = length(out)
            threadlog = 0
            pronies = map(out) do out_structure
                sessionids = getindex.(out_structure, :sessionid)

                P = map(out_structure) do out_session
                    V = out_session[:V] |> ustripall
                    V = bandpass(V, pronyband) # Reduce high-frequency noise
                    V = V[𝑡 = interval]
                    V = V[1:2:end, :, :] # Downsample for speed
                    V = ZScore(V, dims = 1)(V)

                    P = map(eachslice(V, dims = (Depth, :changetime))) do x
                        # * Fit prony
                        P = SpatiotemporalMotifs.prony(x, p; λ = λ) # Real prony fit
                        fs = SpatiotemporalMotifs.frequencies(P)
                        amps = SpatiotemporalMotifs.amplitudes(P)
                        rmses = SpatiotemporalMotifs.errors(P, x)
                        damps = SpatiotemporalMotifs.dampingfactors(P)

                        # * Find strongest frequency
                        rmse, topidx = findmin(rmses)

                        # begin
                        #     lines(x)
                        #     y = SpatiotemporalMotifs.reconstruct(P, x)
                        #     lines!(y, linestyle = :dash)
                        #     y = SpatiotemporalMotifs.reconstruct(P, x; idxs = topidx)
                        #     lines!(y)
                        #     current_figure() |> display
                        # end

                        return amps[topidx], fs[topidx], rmse, damps[topidx]
                    end
                end

                threadlog += 1
                @logprogress threadlog / threadmax
                return ToolsArray(P, (SessionID(sessionids),))
            end
        end
        return pronies
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
        @info "Calculating prony fits for $stimulus"
        pronies = prony_fit(out)

        amplitudes = SpatiotemporalMotifs.ramap(Base.Fix2(getindex, 1), pronies)
        frequencies = SpatiotemporalMotifs.ramap(Base.Fix2(getindex, 2), pronies)
        rmse = SpatiotemporalMotifs.ramap(Base.Fix2(getindex, 3), pronies)
        damps = SpatiotemporalMotifs.ramap(Base.Fix2(getindex, 4), pronies)

        outdict[string(stimulus)] = Dict("prony_fit" => Dict("amplitudes" => amplitudes,
                                                             "frequencies" => frequencies,
                                                             "rmse" => rmse,
                                                             "damps" => damps))
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

        @info "Calculating prony fits for $stimulus"
        pronies = prony_fit(out)

        amplitudes = SpatiotemporalMotifs.ramap(Base.Fix2(getindex, 1), pronies)
        frequencies = SpatiotemporalMotifs.ramap(Base.Fix2(getindex, 2), pronies)
        rmse = SpatiotemporalMotifs.ramap(Base.Fix2(getindex, 3), pronies)
        damps = SpatiotemporalMotifs.ramap(Base.Fix2(getindex, 4), pronies)

        outdict[string(stimulus)] = Dict("prony_fit" => Dict("amplitudes" => amplitudes,
                                                             "frequencies" => frequencies,
                                                             "rmse" => rmse,
                                                             "damps" => damps))
    end
    out = []
    GC.gc()
    return outdict
end

# begin # * Select the peak amplitude and frequency for each channel
#     pronies = plot_data["flash_250ms"]["prony_fit"]
#     af = ramap(pronies) do x
#         frequencies = x.frequencies
#         amplitudes = x.amplitudes

#         idxs = frequencies .> 0
#         frequencies = frequencies[idxs]
#         amplitudes = amplitudes[idxs]

#         idx = findmax(amplitudes) |> last

#         return amplitudes[idx], frequencies[idx]
#     end
#     amplitudes = ramap(first, af)
#     frequencies = ramap(last, af)
# end
function pformat(x)
    bins = intervals(0.0:0.1:0.9)
    x = mean(x, dims = :changetime)
    x = dropdims(x, dims = (𝑡, :changetime))
    x = groupby(x, Depth => Bins(bins))
    x = mean.(x)
    x = set(x, Depth => mean.(lookup(x, Depth)))
end
begin # * Plot normalized amplitude distribution (flashes)
    stimulus = "flash_250ms"

    amps = plot_data[stimulus]["prony_fit"]["amplitudes"]
    bamps = SpatiotemporalMotifs.ramap(pformat, ToolsArray{T, 2} where {T}, amps) .|> stack

    fs = plot_data[stimulus]["prony_fit"]["frequencies"]
    bfs = SpatiotemporalMotifs.ramap(pformat, ToolsArray{T, 2} where {T}, fs) .|> stack

    errors = plot_data[stimulus]["prony_fit"]["rmse"]
    berrors = SpatiotemporalMotifs.ramap(pformat, ToolsArray{T, 2} where {T}, errors) .|>
              stack
end
begin
    f = Figure()
    ax = Axis(f[1, 1])
    map(structures, bamps) do structure, x
        # m = mean(amp, dims = 2)
        # m = dropdims(m, dims = 2)

        μ, (σl, σh) = bootstrapmedian(x, dims = 2)

        band!(ax, Point2f.(collect(σl), lookup(μ, 1)),
              Point2f.(collect(σh), lookup(μ, 1));
              color = (structurecolormap[structure], 0.32), label = structure)
        lines!(ax, parent(μ), lookup(μ, 1); color = structurecolormap[structure])
    end

    f
end

amps = map(structures, amps) do structure, amp
    # * Group across unified depths

end

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
