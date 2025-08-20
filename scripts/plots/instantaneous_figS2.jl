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

plot_data, data_file = produce_or_load(Dict(), calcdir("plots");
                                       filename = savepath("figS2")) do _
    outdict = Dict{String, Dict{String, Any}}()

    interval = 0.3u"s" .. 0.5u"s"
    begin
        stimulus = r"Natural_Images"
        begin # * Load the trial LFP for natural images
            session_table = load(calcdir("plots", "posthoc_session_table.jld2"),
                                 "session_table")
            oursessions = session_table.ecephys_session_id
            path = calcdir("calculations")
            Q = calcquality(path)[Structure = At(structures)]
            Q = Q[SessionID(At(oursessions))]
            @assert mean(Q[stimulus = At(stimulus)]) == 1
            out = load_calculations(Q; stimulus = stimulus, vars = [:ω])
        end

        begin
            bins = intervals(0.0:0.1:0.9)
            ω = map(out) do out_structure
                sessionids = getindex.(out_structure, :sessionid)
                x = map(out_structure) do out_session
                    ω = out_session[:ω]
                    ω = ω[𝑡 = interval]
                    ω = mean(ω, dims = (𝑡, :changetime))
                    ω = dropdims(ω, dims = (𝑡, :changetime))
                    ω = groupby(ω, Depth => Bins(bins))
                    ω = mean.(ω)
                    ω = set(ω, Depth => mean.(lookup(ω, Depth)))
                end
                ToolsArray(x, (SessionID(sessionids),)) |> stack
            end
            ω = ToolsArray(ω, (Structure(structures),)) |> stack
        end
        outdict[string(stimulus)] = Dict("ω" => ω)
    end
    uni = []
    GC.gc()

    begin # * Load the evoked oscillations for the flashes stimulus
        stimulus = "flash_250ms"
        begin # * Load the trial LFP for natural images
            session_table = load(calcdir("plots", "posthoc_session_table.jld2"),
                                 "session_table")
            oursessions = session_table.ecephys_session_id
            path = calcdir("calculations")
            Q = calcquality(path)[Structure = At(structures)]
            Q = Q[SessionID(At(oursessions))]
            @assert mean(Q[stimulus = At(stimulus)]) == 1
            out = load_calculations(Q; stimulus = stimulus, vars = [:ω])
        end

        begin
            bins = intervals(0.0:0.1:0.9)
            ω = map(out) do out_structure
                sessionids = getindex.(out_structure, :sessionid)
                x = map(out_structure) do out_session
                    ω = out_session[:ω]
                    ω = ω[𝑡 = interval]
                    ω = mean(ω, dims = (𝑡, :changetime))
                    ω = dropdims(ω, dims = (𝑡, :changetime))
                    ω = groupby(ω, Depth => Bins(bins))
                    ω = mean.(ω)
                    ω = set(ω, Depth => mean.(lookup(ω, Depth)))
                end
                ToolsArray(x, (SessionID(sessionids),)) |> stack
            end
            ω = ToolsArray(ω, (Structure(structures),)) |> stack
        end
        outdict[string(stimulus)] = Dict("ω" => ω)
    end
    out = []
    GC.gc()
    return outdict
end

begin
    f = TwoPanel()
    gs = subdivide(f, 1, 3)
    layerints = load(calcdir("plots", "grand_unified_layers.jld2"), "layerints")

    titles = ["Natural images post-onset", "Flashes post-onset"]

    for (i, stimulus) in enumerate(["r\"Natural_Images\"", "flash_250ms"])
        ax = Axis(gs[i], xlabel = "Average frequency (Hz)", ylabel = "Cortical depth (%)",
                  title = "$(titles[i])",
                  ytickformat = depthticks, yreversed = true,
                  limits = ((4.7, 7), (0, 1)))
        ω = plot_data[stimulus]["ω"]

        plotlayerints!(ax, layerints; axis = :y, flipside = false, newticks = false,
                       bgcolor = Makie.RGBA(0, 0, 0, 0))

        map(lookup(ω, Structure), eachslice(ω, dims = Structure)) do structure, x
            x = ustripall(x ./ (2π))

            μ, (σl, σh) = bootstrapmedian(x, dims = 2)

            μ, σl, σh = upsample.((μ, σl, σh), 5)

            band!(ax, Point2f.(collect(σl), lookup(μ, 1)),
                  Point2f.(collect(σh), lookup(μ, 1));
                  color = (structurecolormap[structure], 0.32), label = structure)

            lines!(ax, collect(μ), lookup(μ, 1);
                   color = (structurecolormap[structure], 0.8),
                   label = structure)
        end
    end
    Legend(gs[3], contents(gs[2]) |> only; merge = true, title = "Structure")
    addlabels!(f, labelformat)
    display(f)
    wsave(plotdir("figS2", "figS2.pdf"), f)
end
