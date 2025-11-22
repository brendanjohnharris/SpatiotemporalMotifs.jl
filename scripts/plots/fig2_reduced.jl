#! /bin/bash
# -*- mode: julia -*-
#=
exec julia +1.10.10 -t auto --color=yes "${BASH_SOURCE[0]}" "$@"
=#
using DrWatson
@quickactivate "SpatiotemporalMotifs"

using SpatiotemporalMotifs
import TimeseriesTools: freqs
using Peaks
using FileIO
using SpatiotemporalMotifs
import SpatiotemporalMotifs.layers
import SpatiotemporalMotifs.PTHR
using Random
using Distributed
@preamble
set_theme!(foresight(:physics))

stimuli = ["r\"Natural_Images\"", "spontaneous", "flash_250ms"]
xtickformat = terseticks
theta = Interval(SpatiotemporalMotifs.THETA()...)
gamma = Interval(SpatiotemporalMotifs.GAMMA()...)
alpha = 0.8
bandalpha = 0.2
mkpath(plotdir("fig2_reduced"))

data_file = savepath("fig2.jld2")(Dict())
plot_data = load(datadir("plots", data_file))

begin
    structure = "VISp"
    layer = 2
end

for stimulus in stimuli
    @info "Plotting spectra for $stimulus"
    @unpack S, layernames, layernums, layerints, meanlayers, S̄, oursessions, Q = plot_data[string(stimulus)]
    begin
        filebase = stimulus == "spontaneous" ? "" : "_$(val_to_string(stimulus))"
        f = TwoPanel()

        begin # * Mean power spectrum. Bands show 1 S.D.
            axargs = (; xscale = log10, yscale = log10,
                      limits = ((3, 300), (exp10(-15), exp10(-8.5))),
                      yticks = (exp10.([-14, -12, -10]),
                                [
                                    rich("10", superscript("-14")),
                                    rich("10", superscript("-12")),
                                    rich("10", superscript("-10"))
                                ]),
                      xlabel = "Frequency (Hz)",
                      xgridvisible = true,
                      ygridvisible = true,
                      xtickformat,
                      xticks = [3, 10, 30, 100],
                      ylabel = "Mean power spectral density (a.u.)",
                      title = "Power spectral density in $structure, L2/3")
            ax2 = Axis(f[1, 1]; axargs...) # For band annotations
            hideyaxis!(ax2)
            hidexaxis!(ax2)
            hidespines!(ax2)
            hidedecorations!(ax2)
            ax = Axis(f[1, 1]; axargs...)

            # # * Band annotations
            # vspan!(ax2, extrema(theta)..., color = (crimson, 0.22),
            #        label = "𝛉")
            # vlines!(ax2, [7.0], color = (crimson, 0.42), linestyle = :dash, linewidth = 4)
            # vspan!(ax2, extrema(gamma)..., color = (cornflowerblue, 0.22),
            #        label = "𝛄")
            # vlines!(ax2, [54.4], color = (cornflowerblue, 0.42), linestyle = :dash,
            #         linewidth = 4)

            s = S̄[Structure = At(structure)][layer = At([layer])][𝑓 = 3u"Hz" .. 300u"Hz"]
            # s = s ./ 10^((i - 1.5) / 2.5)
            yc = only(mean(s[𝑓 = Near(3u"Hz")]))
            ps, α = plotspectrum!(ax, s; textposition = (3, yc),
                                  color = structurecolormap[structure],
                                  annotations = [:peaks],
                                  label = structure)

            text!(ax, [100], [1e-10]; text = "b = $(round(α; sigdigits=3))", fontsize = 20)

            # ps, α = first.(psa), last.(psa)
            # α = round.(α, sigdigits = 3)
            # axislegend(ax2, position = :lb, labelsize = 12, backgroundcolor = :white,
            #            framevisible = true, padding = (5, 5, 5, 5))

            # axislegend(ax, position = :rt, labelsize = 12, merge = true,
            #            backgroundcolor = :white,
            #            framevisible = true, padding = (5, 5, 5, 5),
            #            nbanks = 3, patchsize = (15, 15), rowgap = 2)
            f
        end

        _b = map(2:5) do layer
            a = plot_data[stimulus]["fooof"]["χ"][Structure = At(structure)]
            S = plot_data[stimulus]["S"][Structure = At(structure)]
            layermap = S.metadata[:layerinfo]
            layermap = S.metadata[:streamlinedepths][layermap[3] .== layer]
            a = a[Depth = Near(layermap)]
            a = mean(a, dims = Depth)
            a = dropdims(a, dims = Depth)
        end
        _b = ToolsArray(_b, (Dim{:layer}(2:5),)) |> stack

        begin # * Plot
            begin
                ax = Axis(f[1, 2]; ylabel = "Cortical layer",
                          xlabel = "Mean spectral exponent",
                          title = "Spectral exponent in VISp", yreversed = true,
                          yticks = (2:5, ["L2/3", "L4", "L5", "L6"]))

                vlines!(ax, 2.0; color = :gray, linewidth = 3, linestyle = :dash)
            end

            begin # * Do statistical test against value of 2.0
                𝑝 = map(eachslice(_b, dims = :layer)) do b
                    𝑝 = HypothesisTests.SignedRankTest(collect(b) .- 2.0)
                    𝑝 = HypothesisTests.pvalue(𝑝, tail = :left)
                end
                𝑝[:] .= MultipleTesting.adjust(collect(𝑝)[:], BenjaminiHochberg())
            end
            for l in lookup(_b, :layer)
                if l > 1
                    x = _b[layer = At(l)]
                    p = 𝑝[layer = At(l)]
                    x = filter(!isnan, x)
                    if p < SpatiotemporalMotifs.PTHR
                        boxplot!(ax, fill(l, length(x)),
                                 collect(x)[:];
                                 orientation = :horizontal, show_outliers = false,
                                 color = Foresight.colororder[l - 1],
                                 strokecolor = Foresight.colororder[l - 1], strokewidth = 3)
                    else
                        boxplot!(ax, fill(l, length(x)),
                                 collect(x)[:];
                                 orientation = :horizontal, show_outliers = false,
                                 color = :transparent,
                                 strokecolor = Foresight.colororder[l - 1], strokewidth = 3)
                    end
                end
            end
        end
        display(f)
    end
    wsave(plotdir("fig2_reduced", "power_spectra$filebase.pdf"), f)
end
