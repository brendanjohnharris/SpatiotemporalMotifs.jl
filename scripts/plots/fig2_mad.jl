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

stimuli = [r"Natural_Images", "spontaneous", "flash_250ms"]
xtickformat = terseticks
theta = Interval(SpatiotemporalMotifs.THETA()...)
gamma = Interval(SpatiotemporalMotifs.GAMMA()...)
alpha = 0.8
bandalpha = 0.2
mkpath(plotdir("madev"))

plot_data, data_file = produce_or_load(Dict(), calcdir("plots");
                                       filename = savepath("madev")) do _
    session_table = load(calcdir("plots", "posthoc_session_table.jld2"), "session_table")
    oursessions = session_table.ecephys_session_id
    path = calcdir("madev")
    QQ = calcquality(path)
    plot_data = map(stimuli) do stimulus
        Q = QQ[stimulus = At(stimulus),
               Structure = At(structures),
               SessionID(At(oursessions))]
        @assert mean(Q) > 0.9

        begin # * Load data
            M = map(lookup(Q, Structure)) do structure
                out = map(lookup(Q, SessionID)) do sessionid
                    if Q[SessionID = At(sessionid), Structure = At(structure)] == 0
                        return nothing
                    end
                    filename = savepath((@strdict sessionid structure stimulus), "jld2",
                                        path)
                    mad = load(filename, "mad")
                    coeffs = load(filename, "coeffs")
                    return mad, coeffs
                end
                out = filter(!isnothing, out)
                out = filter(out) do x # Remove sessions that don't have data to a reasonable depth
                    maximum(DimensionalData.metadata(x[1])[:streamlinedepths]) > 0.90
                end

                coeffs = last.(out)
                out = first.(out)

                m = DimensionalData.metadata.(out)
                sessions = getindex.(m, :sessionid)

                streamlinedepths = getindex.(m, :streamlinedepths)
                layerinfo = getindex.(m, :layerinfo)

                unidepths = commondepths(streamlinedepths)
                out = map(out, streamlinedepths, layerinfo) do o, s, l
                    o = set(o, Depth => s)
                    layernames = ToolsArray(l[1], (Depth(lookup(o, Depth)),))
                    layernums = ToolsArray(l[3], (Depth(lookup(o, Depth)),))
                    o = o[Depth(Near(unidepths))]
                    layernames = layernames[Depth(Near(unidepths))]
                    layernums = layernums[Depth(Near(unidepths))]
                    # @assert length(unique(lookup(o, Depth))) == length(unidepths)
                    @assert issorted(lookup(o, Depth))
                    push!(o.metadata, :layernames => layernames)
                    push!(o.metadata, :layernums => layernums)
                    o = set(o, Depth => unidepths)
                end
                coeffs = map(coeffs, streamlinedepths, layerinfo) do o, s, l
                    o = set(o, Depth => s)
                    layernames = ToolsArray(l[1], (Depth(lookup(o, Depth)),))
                    layernums = ToolsArray(l[3], (Depth(lookup(o, Depth)),))
                    o = o[Depth(Near(unidepths))]
                    layernames = layernames[Depth(Near(unidepths))]
                    layernums = layernums[Depth(Near(unidepths))]
                    # @assert length(unique(lookup(o, Depth))) == length(unidepths)
                    @assert issorted(lookup(o, Depth))
                    # push!(o.metadata, :layernames => layernames)
                    # push!(o.metadata, :layernums => layernums)
                    o = set(o, Depth => unidepths)
                end
                layernames = ToolsArray(stack(getindex.(DimensionalData.metadata.(out),
                                                        :layernames)),
                                        (Dim{:depths}(unidepths), SessionID(sessions)))
                layernums = ToolsArray(stack(getindex.(DimensionalData.metadata.(out),
                                                       :layernums)),
                                       (Dim{:depths}(unidepths), SessionID(sessions)))
                mad = stack(SessionID(sessions), out, dims = 3) .|> Float32
                coeffs = stack(SessionID(sessions), coeffs, dims = 2) .|> Float32
                layernums = parselayernum.(layernames)
                return mad, coeffs, layernames, layernums
            end
            M, coeffs, layernames, layernums = zip(M...)
        end
        begin # * Format layers
            meanlayers = map(layernums) do l
                round.(Int, mean(l, dims = 2))
            end
            M̄ = map(M, meanlayers) do m, l
                # s = set(s, Depth => Dim{:layer}(layernum2name.(parent(l)[:])))
                m = set(m, Depth => Dim{:layer}(parent(l)[:]))
                m = set(m, :layer => DimensionalData.Unordered)
            end
            M̄ = map(M̄) do s
                ss = map(unique(lookup(s, :layer))) do l
                    ls = s[Dim{:layer}(At(l))]
                    if hasdim(ls, :layer)
                        ls = mean(ls, dims = :layer)
                    end
                    ls
                end
                cat(ss..., dims = :layer)
            end
            M̄ = ToolsArray(M̄ |> collect, (Structure(lookup(Q, Structure)),))
            M = ToolsArray(M |> collect, (Structure(lookup(Q, Structure)),))
            # S̄ = map(S̄) do s
            #     N = UnitEnergy(s, dims = 1)
            #     N(s) .|> Float32
            # end # ! No normalization
            layerints = load(calcdir("plots", "grand_unified_layers.jld2"), "layerints")
        end

        coeffs = ToolsArray(collect(coeffs), (Structure(structures),))

        plot_data = @strdict M M̄ coeffs layernames layernums layerints meanlayers oursessions Q
        return plot_data
    end

    return Dict(string.(stimuli) .=> plot_data)
end

for stimulus in stimuli
    @info "Plotting spectra for $stimulus"
    @unpack M, coeffs, layernames, layernums, layerints, meanlayers, M̄, oursessions, Q = plot_data[string(stimulus)]
    begin
        filebase = stimulus == "spontaneous" ? "" : "_$(val_to_string(stimulus))"
        statsfile = plotdir("madev", "madev$filebase.txt")
        close(open(statsfile, "w")) # Create the file or clear it
        f = TwoPanel()

        # begin # * Mean power spectrum. Bands show 1 S.D.
        #     axargs = (; xscale = log10, yscale = log10,
        #               limits = ((3, 300), (exp10(-15), exp10(-8.5))),
        #               yticks = (exp10.([-14, -12, -10]),
        #                         [
        #                             rich("10", superscript("-14")),
        #                             rich("10", superscript("-12")),
        #                             rich("10", superscript("-10"))
        #                         ]),
        #               xlabel = "Frequency (Hz)",
        #               xgridvisible = true,
        #               ygridvisible = true,
        #               xgridstyle = :dash,
        #               ygridstyle = :dash,
        #               xtickformat,
        #               xticks = [3, 10, 30, 100],
        #               ylabel = "Mean power spectral density (a.u.)",
        #               title = "Power spectral density")
        #     ax2 = Axis(f[1, 1]; axargs...) # For band annotations
        #     hideyaxis!(ax2)
        #     hidexaxis!(ax2)
        #     hidespines!(ax2)
        #     hidedecorations!(ax2)
        #     ax = Axis(f[1, 1]; axargs...)

        #     # * Band annotations
        #     vspan!(ax2, extrema(theta)..., color = (crimson, 0.22),
        #            label = "𝛉")
        #     vlines!(ax2, [7.0], color = (crimson, 0.42), linestyle = :dash, linewidth = 4)
        #     vspan!(ax2, extrema(gamma)..., color = (cornflowerblue, 0.22),
        #            label = "𝛄")
        #     vlines!(ax2, [54.4], color = (cornflowerblue, 0.42), linestyle = :dash,
        #             linewidth = 4)

        #     psa = map(enumerate(structures)) do (i, structure)
        #         s = M̄[Structure = At(structure)][𝑡 = 1e-3 .. 1e-2]
        #         # s = s ./ 10^((i - 1.5) / 2.5)
        #         yc = only(mean(s[𝑓 = Near(3u"Hz")]))
        #         if i == 1 && stimulus != r"Natural_Images"
        #             plotspectrum!(ax, s; textposition = (3, yc),
        #                           color = structurecolors[i], annotations = [:peaks],
        #                           label = structure)
        #         else
        #             plotspectrum!(ax, s; textposition = (14, yc),
        #                           color = structurecolors[i], annotations = [],
        #                           label = structure)
        #         end
        #     end
        #     ps, α = first.(psa), last.(psa)
        #     α = round.(α, sigdigits = 3)
        #     axislegend(ax2, position = :lb, labelsize = 12, backgroundcolor = :white,
        #                framevisible = true, padding = (5, 5, 5, 5))
        #     # leg = ["$s (α = $α)" for (s, α) in zip(structures, α)]
        #     # map(enumerate(leg)) do (i, s)
        #     #     if i > 3
        #     #         text!(ax, 290, exp10(-1.1 - ((i - 3) / 3 - 1 - 0.1)); text = s,
        #     #               color = structurecolors[i],
        #     #               fontsize = 14,
        #     #               align = (:right, :bottom))
        #     #     else
        #     #         text!(ax, 50, exp10(-1.1 - (i / 3 - 1 - 0.1)); text = s,
        #     #               color = structurecolors[i],
        #     #               fontsize = 14,
        #     #               align = (:right, :bottom))
        #     #     end
        #     # end

        #     axislegend(ax, position = :rt, labelsize = 12, merge = true,
        #                backgroundcolor = :white,
        #                framevisible = true, padding = (5, 5, 5, 5),
        #                nbanks = 3, patchsize = (15, 15), rowgap = 2)
        #     f
        # end

        begin # * Load the channel-wise fits
            χ = plot_data[string(stimulus)]["coeffs"]

            # * Ensure structures and sessions match
            coeffs = getindex.(coeffs, [SessionID(At(oursessions))])
            coeffs = coeffs[Structure = At(structures)]
        end

        # begin # * Plot the exponent for each subject in VISl (as a check)
        #     chi = b[2] # VISl
        #     chi = chi[:, sortperm(eachcol(chi))][3:end, 1:5:end]
        #     chi = chi ./ median(chi, dims = Depth)
        #     ff = Figure()
        #     ax = Axis(ff[1, 1])
        #     scatter!.([ax], eachcol(chi))
        #     ff
        # end

        begin # * Plot the intercept
            ax = Axis(f[1, 1]; ylabel = "Cortical depth (%)",
                      xlabel = "Median diffusion exponent",
                      limits = ((0.25, 0.75), (0, 1)), ytickformat = depthticks,
                      title = "Diffusion exponent", yreversed = true)
            vlines!(ax, 0.5; color = :gray, linewidth = 3, linestyle = :dash)

            for (i, _b) in coeffs |> enumerate |> collect |> reverse
                μ, (σl, σh) = bootstrapmedian(_b .+ eps() .* randn(size(_b)),
                                              dims = SessionID)
                μ, σl, σh = upsample.((μ, σl, σh), 5)

                band!(ax, Point2f.(collect(σl), lookup(μ, 1)),
                      Point2f.(collect(σh), lookup(μ, 1));
                      color = (structurecolors[i], 0.32), label = structures[i])
                lines!(ax, collect(μ), lookup(μ, 1); color = (structurecolors[i], alpha),
                       label = structures[i])
            end
            l = axislegend(ax, position = :lb, nbanks = 2, labelsize = 12, merge = true,
                           margin = (30, 0, 10, 0))
            reverselegend!(l)
            plotlayerints!(ax, layerints; axis = :y, newticks = false)
        end

        begin # * Is the exponent correlated to depth?
            tps = [SpatiotemporalMotifs.mediankendallpvalue(lookup(x, Depth), x) for x in χ]
            τs = cat(first.(tps), last.(tps); dims = Var([:kendall, :pvalue]))
            mtau = median(τs[2:end, 1]) # Median excluding VISp
            open(statsfile, "a+") do file
                write(file, "\n## Diffusion exponent\n")
                show(file, DataFrame(DimTable(τs; layersfrom = Var)))
                write(file, "\nMedian τ (no VISp) = $mtau\n")
            end
        end

        # begin # * Plot fooof residuals. Bands are 1 S.D.
        #     # f = Figure()
        #     Sr_log = map(ustripall.(S), L, meanlayers) do s, l, m
        #         s = deepcopy(s)
        #         l = deepcopy(l)
        #         idxs = indexin(lookup(s, SessionID), lookup(l, SessionID))
        #         l = l[:, idxs] # Match sessions just in case
        #         map(eachslice(s, dims = (Depth, SessionID)), l) do s, l
        #             _s = log10.(ustripall(s))
        #             s .= _s .- (_s |> freqs .|> l .|> log10)
        #         end
        #         s = set(s, Depth => Dim{:layer}(layernum2name.(parent(m)[:])))
        #         s = set(s, :layer => DimensionalData.Irregular)
        #     end

        #     Sr_log = ToolsArray(Sr_log |> collect,
        #                         (Structure(lookup(Q, Structure)),))

        #     for (i, structure) in enumerate(["VISl"])
        #         ax2 = Axis(f[1, i + 1]; xscale = log10,
        #                    limits = ((3, 300), (-0.1, 3.5)), xtickformat,
        #                    xlabel = "Frequency (Hz)",
        #                    ylabel = "Residual spectral density (dB)",
        #                    title = "Residual spectral density in " * structure,
        #                    xticks = [3, 10, 30, 100]) # xticksvisible = false, yaxisposition = :right,
        #         #    xticklabelsvisible = false,
        #         vlines!(ax2, [7.0], color = (crimson, 0.42), linestyle = :dash,
        #                 linewidth = 4)
        #         vlines!(ax2, [54.4], color = (cornflowerblue, 0.42), linestyle = :dash,
        #                 linewidth = 4)
        #         vspan!(ax, extrema(theta)..., color = (crimson, 0.22),
        #                label = "𝛉 ($(theta.left) – $(theta.right) Hz)")
        #         vspan!(ax, extrema(gamma)..., color = (cornflowerblue, 0.22),
        #                label = "𝛄 ($(gamma.left) – $(gamma.right) Hz)")

        #         for (i, (c, l)) in (reverse ∘ collect ∘ enumerate ∘ zip)(layercolors,
        #                                                                  layers)
        #             s = Sr_log[Structure = At(structure)][Freq(3 .. 300)]
        #             s = s[layer = (lookup(s, :layer) .== [l])]
        #             s = dropdims(nansafe(mean; dims = :layer)(s), dims = :layer)
        #             d = (length(layers) - i + 1) / 2
        #             hlines!(ax2, [d]; color = (c, 0.22), linestyle = :dash)
        #             μ = dropdims(mean(s, dims = SessionID), dims = SessionID) .+
        #                 d
        #             σ = dropdims(std(s, dims = SessionID), dims = SessionID)
        #             band!(ax2, TimeseriesTools.freqs(μ), collect(μ .- σ), collect(μ .+ σ);
        #                   color = (c, bandalpha))
        #             lines!(ax2, TimeseriesTools.freqs(μ), collect(μ); color = (c, alpha))
        #         end
        #         C = Colorbar(f[1, i + 1][1, 2],
        #                      colormap = reverse(cgrad(layercolors, categorical = true)),
        #                      ticks = (range(0, 1, length = (2 * length(layercolors) + 1))[2:2:end],
        #                               ["L"] .* reverse(layers)),
        #                      ticklabelrotation = 0,# π / 2,
        #                      ticklabelsize = 13)
        #         # linkxaxes!(ax, ax2)
        #     end
        #     f
        # end

        begin # * Plot the hierarchical correlation across layers
            N = 10000
            method = :group
            unidepths = commondepths(lookup.(χ, [Depth]))
            x = getindex.([SpatiotemporalMotifs.hierarchy_scores], structures)

            unichi = getindex.(coeffs, [Depth(Near(unidepths))])
            unichi = set.(unichi, [Depth => unidepths])
            y = stack(Structure(structures), unichi)

            μ, σ, 𝑝 = hierarchicalkendall(x, y, method; N)
        end
        begin
            # μ[𝑝 .> PTHR] .= NaN
            # μt[𝑝t .> PTHR] .= NaN
            # μg[𝑝g .> PTHR] .= NaN

            markersize = 10

            ax = Axis(f[1, 2]; ylabel = "Cortical depth (%)",
                      xlabel = "Kendall's 𝜏",
                      ytickformat = depthticks,
                      xtickformat,
                      title = "Diffusion exponent gradient",
                      limits = ((-0.6, 0.4), (0, 1)),
                      yreversed = true,
                      yticklabelsvisible = false)

            vlines!(ax, 0; color = :gray, linewidth = 3, linestyle = :dash)

            band!(ax, Point2f.(collect(first.(σ)), unidepths),
                  Point2f.(collect(last.(σ)), unidepths);
                  color = (:cornflowerblue, bandalpha),
                  label = "Diffusion exponent")
            # lines!(ax, unidepths, collect(μ); alpha = bandalpha,
            #        label = "1/𝑓 exponent", color = crimson)
            scatter!(ax, collect(μ[𝑝 .< PTHR]), unidepths[𝑝 .< PTHR];
                     label = "Diffusion exponent", color = :cornflowerblue, markersize)
            scatter!(ax, collect(μ[𝑝 .≥ PTHR]), unidepths[𝑝 .≥ PTHR]; color = :transparent,
                     strokecolor = :cornflowerblue,
                     strokewidth = 1, markersize)

            # axislegend(ax, position = :lb, merge = true, labelsize = 12, nbanks = 1)

            plotlayerints!(ax, layerints; axis = :y, newticks = false)
        end
        # addlabels!(f, ["a", "d", "f", "e", "c", "b", "g", "h"])
        f |> display
    end
    wsave(plotdir("madev", "madev$filebase.pdf"), f)
end
