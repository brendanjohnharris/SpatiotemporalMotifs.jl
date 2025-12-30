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
mkpath(plotdir("fig2"))

if !isfile(calcdir("plots", savepath("fig2", Dict(), "jld2")))
    if nprocs() == 1
        if SpatiotemporalMotifs.CLUSTER()
            using USydClusters
            ourprocs = USydClusters.Physics.addprocs(30; mem = 6, ncpus = 1,
                                                     project = projectdir(),
                                                     queue = "l40s")
        else
            addprocs(19)
        end
    end
    @everywhere using SpatiotemporalMotifs
    @everywhere SpatiotemporalMotifs.@preamble
end

plot_data, data_file = produce_or_load(Dict(), calcdir("plots");
                                       filename = savepath("fig2")) do _
    session_table = load(calcdir("plots", "posthoc_session_table.jld2"), "session_table")
    oursessions = session_table.ecephys_session_id
    path = calcdir("power_spectra")
    QQ = calcquality(path)
    plot_data = map(stimuli) do stimulus
        Q = QQ[stimulus = At(stimulus),
               Structure = At(structures),
               SessionID(At(oursessions))]
        @assert mean(Q) > 0.9

        begin # * Load data
            S = map(lookup(Q, Structure)) do structure
                out = map(lookup(Q, SessionID)) do sessionid
                    if Q[SessionID = At(sessionid), Structure = At(structure)] == 0
                        return nothing
                    end
                    filename = savepath((@strdict sessionid structure stimulus), "jld2",
                                        path)
                    S = load(filename, "S")
                end
                out = filter(!isnothing, out)
                out = filter(out) do x # Remove sessions that don't have data to a reasonable depth
                    maximum(DimensionalData.metadata(x)[:streamlinedepths]) > 0.90
                end

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
                layernames = ToolsArray(stack(getindex.(DimensionalData.metadata.(out),
                                                        :layernames)),
                                        (Dim{:depths}(unidepths), SessionID(sessions)))
                layernums = ToolsArray(stack(getindex.(DimensionalData.metadata.(out),
                                                       :layernums)),
                                       (Dim{:depths}(unidepths), SessionID(sessions)))
                S = stack(SessionID(sessions), out, dims = 3) .|> Float32
                layernums = parselayernum.(layernames)
                return S, layernames, layernums
            end
            S, layernames, layernums = zip(S...)
        end
        begin # * Format layers
            meanlayers = map(layernums) do l
                round.(Int, mean(l, dims = 2))
            end
            S̄ = map(S, meanlayers) do s, l
                # s = set(s, Depth => Dim{:layer}(layernum2name.(parent(l)[:])))
                s = set(s, Depth => Dim{:layer}(parent(l)[:]))
                s = set(s, :layer => DimensionalData.Unordered)
            end
            S̄ = map(S̄) do s
                ss = map(unique(lookup(s, :layer))) do l
                    ls = s[Dim{:layer}(At(l))]
                    if hasdim(ls, :layer)
                        ls = mean(ls, dims = :layer)
                    end
                    ls
                end
                cat(ss..., dims = :layer)
            end
            S̄ = ToolsArray(S̄ |> collect, (Structure(lookup(Q, Structure)),))
            S = ToolsArray(S |> collect, (Structure(lookup(Q, Structure)),))
            # S̄ = map(S̄) do s
            #     N = UnitEnergy(s, dims = 1)
            #     N(s) .|> Float32
            # end # ! No normalization
            layerints = load(calcdir("plots", "grand_unified_layers.jld2"), "layerints")
        end

        L = @withprogress name="Fooof $stimulus" begin
            threadlog = -1
            threadmax = length(S)
            map(S) do s # If you can set up the cluster workers as above, each stimulus should take about 5 minutes. Otherwise, 30 minutes
                threadlog += 1
                @logprogress threadlog / threadmax
                pmap(SpatiotemporalMotifs.fooof, eachslice(ustripall(s), dims = (2, 3)))
            end
        end
        χ = [getindex.(last.(l), :χ) for l in L]
        b = [getindex.(last.(l), :b) for l in L]
        k = [getindex.(last.(l), :k) for l in L]
        χ, b, k = ToolsArray.([χ, b, k], [(Structure(structures),)])
        fooof = Dict("χ" => χ, "b" => b, "k" => k)

        plot_data = @strdict S layernames layernums layerints meanlayers S̄ oursessions Q fooof
        return plot_data
    end

    return Dict(string.(stimuli) .=> plot_data)
end

for stimulus in reverse(stimuli)
    @info "Plotting spectra for $stimulus"
    @unpack S, layernames, layernums, layerints, meanlayers, S̄, oursessions, Q = plot_data[string(stimulus)]
    is_flash = stimulus == "flash_250ms" # ?
    is_flash && mkpath(datadir("source_data", "fig2")) # ?
    begin
        filebase = stimulus == "spontaneous" ? "" : "_$(val_to_string(stimulus))"
        statsfile = plotdir("fig2", "power_spectra$filebase.txt")
        close(open(statsfile, "w")) # Create the file or clear it
        f = SixPanel()

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
                      xgridstyle = :dash,
                      ygridstyle = :dash,
                      xtickformat,
                      xticks = [3, 10, 30, 100],
                      ylabel = "Mean power spectral density (arb. units)",
                      title = "Power spectral density")
            ax2 = Axis(f[1, 1]; axargs...) # For band annotations
            hideyaxis!(ax2)
            hidexaxis!(ax2)
            hidespines!(ax2)
            hidedecorations!(ax2)
            ax = Axis(f[1, 1]; axargs...)

            # * Band annotations
            vspan!(ax2, extrema(theta)..., color = (crimson, 0.22),
                   label = "𝛉")
            vlines!(ax2, [7.0], color = (crimson, 0.42), linestyle = :dash, linewidth = 4)
            vspan!(ax2, extrema(gamma)..., color = (cornflowerblue, 0.22),
                   label = "𝛄")
            vlines!(ax2, [54.4], color = (cornflowerblue, 0.42), linestyle = :dash,
                    linewidth = 4)

            psa = map(enumerate(structures)) do (i, structure)
                s = S̄[Structure = At(structure)][𝑓 = 3u"Hz" .. 300u"Hz"]
                s = s ./ 10^((i - 1.5) / 2.5)
                yc = only(mean(s[𝑓 = Near(3u"Hz")]))
                if i == 1 && stimulus != r"Natural_Images"
                    plotspectrum!(ax, s; textposition = (3, yc),
                                  color = structurecolors[i], annotations = [:peaks],
                                  label = structure)
                else
                    plotspectrum!(ax, s; textposition = (14, yc),
                                  color = structurecolors[i], annotations = [],
                                  label = structure)
                end
            end
            ps, α = first.(psa), last.(psa)
            α = round.(α, sigdigits = 3)
            axislegend(ax2, position = :lb, labelsize = 12, backgroundcolor = :white,
                       framevisible = true, padding = (5, 5, 5, 5))
            # leg = ["$s (α = $α)" for (s, α) in zip(structures, α)]
            # map(enumerate(leg)) do (i, s)
            #     if i > 3
            #         text!(ax, 290, exp10(-1.1 - ((i - 3) / 3 - 1 - 0.1)); text = s,
            #               color = structurecolors[i],
            #               fontsize = 14,
            #               align = (:right, :bottom))
            #     else
            #         text!(ax, 50, exp10(-1.1 - (i / 3 - 1 - 0.1)); text = s,
            #               color = structurecolors[i],
            #               fontsize = 14,
            #               align = (:right, :bottom))
            #     end
            # end

            axislegend(ax, position = :rt, labelsize = 12, merge = true,
                       backgroundcolor = :white,
                       framevisible = true, padding = (5, 5, 5, 5),
                       nbanks = 3, patchsize = (15, 15), rowgap = 2)
            f
        end
        if is_flash # ?
            outdir = datadir("source_data", "fig2") # ?
            df = DataFrame() # ?
            for (i, structure) in enumerate(structures) # ?
                sname = replace(structure, "/" => "") # ?
                s = S̄[Structure = At(structure)][𝑓 = 3u"Hz" .. 300u"Hz"] # ?
                s = s ./ 10^((i - 1.5) / 2.5) # ?
                μ = mean(s, dims = (SessionID, :layer)) # ?
                μ = dropdims(μ, dims = (SessionID, :layer)) |> ustripall # ?
                σ = std(s, dims = (SessionID, :layer)) ./ 2 # ?
                σ = dropdims(σ, dims = (SessionID, :layer)) |> ustripall # ?
                if i == 1 # ?
                    df[!, "Frequency (Hz)"] = lookup(μ, 𝑓) # ?
                end # ?
                df[!, "Mean power $sname (arb. units)"] = collect(μ) # ?
                df[!, "SD $sname (arb. units)"] = collect(σ) # ?
            end # ?
            CSV.write(joinpath(outdir, "panel_a_power_spectrum.csv"), df) # ?
        end # ?

        begin # * Load the channel-wise fits
            χ = plot_data[string(stimulus)]["fooof"]["χ"]
            b = plot_data[string(stimulus)]["fooof"]["b"]
            k = plot_data[string(stimulus)]["fooof"]["k"]
            L = SpatiotemporalMotifs.ramap(fooof, χ, b, k) # Fit functions

            # * Ensure structures and sessions match
            L = getindex.(L, [SessionID(At(oursessions))])
            L = L[Structure = At(structures)]
            χ = getindex.(χ, [SessionID(At(oursessions))])
            χ = χ[Structure = At(structures)]
            b = getindex.(b, [SessionID(At(oursessions))])
            b = b[Structure = At(structures)]

            map(b) do _b
                map(eachslice(_b, dims = SessionID)) do a
                    N = nansafe(RobustZScore)
                    N = fit(N, a)
                    normalize!(a, N)
                end
            end

            @assert all(last.(size.(L)) .≥ last.(size.(S))) # Check we have residuals for all sessions we have spectra for
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
            ax = Axis(f[3, 1:2][1, 2]; #ylabel = "Cortical depth (%)",
                      xlabel = "Normalized 1/𝑓 intercept",
                      limits = ((-2.75, 2.75), (0, 1)), ytickformat = depthticks,
                      title = "1/𝑓 intercept", yreversed = true, yticklabelsvisible = false)
            for (i, _b) in b |> enumerate |> collect |> reverse
                μ, (σl, σh) = bootstrapmedian(_b .+ eps() .* randn(size(_b)),
                                              dims = SessionID)
                μ, σl, σh = upsample.((μ, σl, σh), 5)

                band!(ax, Point2f.(collect(σl), lookup(μ, 1)),
                      Point2f.(collect(σh), lookup(μ, 1));
                      color = (structurecolors[i], 0.32), label = structures[i])
                lines!(ax, collect(μ), lookup(μ, 1); color = (structurecolors[i], alpha),
                       label = structures[i])
            end
            # l = axislegend(ax, position = :lb, nbanks = 2, labelsize = 12, merge = true)
            # reverselegend!(l)
            plotlayerints!(ax, layerints; axis = :y, newticks = false, flipside = true)
        end
        if is_flash # ?
            outdir = datadir("source_data", "fig2") # ?
            df = DataFrame() # ?
            for (i, _b) in enumerate(b) # ?
                sname = replace(structures[i], "/" => "") # ?
                μ, (σl, σh) = bootstrapmedian(_b .+ eps() .* randn(size(_b)), # ?
                                              dims = SessionID) # ?
                μ, σl, σh = upsample.((μ, σl, σh), 5) # ?
                if i == 1 # ?
                    df[!, "Cortical depth (%)"] = lookup(μ, 1) # ?
                end # ?
                df[!, "Median intercept $sname (normalized)"] = collect(μ) # ?
                df[!, "CI lower $sname"] = collect(σl) # ?
                df[!, "CI upper $sname"] = collect(σh) # ?
            end # ?
            CSV.write(joinpath(outdir, "panel_c_intercept.csv"), df) # ?
        end # ?

        begin # * Plot the exponent
            ax = Axis(f[3, 1:2][1, 1]; ylabel = "Cortical depth (%)",
                      xlabel = "1/𝑓 exponent",
                      limits = ((0.9, 2.1), (0, 1)), ytickformat = depthticks,
                      xtickformat,
                      title = "1/𝑓 exponent", yreversed = true)
            for (i, chi) in χ |> enumerate |> collect |> reverse
                μ, (σl, σh) = bootstrapmedian(chi, dims = SessionID)
                μ, σl, σh = upsample.((μ, σl, σh), 5)

                band!(ax, Point2f.(collect(σl), lookup(μ, 1)),
                      Point2f.(collect(σh), lookup(μ, 1));
                      color = (structurecolors[i], 0.32), label = structures[i])
                lines!(ax, collect(μ), lookup(μ, 1); color = (structurecolors[i], alpha),
                       label = structures[i])
            end
            # l = axislegend(ax, position = :lb, nbanks = 2, labelsize = 12, merge = true)
            # reverselegend!(l)
            plotlayerints!(ax, layerints; axis = :y, newticks = false, flipside = true)
        end
        if is_flash # ?
            outdir = datadir("source_data", "fig2") # ?
            df = DataFrame() # ?
            for (i, chi) in enumerate(χ) # ?
                sname = replace(structures[i], "/" => "") # ?
                μ, (σl, σh) = bootstrapmedian(chi, dims = SessionID) # ?
                μ, σl, σh = upsample.((μ, σl, σh), 5) # ?
                if i == 1 # ?
                    df[!, "Cortical depth (%)"] = lookup(μ, 1) # ?
                end # ?
                df[!, "Median exponent $sname"] = collect(μ) # ?
                df[!, "CI lower $sname"] = collect(σl) # ?
                df[!, "CI upper $sname"] = collect(σh) # ?
            end # ?
            CSV.write(joinpath(outdir, "panel_b_exponent.csv"), df) # ?
        end # ?
        begin # * Is the exponent correlated to depth?
            tps = [SpatiotemporalMotifs.mediankendallpvalue(lookup(x, Depth), x) for x in χ]
            τs = cat(first.(tps), last.(tps); dims = Var([:kendall, :pvalue]))
            mtau = median(τs[2:end, 1]) # Median excluding VISp
            open(statsfile, "a+") do file
                write(file, "\n## 1/𝑓 exponent\n")
                show(file, DataFrame(DimTable(τs; layersfrom = Var)))
                write(file, "\nMedian τ (no VISp) = $mtau\n")
            end
        end

        begin # * Relative power (supplement)
            fs = SixPanel()
            gs = subdivide(fs, 3, 2)
            map(enumerate(S)) do (i, x)
                x = mapslices(x, dims = (1, 2)) do s
                    s ./ maximum(s, dims = 2)
                end |> ustripall
                x = median(x[𝑓 = 1 .. 300], dims = SessionID)
                x = dropdims(x, dims = SessionID)[:, 2:end]

                ax = Axis(gs[i][1, 1], xlabel = "Frequency (Hz)",
                          ylabel = "Cortical depth (%)",
                          ytickformat = depthticks,
                          limits = ((1, 200), (nothing, nothing)), yreversed = true,
                          aspect = 1,
                          title = metadata(x)[:structure], xticks = [1, 100, 200])
                x = upsample(x, 5, 2)
                p = heatmap!(ax, lookup(x, 1), lookup(x, 2) .* 100, x |> collect,
                             colormap = :viridis, colorrange = (0, 1), rasterize = 5)
                if i % 2 == 0
                    Colorbar(gs[i][1, 2], p, label = "Relative power (arb. units)",
                             tellheight = true)
                    hideyaxis!(ax)
                end
                if i < 5
                    hidexaxis!(ax)
                end
            end
            # rowsize!(fs.layout, 1, Aspect(2, 0.9))
            # rowsize!(fs.layout, 2, Aspect(2, 0.9))
            # rowsize!(fs.layout, 3, Aspect(2, 0.9))
            fs
            wsave(plotdir("fig2",
                          "relative_power_supplement$filebase.pdf"),
                  fs)
        end

        begin # * Plot fooof residuals. Bands are 1 S.D.
            # f = Figure()
            Sr_log = map(ustripall.(S), L, meanlayers) do s, l, m
                s = deepcopy(s)
                l = deepcopy(l)
                idxs = indexin(lookup(s, SessionID), lookup(l, SessionID))
                l = l[:, idxs] # Match sessions just in case
                map(eachslice(s, dims = (Depth, SessionID)), l) do s, l
                    _s = log10.(ustripall(s))
                    s .= _s .- (_s |> freqs .|> l .|> log10)
                end
                s = set(s, Depth => Dim{:layer}(layernum2name.(parent(m)[:])))
                s = set(s, :layer => DimensionalData.Irregular)
            end

            Sr_log = ToolsArray(Sr_log |> collect,
                                (Structure(lookup(Q, Structure)),))

            for (i, structure) in enumerate(["VISl"])
                ax2 = Axis(f[1, i + 1]; xscale = log10,
                           limits = ((3, 300), (-0.1, 3.5)), xtickformat,
                           xlabel = "Frequency (Hz)",
                           ylabel = "Residual spectral density (dB)",
                           title = "Residual spectral density in " * structure,
                           xticks = [3, 10, 30, 100]) # xticksvisible = false, yaxisposition = :right,
                #    xticklabelsvisible = false,
                vlines!(ax2, [7.0], color = (crimson, 0.42), linestyle = :dash,
                        linewidth = 4)
                vlines!(ax2, [54.4], color = (cornflowerblue, 0.42), linestyle = :dash,
                        linewidth = 4)
                vspan!(ax, extrema(theta)..., color = (crimson, 0.22),
                       label = "𝛉 ($(theta.left) – $(theta.right) Hz)")
                vspan!(ax, extrema(gamma)..., color = (cornflowerblue, 0.22),
                       label = "𝛄 ($(gamma.left) – $(gamma.right) Hz)")

                for (i, (c, l)) in (reverse ∘ collect ∘ enumerate ∘ zip)(layercolors,
                                                                         layers)
                    s = Sr_log[Structure = At(structure)][Freq(3 .. 300)]
                    s = s[layer = (lookup(s, :layer) .== [l])]
                    s = dropdims(nansafe(mean; dims = :layer)(s), dims = :layer)
                    d = (length(layers) - i + 1) / 2
                    hlines!(ax2, [d]; color = (c, 0.22), linestyle = :dash)
                    μ = dropdims(mean(s, dims = SessionID), dims = SessionID) .+
                        d
                    σ = dropdims(std(s, dims = SessionID), dims = SessionID)
                    band!(ax2, TimeseriesTools.freqs(μ), collect(μ .- σ), collect(μ .+ σ);
                          color = (c, bandalpha))
                    lines!(ax2, TimeseriesTools.freqs(μ), collect(μ); color = (c, alpha))
                end
                C = Colorbar(f[1, i + 1][1, 2],
                             colormap = reverse(cgrad(layercolors, categorical = true)),
                             ticks = (range(0, 1, length = (2 * length(layercolors) + 1))[2:2:end],
                                      ["L"] .* reverse(layers)),
                             ticklabelrotation = 0,# π / 2,
                             ticklabelsize = 13)
                # linkxaxes!(ax, ax2)
            end
            f
        end
        if is_flash # ?
            outdir = datadir("source_data", "fig2") # ?
            structure = "VISl" # ?
            df = DataFrame() # ?
            for (i, l) in enumerate(layers) # ?
                lname = replace(l, "/" => "") # ?
                s = Sr_log[Structure = At(structure)][Freq(3 .. 300)] # ?
                s = s[layer = (lookup(s, :layer) .== [l])] # ?
                s = dropdims(nansafe(mean; dims = :layer)(s), dims = :layer) # ?
                μ = dropdims(mean(s, dims = SessionID), dims = SessionID) # ?
                σ = dropdims(std(s, dims = SessionID), dims = SessionID) # ?
                if i == 1 # ?
                    df[!, "Frequency (Hz)"] = TimeseriesTools.freqs(μ) # ?
                end # ?
                df[!, "Mean residual L$lname (dB)"] = collect(μ) # ?
                df[!, "SD L$lname (dB)"] = collect(σ) # ?
            end # ?
            CSV.write(joinpath(outdir, "panel_d_residual_spectrum_VISl.csv"), df) # ?
        end # ?

        begin # * Residual power supplement
            sf = SixPanel()
            gs = subdivide(sf, 3, 2)

            for (i, structure) in enumerate(structures)
                ax2 = Axis(gs[i]; xscale = log10,
                           limits = ((3, 300), (-0.1, 3.5)), xtickformat,
                           xlabel = "Frequency (Hz)",
                           ylabel = "Residual power (dB)", title = structure,
                           xticks = [3, 10, 30, 100]) # xticksvisible = false, yaxisposition = :right,
                #    xticklabelsvisible = false,
                vlines!(ax2, [7.0], color = (crimson, 0.42), linestyle = :dash,
                        linewidth = 4)
                vlines!(ax2, [54.4], color = (cornflowerblue, 0.42), linestyle = :dash,
                        linewidth = 4)

                for (i, (c, l)) in (reverse ∘ collect ∘ enumerate ∘ zip)(layercolors,
                                                                         layers)
                    s = Sr_log[Structure = At(structure)][Freq(3 .. 300)]
                    s = s[layer = (lookup(s, :layer) .== [l])]
                    s = dropdims(mean(s, dims = :layer), dims = :layer)
                    d = (length(layers) - i + 1) / 2
                    hlines!(ax2, [d]; color = (c, 0.22), linestyle = :dash)
                    μ = dropdims(mean(s, dims = SessionID), dims = SessionID) .+
                        d
                    σ = dropdims(std(s, dims = SessionID), dims = SessionID)
                    band!(ax2, TimeseriesTools.freqs(μ), collect(μ .- σ), collect(μ .+ σ);
                          color = (c, bandalpha))
                    lines!(ax2, TimeseriesTools.freqs(μ), collect(μ); color = (c, alpha))
                end
                C = Colorbar(gs[i][1, 2],
                             colormap = reverse(cgrad(layercolors, categorical = true)),
                             ticks = (range(0, 1, length = (2 * length(layercolors) + 1))[2:2:end],
                                      ["L"] .* reverse(layers)),
                             ticklabelsize = 14)
                # linkxaxes!(ax, ax2)
            end
            addlabels!(sf, labelformat)
            wsave(plotdir("fig2",
                          "residual_power_supplement$filebase.pdf"),
                  sf)
            sf
        end

        begin # * Recalculate the residual power in each band. Still in decibels, but labelled by depth
            Sr = deepcopy(ustripall.(S))
            map(Sr, L) do s, l # Map over structures
                idxs = indexin(lookup(s, SessionID), lookup(l, SessionID))
                l = l[:, idxs] # Match sessions just in case
                for i in CartesianIndices(l)
                    s[:, i] .= (s[:, i]) .- (l[i].(freqs(s[:, i]))) # Units of power spectral density
                end
            end
        end

        begin # * Plot the total residual theta power across channels
            ax = Axis(f[2, 1:2][1, 1]; ylabel = "Cortical depth (%)",
                      #   xticks = WilkinsonTicks(4),
                      xtickformat = depthticks,
                      ytickformat = depthticks,
                      xlabel = "Residual θ power (%)",
                      yreversed = true,
                      title = "Residual θ power") # [$(unit(eltype(S[1][1])))]

            θr = map(structures) do s
                ss = Sr[Structure = At(s)][Freq(theta)]
                # ss = upsample(ss, 5, 2)
                # no = mean(sum(ustripall(S[Structure = At(s)]), dims = Freq);
                #           dims = Depth)
                no = sum(ustripall(S[Structure = At(s)]), dims = Freq) # Total power of each channel
                x = sum(ss, dims = Freq) ./ no # The fraction of power above the 1/𝑓 component in a given frequency band
                x = dropdims(x, dims = Freq)
            end
            θr = ToolsArray(θr, (Structure(structures),))
            for (i, s) in reverse(collect(enumerate(structures)))
                x = θr[Structure = At(s)]
                μ = dropdims(mean(x, dims = SessionID), dims = SessionID)
                # σl = dropdims(quantile(x, 0.25, dims = SessionID), dims = SessionID)
                # σh = dropdims(quantile(x, 0.75, dims = SessionID), dims = SessionID)
                # σ = dropdims(std(x, dims = SessionID), dims = SessionID) ./ 2
                # σl = μ .- σ
                # σh = μ .+ σ
                μ, (σl, σh) = bootstrapmedian(x, dims = SessionID)
                μ, σl, σh = upsample.((μ, σl, σh), 5)

                band!(ax, Point2f.(collect(σl), lookup(μ, 1)),
                      Point2f.(collect(σh), lookup(μ, 1));
                      color = (structurecolors[i], 0.32), label = structures[i])
                lines!(ax, collect(μ), lookup(μ, 1); color = (structurecolors[i], alpha),
                       label = structures[i])
            end
            # leg = axislegend(ax, position = :lt, nbanks = 3, labelsize = 12, merge = true)
            # reverselegend!(leg)
            plotlayerints!(ax, layerints; axis = :y, newticks = false, flipside = true)
            ax.limits = ((-0.14, 0.55), (0, 1))
            display(f)
        end
        if is_flash # ?
            outdir = datadir("source_data", "fig2") # ?
            df = DataFrame() # ?
            for (i, s) in enumerate(structures) # ?
                sname = replace(s, "/" => "") # ?
                x = θr[Structure = At(s)] # ?
                μ, (σl, σh) = bootstrapmedian(x, dims = SessionID) # ?
                μ, σl, σh = upsample.((μ, σl, σh), 5) # ?
                if i == 1 # ?
                    df[!, "Cortical depth (%)"] = lookup(μ, 1) # ?
                end # ?
                df[!, "Median residual θ $sname (%)"] = collect(μ) # ?
                df[!, "CI lower $sname"] = collect(σl) # ?
                df[!, "CI upper $sname"] = collect(σh) # ?
            end # ?
            CSV.write(joinpath(outdir, "panel_e_residual_theta.csv"), df) # ?
        end # ?
        begin # * Does residual theta increase along layers
            tps = [SpatiotemporalMotifs.mediankendallpvalue(lookup(x, Depth), x)
                   for x in θr]
            τs = cat(first.(tps), last.(tps); dims = Var([:kendall, :pvalue]))
            mtau = median(τs[:, 1])
            open(statsfile, "a+") do file
                write(file, "\n## Residual θ\n")
                show(file, DataFrame(DimTable(τs; layersfrom = Var)))
                write(file, "\nMedian τ = $mtau\n")
            end
        end

        begin # * Residual gamma power across channels
            # f = Figure()
            ax = Axis(f[2, 1:2][1, 2]; #ylabel = "Cortical depth (%)",
                      #   xticks = WilkinsonTicks(4),
                      xtickformat = depthticks,
                      ytickformat = depthticks,
                      xlabel = "Residual γ power (%)",
                      yreversed = true,
                      title = "Residual γ power",
                      yticklabelsvisible = false) # [$(unit(eltype(S[1][1])))]

            γr = map(structures) do s
                ss = Sr[Structure = At(s)][Freq(gamma)]
                # no = mean(sum(ustripall(S[Structure = At(s)]), dims = Freq);
                #           dims = Depth)
                no = sum(ustripall(S[Structure = At(s)]), dims = Freq) # Total power of each channel
                ss = sum(ss, dims = Freq) ./ no # step(lookup(ss, Freq))
                ss = dropdims(ss, dims = Freq)
            end
            γr = ToolsArray(γr, (Structure(structures),))

            for (i, s) in structures |> enumerate |> collect |> reverse
                ss = γr[Structure = At(s)]
                μ, (σl, σh) = bootstrapmedian(ss, dims = SessionID)
                μ, σl, σh = upsample.((μ, σl, σh), 5)
                band!(ax, Point2f.(collect(σl), lookup(μ, 1)),
                      Point2f.(collect(σh), lookup(μ, 1));
                      color = (structurecolors[i], 0.32), label = structures[i])
                lines!(ax, collect(μ), lookup(μ, 1); color = (structurecolors[i], alpha),
                       label = structures[i])
            end

            # l = axislegend(ax, position = :rt, nbanks = 2, labelsize = 12, merge = true)
            # reverselegend!(l)
            plotlayerints!(ax, layerints; axis = :y, newticks = false, flipside = true)
            ax.limits = ((nothing, nothing), (0, 1))
            display(f)
        end
        if is_flash # ?
            outdir = datadir("source_data", "fig2") # ?
            df = DataFrame() # ?
            for (i, s) in enumerate(structures) # ?
                sname = replace(s, "/" => "") # ?
                ss = γr[Structure = At(s)] # ?
                μ, (σl, σh) = bootstrapmedian(ss, dims = SessionID) # ?
                μ, σl, σh = upsample.((μ, σl, σh), 5) # ?
                if i == 1 # ?
                    df[!, "Cortical depth (%)"] = lookup(μ, 1) # ?
                end # ?
                df[!, "Median residual γ $sname (%)"] = collect(μ) # ?
                df[!, "CI lower $sname"] = collect(σl) # ?
                df[!, "CI upper $sname"] = collect(σh) # ?
            end # ?
            CSV.write(joinpath(outdir, "panel_f_residual_gamma.csv"), df) # ?
        end # ?

        begin # * Plot the hierarchical correlation across layers
            N = 10000
            method = :group
            unidepths = commondepths(lookup.(χ, [Depth]))
            x = getindex.([SpatiotemporalMotifs.hierarchy_scores], structures)

            unichi = getindex.(χ, [Depth(Near(unidepths))])
            unichi = set.(unichi, [Depth => unidepths])
            y = stack(Structure(structures), unichi)

            unib = getindex.(b, [Depth(Near(unidepths))])
            unib = set.(unib, [Depth => unidepths])
            z = stack(Structure(structures), unib)

            @assert all(dims(y, Structure) .== structures)
            thet = getindex.(θr, [Depth(Near(unidepths))]) |> stack
            gamm = getindex.(γr, [Depth(Near(unidepths))]) |> stack

            μ, σ, 𝑝 = hierarchicalkendall(x, y, method; N)
            μb, σb, 𝑝b = hierarchicalkendall(x, z, method; N)
            μt, σt, 𝑝t = hierarchicalkendall(x, thet, method; N)
            μg, σg, 𝑝g = hierarchicalkendall(x, gamm, method; N)
        end
        begin
            # μ[𝑝 .> PTHR] .= NaN
            # μt[𝑝t .> PTHR] .= NaN
            # μg[𝑝g .> PTHR] .= NaN

            markersize = 10

            ax = Axis(f[3, 1:2][1, 3]; #ylabel = "Cortical depth (%)",
                      xlabel = "Kendall's 𝜏",
                      ytickformat = depthticks,
                      xtickformat,
                      title = "1/𝑓 gradients", limits = ((-0.73, 0.73), (0, 1)),
                      yreversed = true,
                      yticklabelsvisible = false)

            vlines!(ax, 0; color = :gray, linewidth = 3, linestyle = :dash)

            band!(ax, Point2f.(collect(first.(σ)), unidepths),
                  Point2f.(collect(last.(σ)), unidepths);
                  color = (crimson, bandalpha),
                  label = "1/𝑓 exponent")
            # lines!(ax, unidepths, collect(μ); alpha = bandalpha,
            #        label = "1/𝑓 exponent", color = crimson)
            scatter!(ax, collect(μ[𝑝 .< PTHR]), unidepths[𝑝 .< PTHR];
                     label = "1/𝑓 exponent", color = crimson, markersize)
            scatter!(ax, collect(μ[𝑝 .≥ PTHR]), unidepths[𝑝 .≥ PTHR]; color = :transparent,
                     strokecolor = crimson,
                     strokewidth = 1, markersize)

            band!(ax, Point2f.(collect(first.(σb)), unidepths),
                  Point2f.(collect(last.(σb)), unidepths);
                  color = (cornflowerblue, bandalpha),
                  label = "1/𝑓 intercept")
            # lines!(ax, unidepths, collect(μ); alpha = bandalpha,
            #        label = "1/𝑓 exponent", colorcornflowerblue)
            scatter!(ax, collect(μb[𝑝b .< PTHR]), unidepths[𝑝b .< PTHR];
                     label = "1/𝑓 intercept", color = cornflowerblue, markersize)
            scatter!(ax, collect(μb[𝑝b .≥ PTHR]), unidepths[𝑝b .≥ PTHR];
                     color = :transparent,
                     strokecolor = cornflowerblue,
                     strokewidth = 1, markersize)

            axislegend(ax, position = :lb, merge = true, labelsize = 12, nbanks = 1)

            plotlayerints!(ax, layerints; axis = :y, newticks = false, flipside = true)
        end
        if is_flash # ?
            outdir = datadir("source_data", "fig2") # ?
            df = DataFrame(; Depth = unidepths, Exponent_tau = collect(μ),
                           Exponent_CI_lower = collect(first.(σ)),
                           Exponent_CI_upper = collect(last.(σ)), Exponent_p = collect(𝑝),
                           Intercept_tau = collect(μb),
                           Intercept_CI_lower = collect(first.(σb)),
                           Intercept_CI_upper = collect(last.(σb)),
                           Intercept_p = collect(𝑝b)) # ?
            rename!(df, "Depth" => "Cortical depth (%)", "Exponent_tau" => "Exponent τ",
                    "Exponent_CI_lower" => "Exponent CI lower",
                    "Exponent_CI_upper" => "Exponent CI upper",
                    "Exponent_p" => "Exponent p", "Intercept_tau" => "Intercept τ",
                    "Intercept_CI_lower" => "Intercept CI lower",
                    "Intercept_CI_upper" => "Intercept CI upper",
                    "Intercept_p" => "Intercept p") # ?
            CSV.write(joinpath(outdir, "panel_g_1f_gradients.csv"), df) # ?
        end # ?
        begin
            # μ[𝑝 .> PTHR] .= NaN
            # μt[𝑝t .> PTHR] .= NaN
            # μg[𝑝g .> PTHR] .= NaN

            ax = Axis(f[2, 1:2][1, 3]; # ylabel = "Cortical depth (%)",
                      xlabel = "Kendall's 𝜏",
                      ytickformat = depthticks, xtickformat,
                      title = "Timescale gradients", limits = ((-0.73, 0.73), (0, 1)),
                      yreversed = true,
                      yticklabelsvisible = false)

            vlines!(ax, 0; color = :gray, linewidth = 3)

            band!(ax, Point2f.(collect(first.(σt)), unidepths),
                  Point2f.(collect(last.(σt)), unidepths);
                  color = (crimson, bandalpha), label = "Residual θ")
            # lines!(ax, unidepths, collect(μt); alpha = bandalpha, label = "Residual θ power",
            #        color = crimson)
            scatter!(ax, collect(μt[𝑝t .< PTHR]), unidepths[𝑝t .< PTHR];
                     label = "Residual θ", color = crimson, markersize)
            scatter!(ax, collect(μt[𝑝t .≥ PTHR]), unidepths[𝑝t .≥ PTHR];
                     color = :transparent, strokecolor = crimson,
                     strokewidth = 1, markersize)

            band!(ax, Point2f.(collect(first.(σg)), unidepths),
                  Point2f.(collect(last.(σg)), unidepths);
                  color = (cornflowerblue, bandalpha),
                  label = "Residual γ")
            # lines!(ax, unidepths, collect(μg); alpha = bandalpha, label = "Residual γ power",
            #    color = cornflowerblue)
            scatter!(ax, collect(μg[𝑝g .< PTHR]), unidepths[𝑝g .< PTHR];
                     label = "Residual γ", color = cornflowerblue, markersize)
            scatter!(ax, collect(μg[𝑝g .≥ PTHR]), unidepths[𝑝g .≥ PTHR];
                     color = :transparent, strokecolor = cornflowerblue,
                     strokewidth = 1, markersize)

            axislegend(ax, position = :lb, merge = true, labelsize = 12, nbanks = 1)

            plotlayerints!(ax, layerints; axis = :y, newticks = false, flipside = true)
        end
        if is_flash # ?
            outdir = datadir("source_data", "fig2") # ?
            df = DataFrame(; Depth = unidepths, Theta_tau = collect(μt),
                           Theta_CI_lower = collect(first.(σt)),
                           Theta_CI_upper = collect(last.(σt)), Theta_p = collect(𝑝t),
                           Gamma_tau = collect(μg), Gamma_CI_lower = collect(first.(σg)),
                           Gamma_CI_upper = collect(last.(σg)), Gamma_p = collect(𝑝g)) # ?
            rename!(df, "Depth" => "Cortical depth (%)", "Theta_tau" => "θ τ",
                    "Theta_CI_lower" => "θ CI lower", "Theta_CI_upper" => "θ CI upper",
                    "Theta_p" => "θ p", "Gamma_tau" => "γ τ",
                    "Gamma_CI_lower" => "γ CI lower", "Gamma_CI_upper" => "γ CI upper",
                    "Gamma_p" => "γ p") # ?
            CSV.write(joinpath(outdir, "panel_h_timescale_gradients.csv"), df) # ?
        end # ?
        addlabels!(f, ["a", "g", "c", "d", "e", "h", "f", "b"])
        f |> display
    end
    wsave(plotdir("fig2", "power_spectra$filebase.pdf"), f)
end
