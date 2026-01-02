#! /bin/bash
#=
exec julia -t auto --color=yes "${BASH_SOURCE[0]}" "$@"
=#
using DrWatson
@quickactivate "SpatiotemporalMotifs"

import DimensionalData: metadata
using SpatiotemporalMotifs
@preamble
import TimeseriesTools: freqs
using GraphMakie
using GraphMakie.Graphs
using SimpleWeightedGraphs
using LinearAlgebra
using Distributions
set_theme!(foresight(:physics))
Random.seed!(42)

begin # ! Don't change
    stimulus = r"Natural_Images" # "flash_250ms" # To do flash, you need to regenerate the plot data
    vars = [:ϕ, :ω]
end

if !isfile(savepath("fig4", Dict(), "jld2", calcdir("plots"))) # * Use extra workers if we can
    if SpatiotemporalMotifs.CLUSTER()
        using USydClusters
        USydClusters.Physics.addprocs(8; mem = 16, ncpus = 4,
                                      project = projectdir())
        @everywhere using SpatiotemporalMotifs
        @everywhere SpatiotemporalMotifs.@preamble
    end
end

plot_data, data_file = produce_or_load(Dict(), calcdir("plots");
                                       filename = savepath("fig4")) do config
    session_table = load(calcdir("plots", "posthoc_session_table.jld2"), "session_table")
    oursessions = session_table.ecephys_session_id

    path = calcdir("calculations")
    Q = calcquality(path)[Structure = At(structures)]
    out = load_calculations(Q; stimulus, vars)
    out = map(out) do O
        filter(o -> (o[:sessionid] in oursessions), O)
    end

    begin # * Functional hierarchy scores
        # ? See https://github.com/AllenInstitute/neuropixels_platform_paper/blob/master/Figure2/comparison_anatomical_functional_connectivity_final.ipynb
        dir = tempdir()
        baseurl = "https://github.com/AllenInstitute/neuropixels_platform_paper/raw/master/data/processed_data"
        f = "FFscore_grating_10.npy"
        Downloads.download(joinpath(baseurl, f), joinpath(dir, f))
        @assert SpatiotemporalMotifs.structures ==
                ["VISp", "VISl", "VISrl", "VISal", "VISpm", "VISam"] # The order of this FC matrix
        FF_score, FF_score_b, FF_score_std = eachslice(load(joinpath(dir, f)), dims = 1)
        FF_score = ToolsArray(FF_score', # Original matrix has higher-order structures on rows
                              (Dim{:structure_1}(SpatiotemporalMotifs.structures),
                               Dim{:structure_2}(SpatiotemporalMotifs.structures)))
    end

    begin # * Subject-by-subject phasedelays
        unidepths = range(0.05, 0.95, length = 19)
        pxy = map(out...) do o...
            ϕs = map(o) do _o
                ϕ = _o[:ϕ]
                ω = _o[:ω]
                odepths = lookup(ϕ, Depth)
                ϕ = ϕ[Depth(Near(unidepths))]
                ω = ω[Depth(Near(unidepths))]
                (ϕ, ω)
            end
            [getindex.(ϕs, i) for i in 1:2]
        end
        ϕs = getindex.(pxy, 1)
        ωs = getindex.(pxy, 2)
    end

    uniphis = map(ϕs, ωs) do ϕ, ω # Takes about 5 mins over 32 cores, 180 Gb
        changetimes = lookup.(ϕ, :changetime)
        latency = [maximum(abs.(a .- b))
                   for (a, b) in zip(changetimes[2:end], changetimes[1:(end - 1)])]
        @assert all(latency .< 0.05u"s")
        changetimes = first(changetimes)
        ts = times.(ϕ)
        @assert abs(minimum(first.(ts)) - maximum(first.(ts))) < 0.016u"s"
        ϕ = map(ϕ, ω) do x, w
            x[w .< 0u"Hz"] .= NaN # ! Mask negative frequency periods
            x[1:1562, :, :] # Some trials are a tiny bit shorter
        end
        ϕ = set.(ϕ, [𝑡 => times(ϕ[1])])
        ϕ = set.(ϕ, [Depth => Depth(unidepths)])
        ϕ = set.(ϕ, [Dim{:changetime} => Dim{:changetime}(changetimes)])
        stack(Structure(structures), ϕ)
    end
    Δϕ = map(uniphis) do uniphi
        Δ = [(a, b)
             for a in eachslice(uniphi, dims = 4), b in eachslice(uniphi, dims = 4)]
        Δ = Δ[filter(!=(0), triu(LinearIndices(Δ), 1))]
        Δ = map(Δ) do (a, b)
            mod.(b .- a .+ π, eltype(b)(2π)) .- π
        end
        Δ = stack(Dim{:pair}(eachindex(Δ)), Δ, dims = 4)
        Δ = permutedims(Δ, (4, 1, 2, 3))
    end
    plate = getindex.(Δϕ, [:], [1], [:], [1])

    begin # * Hierarchical distances
        @assert structures ==
                [DimensionalData.metadata(phi)[:structure] for phi in ϕs[1]]
        hs = getindex.([hierarchy_scores], structures) .|> Float32
        Δhs = [b - a for a in hs, b in hs] # Make a distance matrix
        Δhs = Δhs[filter(!=(0), triu(LinearIndices(Δhs), 1))]
        Δhs = Δhs ./ maximum(abs.(Δhs)) # Normalize to that the maximum distance between hierarchichal scores is 1
        Δhs = map(plate) do Δx # Copy onto shape of Δxs
            h = deepcopy(Δx)
            for _h in eachslice(h, dims = Depth)
                _h .= Δhs
            end
            return h
        end
        @assert maximum(maximum.(Δhs)) ≤ 1
    end

    begin # * Functional distances
        Δfs = FF_score[filter(!=(0), triu(LinearIndices(FF_score), 1))] # ? Δhs and Δfs match data used for siegle2021 figure 2f, without same-region FC. FF_score is already a 'distance' matrix.
        Δfs = Δfs ./ maximum(abs.(Δfs))
        Δfs = map(plate) do Δx
            f = deepcopy(Δx)
            for _f in eachslice(f, dims = Depth)
                _f .= Δfs
            end
            return f
        end
        @assert maximum(maximum.(Δfs)) ≤ 1
    end

    ∂h = map(Δϕ, Δhs) do Δ, Δh
        mapslices(Δ; dims = (:pair, Depth)) do Δ
            .-sign.(Δ) ./ (Δh) # Minus because phase increases over time
        end
    end
    ∂f = map(Δϕ, Δfs) do Δ, Δf
        mapslices(Δ; dims = (:pair, Depth)) do Δ
            .-sign.(Δ) ./ (Δf) # Minus because phase increases over time
        end
    end

    ∂h = map(∂h) do ∂h
        dropdims(nansafe(mean, dims = :pair)(∂h), dims = :pair)
    end
    ∂f = map(∂f) do ∂f
        dropdims(nansafe(mean, dims = :pair)(∂f), dims = :pair)
    end

    begin # * Average quantities
        ∂h̄ = dropdims(mean(nansafe(mean, dims = :changetime).(∂h)); dims = :changetime)
        ∂f̄ = dropdims(mean(nansafe(mean, dims = :changetime).(∂f)); dims = :changetime)
    end

    begin # * Probe-shuffled nulls
        Δϕ_sur = map(uniphis) do uniphi
            idxs = randperm(size(uniphi, 4)) # Random probe reshuffling
            uniphi = uniphi[:, :, :, idxs]
            Δ = [(a, b)
                 for a in eachslice(uniphi, dims = 4),
                     b in eachslice(uniphi, dims = 4)]
            Δ = Δ[filter(!=(0), triu(LinearIndices(Δ), 1))]
            Δ = map(Δ) do (a, b)
                mod.(b .- a .+ π, eltype(b)(2π)) .- π
            end
            Δ = stack(Dim{:pair}(eachindex(Δ)), Δ, dims = 4)
            Δ = permutedims(Δ, (4, 1, 2, 3))
        end
        ∂h_sur = map(Δϕ_sur, Δhs) do Δ, Δh
            mapslices(Δ; dims = (:pair, Depth)) do Δ
                .-sign.(Δ) ./ (Δh) # Minus because phase increases over time
            end
        end
        ∂f_sur = map(Δϕ_sur, Δfs) do Δ, Δf
            mapslices(Δ; dims = (:pair, Depth)) do Δ
                .-sign.(Δ) ./ (Δf) # Minus because phase increases over time
            end
        end
        ∂h_sur = map(∂h_sur) do ∂h
            dropdims(nansafe(mean, dims = :pair)(∂h), dims = :pair)
        end
        ∂f_sur = map(∂f_sur) do ∂f
            dropdims(nansafe(mean, dims = :pair)(∂f), dims = :pair)
        end
        ∂h̄_sur = dropdims(mean(nansafe(mean, dims = :changetime).(∂h_sur));
                           dims = :changetime)
        ∂f̄_sur = dropdims(mean(nansafe(mean, dims = :changetime).(∂f_sur));
                           dims = :changetime)
    end

    begin # * Calculate p-values of for the averages by comparing to the surrogate distribution
        # ? We want to compare the distribution of real order parameters to the shuffled
        # distribution; we have one shuffled value for each real value, so we can use the
        # paired Wilcoxon test on the non-NaN values.
        begin # * anatomical
            𝑝_h = deepcopy(∂h̄)
            # super = cat(∂h...; dims = 3)
            # super_sur = cat(parent.(∂h_sur)...; dims = 3)
            super = mean.(∂h; dims = 3)
            super_sur = mean.(∂h_sur; dims = 3)
            super = cat(super..., dims = 3)
            super_sur = cat(super_sur..., dims = 3)
            𝑝_h .= map(eachslice(super, dims = (1, 2)),
                       eachslice(super_sur, dims = (1, 2))) do x, y
                idxs = (.!isnan.(x)) .& (.!isnan.(y))
                # 𝑝 = pvalue(HypothesisTests.MannWhitneyUTest(x[idxs], y[idxs]))
                _x = Float64.(x[idxs]) |> collect
                _y = Float64.(y[idxs]) |> collect
                𝑝 = pvalue(HypothesisTests.SignedRankTest(_x, _y))
            end
        end

        begin # * functional
            𝑝_f = deepcopy(∂f̄)
            # super = cat(∂f...; dims = 3)
            # super_sur = cat(parent.(∂f_sur)...; dims = 3)
            super = mean.(∂f; dims = 3)
            super_sur = mean.(∂f_sur; dims = 3)
            super = cat(super..., dims = 3)
            super_sur = cat(super_sur..., dims = 3)
            𝑝_f .= map(eachslice(super, dims = (1, 2)),
                       eachslice(super_sur, dims = (1, 2))) do x, y
                idxs = (.!isnan.(x)) .& (.!isnan.(y))
                # 𝑝 = pvalue(HypothesisTests.MannWhitneyUTest(x[idxs], y[idxs]))
                _x = Float64.(x[idxs]) |> collect
                _y = Float64.(y[idxs]) |> collect
                𝑝 = pvalue(HypothesisTests.SignedRankTest(_x, _y))
            end
        end
    end

    return (@strdict unidepths FF_score ∂h̄ ∂f̄ ∂h̄_sur ∂f̄_sur 𝑝_h 𝑝_f)
end

begin # * Plots
    @info "Plotting inter-areal phase delays"
    @unpack unidepths, FF_score, ∂h̄, ∂f̄, ∂h̄_sur, ∂f̄_sur, 𝑝_h, 𝑝_f = plot_data
    fullplot = (SpatiotemporalMotifs.THETA() == (3, 10)) &&
               (SpatiotemporalMotifs.GAMMA() == (30, 100)) &&
               !SpatiotemporalMotifs.CAUSAL_FILTER() &&
               stimulus == r"Natural_Images"

    if fullplot
        f = Foresight._panels(; size = (720, 460))
        gs = subdivide(f, 2, 2)
    else
        f = Foresight.TwoPanel()
        gs = subdivide(f, 1, 2)
        gs = permutedims([gs gs], (2, 1)) # So indexing is the same
    end

    layerints = load(calcdir("plots", "grand_unified_layers.jld2"), "layerints")
    if fullplot # * Schematic of hierarchy
        ax = Axis(gs[1], yreversed = true, aspect = DataAspect(),
                  title = "Anatomical hierarchy")
        ag = getindex.([hierarchy_scores], structures)
        ag = [b - a for (a, b) in Iterators.product(ag, ag)]
        ag[ag .< 0.0] .= NaN
        ag = 1.0 .- abs.(ag)
        ag[isnan.(ag)] .= 0.0
        ag[ag .< 0.6] .= 0.0
        ag = ag .- Diagonal(ag)
        ag[ag .> 0] = MinMax(ag[ag .> 0])(ag[ag .> 0])
        ag = SimpleWeightedDiGraph(ag .* 3)

        # fg = collect(FF_score[At(structures), At(structures)])
        # fg = fg ./ maximum(abs.(fg))
        # fg = (fg - fg') / 2
        # fg[fg .< 0.0] .= NaN
        # fg = 1.0 .- abs.(fg)
        # fg[isnan.(fg)] .= 0.0
        # fg[fg .< 0.3] .= 0.0
        # fg = fg .- Diagonal(fg)
        # fg[fg .> 0] = MinMax(fg[fg .> 0])(fg[fg .> 0])
        # fg = SimpleWeightedDiGraph(fg .* 3)

        hidedecorations!(ax)
        hidespines!(ax)
        outline, infill, l = plot_visual_cortex!(ax; colormap = structurecolors,
                                                 fillalpha = 0.2,
                                                 strokealpha = 0.8)
        delete!(l)
        plotstructurecenters!(ax, ag; curve_distance = -20, node_size = 25,
                              ilabels_fontsize = 7)
    end
    if fullplot # * Schematic of hierarchy
        ax = Axis(gs[3], yreversed = true, aspect = DataAspect(),
                  title = "Functional hierarchy")

        fg = collect(FF_score[At(structures), At(structures)])
        fg = fg ./ maximum(abs.(fg))
        fg = (fg - fg') / 2
        fg[fg .< 0.0] .= NaN
        fg = 1.0 .- abs.(fg)
        fg[isnan.(fg)] .= 0.0
        fg[fg .< 0.3] .= 0.0
        fg = fg .- Diagonal(fg)
        fg[fg .> 0] = MinMax(fg[fg .> 0])(fg[fg .> 0])
        fg = SimpleWeightedDiGraph(fg .* 3)

        hidedecorations!(ax)
        hidespines!(ax)
        outline, infill, l = plot_visual_cortex!(ax; colormap = structurecolors,
                                                 fillalpha = 0.2,
                                                 strokealpha = 0.8)
        delete!(l)
        plotstructurecenters!(ax, fg; edge_color = :crimson, node_size = 25,
                              ilabels_fontsize = 7)
    end
    begin # * Correlation to hierarchy score
        ax = Axis(gs[2][1, 1], yreversed = true,
                  ylabel = "Cortical depth (%)", title = " ", ytickformat = depthticks,
                  xticks = -0.25:0.25:0.75, xtickformat = terseticks)
        if !fullplot
            ax.xlabel = "Time (s)"
        end
        plevels = [-3.0, -5.0]
        pdiff = mean(diff(plevels))
        prange = extrema(plevels) .+ [pdiff, -pdiff] ./ 2
        prange = Tuple(prange)
        levelmap = cgrad(:binary, range(0, 1, length = length(plevels) + 1);
                         categorical = true)
        # ∂̄ = dropdims(mean(∂h, dims = Trial), dims = Trial)
        H = deepcopy(𝑝_h)
        H[:] .= adjust(H[:], BenjaminiHochberg())
        H = log10.(H)
        p, _ = plotlayermap!(ax, ∂h̄[𝑡(SpatiotemporalMotifs.INTERVAL)] |> ustripall,
                             colormap = binarysunset,
                             colorrange = (-1.0, 1.0))
        contour!(ax, H[𝑡(SpatiotemporalMotifs.INTERVAL)] |> ustripall,
                 colormap = levelmap, levels = plevels, linewidth = 1.5,
                 linestyle = :dash, colorrange = prange)
        Colorbar(gs[2][1, 2], p;
                 label = rich("Mean order parameter ",
                              rich("A", subscript("θ"), font = "Times Italic")),
                 tickformat = terseticks)
        Colorbar(gs[2][1, 1]; colormap = levelmap,
                 ticks = plevels,
                 colorrange = extrema(plevels) .+
                              [mean(diff(plevels)), -mean(diff(plevels))] ./ 2,
                 tickformat = X -> [L"10^{%$(round(Int, x))}" for x in X],
                 vertical = false,
                 flipaxis = true, label = "Corrected 𝑝-value", tellheight = false,
                 valign = :top)
        plotlayerints!(ax, layerints; flipside = false, newticks = false,
                       bgcolor = Makie.RGBA(0, 0, 0, 0))
        ax.limits = (nothing, (0.05, 0.95))
        f
    end
    begin # * Correlation to functional hierarchy score
        ax = Axis(gs[4][1, 1], yreversed = true, xlabel = "Time (s)",
                  ylabel = "Cortical depth (%)", title = " ", ytickformat = depthticks,
                  xticks = -0.25:0.25:0.75, xtickformat = terseticks)
        if !fullplot
            ax.ylabel = ""
            ax.yticklabelsvisible = false
        end
        # ∂̄ = dropdims(mean(∂h, dims = Trial), dims = Trial)
        H = deepcopy(𝑝_f)
        H[:] .= adjust(H[:], BenjaminiHochberg())
        H = log10.(H)
        p, _ = plotlayermap!(ax, ∂f̄[𝑡(SpatiotemporalMotifs.INTERVAL)] |> ustripall,
                             colormap = binarysunset,
                             colorrange = (-0.74, 0.74))
        contour!(ax, H[𝑡(SpatiotemporalMotifs.INTERVAL)] |> ustripall,
                 colormap = levelmap, levels = [0, plevels...], linewidth = 1.5,
                 linestyle = :dash, colorrange = prange)
        Colorbar(gs[4][1, 2], p;
                 label = rich("Mean order parameter ",
                              rich("F", subscript("θ"), font = "Times Italic")),
                 tickformat = terseticks)
        if fullplot
            clabel = ""
        else
            clabel = "Corrected 𝑝-value"
        end
        Colorbar(gs[4][1, 1]; colormap = levelmap, ticks = plevels,
                 colorrange = extrema(plevels) .+
                              [mean(diff(plevels)), -mean(diff(plevels))] ./ 2,
                 tickformat = X -> [L"10^{%$(round(Int, x))}" for x in X],
                 vertical = false,
                 flipaxis = true,
                 tellheight = false,
                 valign = :top,
                 label = clabel)
        plotlayerints!(ax, layerints; flipside = false, newticks = false,
                       bgcolor = Makie.RGBA(0, 0, 0, 0))
        ax.limits = (nothing, (0.05, 0.95))
        f
    end
    if fullplot
        colsize!(f.layout, 1, Relative(0.35))
        addlabels!(f, labelformat)
    else
        colgap!(f.layout, Relative(0.075))
        if stimulus === r"Natural_Images"
            addlabels!(f, ["h", "i"])
        else
            addlabels!(f, ["a", "b"])
        end
    end
    wsave(plotdir("fig4", "interareal_phasedelays.pdf"), f)
    if fullplot # ?
        begin # ? Save source data
            mkpath(datadir("source_data", "fig4")) # ?
            outdir = datadir("source_data", "fig4") # ?
            # ? Anatomical hierarchy order parameter
            data_h = ∂h̄[𝑡(SpatiotemporalMotifs.INTERVAL)] |> ustripall # ?
            H_h = deepcopy(𝑝_h) # ?
            H_h[:] .= adjust(H_h[:], BenjaminiHochberg()) # ?
            df_h = DataFrame(data_h) # ?
            df_h_p = DataFrame(H_h[𝑡(SpatiotemporalMotifs.INTERVAL)]) # ?
            rename!(df_h, "𝑡" => "Time (s)", "Depth" => "Cortical depth (%)", # ?
                    last(names(df_h)) => "A_θ (anatomical order parameter)") # ?
            df_h[!, "Corrected p-value"] = df_h_p[!, last(names(df_h_p))] # ?
            CSV.write(joinpath(outdir, "panel_b_anatomical_order_parameter.csv"), df_h) # ?
            # ? Functional hierarchy order parameter
            data_f = ∂f̄[𝑡(SpatiotemporalMotifs.INTERVAL)] |> ustripall # ?
            H_f = deepcopy(𝑝_f) # ?
            H_f[:] .= adjust(H_f[:], BenjaminiHochberg()) # ?
            df_f = DataFrame(data_f) # ?
            df_f_p = DataFrame(H_f[𝑡(SpatiotemporalMotifs.INTERVAL)]) # ?
            rename!(df_f, "𝑡" => "Time (s)", "Depth" => "Cortical depth (%)", # ?
                    last(names(df_f)) => "F_θ (functional order parameter)") # ?
            df_f[!, "Corrected p-value"] = df_f_p[!, last(names(df_f_p))] # ?
            CSV.write(joinpath(outdir, "panel_d_functional_order_parameter.csv"), df_f) # ?
        end # ?
    elseif stimulus == "flash_250ms" # ?
        begin # ? Save source data for figS16 (flash stimulus hierarchical order parameters)
            mkpath(datadir("source_data", "figS16")) # ?
            outdir = datadir("source_data", "figS16") # ?
            # ? Anatomical hierarchy order parameter
            data_h = ∂h̄[𝑡(SpatiotemporalMotifs.INTERVAL)] |> ustripall # ?
            H_h = deepcopy(𝑝_h) # ?
            H_h[:] .= adjust(H_h[:], BenjaminiHochberg()) # ?
            df_h = DataFrame(data_h) # ?
            df_h_p = DataFrame(H_h[𝑡(SpatiotemporalMotifs.INTERVAL)]) # ?
            rename!(df_h, "𝑡" => "Time (s)", "Depth" => "Cortical depth (%)", # ?
                    last(names(df_h)) => "A_θ (anatomical order parameter)") # ?
            df_h[!, "Corrected p-value"] = df_h_p[!, last(names(df_h_p))] # ?
            CSV.write(joinpath(outdir, "panel_a_anatomical_order_parameter.csv"), df_h) # ?
            # ? Functional hierarchy order parameter
            data_f = ∂f̄[𝑡(SpatiotemporalMotifs.INTERVAL)] |> ustripall # ?
            H_f = deepcopy(𝑝_f) # ?
            H_f[:] .= adjust(H_f[:], BenjaminiHochberg()) # ?
            df_f = DataFrame(data_f) # ?
            df_f_p = DataFrame(H_f[𝑡(SpatiotemporalMotifs.INTERVAL)]) # ?
            rename!(df_f, "𝑡" => "Time (s)", "Depth" => "Cortical depth (%)", # ?
                    last(names(df_f)) => "F_θ (functional order parameter)") # ?
            df_f[!, "Corrected p-value"] = df_f_p[!, last(names(df_f_p))] # ?
            CSV.write(joinpath(outdir, "panel_b_functional_order_parameter.csv"), df_f) # ?
        end # ?
    end # ?
    f
end
