#! /bin/bash
# -*- mode: julia -*-
#=
exec julia -t auto --color=yes "${BASH_SOURCE[0]}" "$@"
=#
using DrWatson
@quickactivate "SpatiotemporalMotifs"

using SpatiotemporalMotifs
import TimeseriesTools: freqs
using Peaks
using FileIO
using SpatiotemporalMotifs
import SpatiotemporalMotifs.layers
@preamble
set_theme!(foresight(:physics))

stimulus = "spontaneous" # "flash_250ms"
xtickformat = x -> string.(round.(Int, x))
theta = 5 .. 10
gamma = 30 .. 90
alpha = 0.8
bandalpha = 0.2

session_table = load(datadir("posthoc_session_table.jld2"), "session_table")
oursessions = session_table.ecephys_session_id

path = datadir("power_spectra")
Q = calcquality(path)[stimulus = At(stimulus), structure = At(structures)]
Q = Q[Dim{:sessionid}(At(oursessions))]
f = Figure(size = (900, 1080))

begin # * Load data
    S = map(lookup(Q, :structure)) do structure
        out = map(lookup(Q, :sessionid)) do sessionid
            if Q[sessionid = At(sessionid), structure = At(structure)] == 0
                return nothing
            end
            filename = savepath((@strdict sessionid structure stimulus), "jld2", path)
            S = load(filename, "S")
        end
        out = filter(!isnothing, out)
        out = filter(x -> maximum(DimensionalData.metadata(x)[:streamlinedepths]) > 0.90,
                     out) # Remove sessions that don't have data to a reasonable depth

        m = DimensionalData.metadata.(out)
        sessions = getindex.(m, :sessionid)

        streamlinedepths = getindex.(m, :streamlinedepths)
        layerinfo = getindex.(m, :layerinfo)

        unidepths = commondepths(streamlinedepths)
        out = map(out, streamlinedepths, layerinfo) do o, s, l
            o = set(o, Dim{:depth} => s)
            layernames = DimArray(l[1], Dim{:depth}(lookup(o, Dim{:depth})))
            layernums = DimArray(l[3], Dim{:depth}(lookup(o, Dim{:depth})))
            o = o[Dim{:depth}(Near(unidepths))]
            layernames = layernames[Dim{:depth}(Near(unidepths))]
            layernums = layernums[Dim{:depth}(Near(unidepths))]
            # @assert length(unique(lookup(o, Dim{:depth}))) == length(unidepths)
            @assert issorted(lookup(o, Dim{:depth}))
            push!(o.metadata, :layernames => layernames)
            push!(o.metadata, :layernums => layernums)
            o = set(o, Dim{:depth} => unidepths)
        end
        layernames = DimArray(stack(getindex.(DimensionalData.metadata.(out), :layernames)),
                              (Dim{:depths}(unidepths), Dim{:sessionid}(sessions)))
        layernums = DimArray(stack(getindex.(DimensionalData.metadata.(out), :layernums)),
                             (Dim{:depths}(unidepths), Dim{:sessionid}(sessions)))
        S = stack(Dim{:sessionid}(sessions), out, dims = 3)
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
        # s = set(s, Dim{:depth} => Dim{:layer}(layernum2name.(parent(l)[:])))
        s = set(s, Dim{:depth} => Dim{:layer}(parent(l)[:]))
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
    S̄ = DimArray(S̄ |> collect, (Dim{:structure}(lookup(Q, :structure)),))
    S = DimArray(S |> collect, (Dim{:structure}(lookup(Q, :structure)),))
    S̄ = map(S̄) do s
        N = UnitEnergy(s, dims = 1)
        N(s)
    end
    layerints = load(datadir("grand_unified_layers.jld2"), "layerints")
end;

begin # * Mean power spectrum in VISp and VISam. Bands show 1 S.D.
    ax = Axis(f[1, 1]; xscale = log10, yscale = log10,
              limits = ((3, 300), (1e-5, 1)),
              xlabel = "Frequency (Hz)",
              xgridvisible = true,
              ygridvisible = true,
              xgridstyle = :dash,
              ygridstyle = :dash,
              xtickformat,
              ylabel = "Mean power spectral density (a.u.)",
              title = "Power spectral density")

    # * Band annotations
    vspan!(ax, extrema(theta)..., color = (crimson, 0.22),
           label = "𝛉 ($(theta.left) – $(theta.right) Hz)")
    vlines!(ax, [7.0], color = (crimson, 0.42), linestyle = :dash, linewidth = 4)
    vspan!(ax, extrema(gamma)..., color = (cornflowerblue, 0.22),
           label = "𝛄 ($(gamma.left) – $(gamma.right) Hz)")
    vlines!(ax, [54.4], color = (cornflowerblue, 0.42), linestyle = :dash, linewidth = 4)

    psa = map(enumerate(structures)) do (i, s)
        s = S̄[structure = At(s)][Freq = 3u"Hz" .. 300u"Hz"]
        s = s ./ 10^((i - 1.5) / 2.5)
        yc = only(mean(s[Freq = Near(3u"Hz")]))
        if i == 1
            plotspectrum!(ax, s; textposition = (3, yc),
                          color = structurecolors[i], annotations = [:peaks])
        else
            plotspectrum!(ax, s; textposition = (14, yc),
                          color = structurecolors[i], annotations = [])
        end
    end
    ps, α = first.(psa), last.(psa)
    α = round.(α, sigdigits = 3)
    axislegend(ax, position = :lb, labelsize = 12, backgroundcolor = :white,
               framevisible = true, padding = (5, 5, 5, 5))
    leg = ["$s (α = $α)" for (s, α) in zip(structures, α)]
    map(enumerate(leg)) do (i, s)
        text!(ax, 290, exp10(-1.1 - (i / 3 - 1 - 0.1)); text = s,
              color = structurecolors[i],
              fontsize = 14,
              align = (:right, :bottom))
    end
    # axislegend(ax, ps, leg, padding = 5, framevisible = true, labelsize = 12)
    f
end

begin # * Calculate the channel-wise fits
    file = datadir("fooof.jld2")
    if isfile(file)
        χ, L = load(file, "χ", "L")
    else
        L = map(S) do s
            map(fooof, eachslice(ustripall(s), dims = (2, 3)))
        end
        χ = [getindex.(last.(l), :χ) for l in L]
        L = [first.(l) for l in L]
        save(file, Dict("χ" => χ, "L" => L))
    end
    L = getindex.(L, [Dim{:sessionid}(At(oursessions))])
    L = L[structure = At(structures)]
    @assert size(L[1])[end] == size(S[1])[end]
end

begin # * Plot the exponent
    # f = Figure()
    ax = Axis(f[1, 2]; xlabel = "Cortical depth (%)", ylabel = "1/f exponent",
              limits = ((0, 1), (0.9, 2.1)), xtickformat = depthticks,
              title = "1/f exponent")
    for (i, chi) in χ |> enumerate |> collect |> reverse
        μ, (σl, σh) = bootstrapmedian(chi, dims = :sessionid)
        μ, σl, σh = upsample.((μ, σl, σh), 5)

        band!(ax, lookup(μ, 1), collect(σl), collect(σh);
              color = (structurecolors[i], 0.32), label = structures[i])
        lines!(ax, lookup(μ, 1), μ; color = (structurecolors[i], alpha),
               label = structures[i])
    end
    l = axislegend(ax, position = :rb, nbanks = 2, labelsize = 12, merge = true)
    reverselegend!(l)
    plotlayerints!(ax, layerints; axis = :x, newticks = false, flipside = true)
    display(f)
end

begin # * Relative power (supplement)
    fs = SixPanel()
    gs = subdivide(fs, 3, 2)
    map(enumerate(S)) do (i, x)
        x = mapslices(x, dims = (1, 2)) do s
            s ./ maximum(s, dims = 2)
        end |> ustripall
        x = dropdims(median(x[Freq = 1 .. 300], dims = :sessionid), dims = :sessionid)[:,
                                                                                       2:end]

        ax = Axis(gs[i][1, 1], xlabel = "Frequency (Hz)", ylabel = "Cortical depth (%)",
                  xtickformat = depthticks,
                  limits = ((1, 200), (nothing, nothing)), yreversed = true, aspect = 1,
                  title = metadata(x)[:structure], xticks = [1, 100, 200])
        x = upsample(x, 5, 2)
        p = heatmap!(ax, lookup(x, 1), lookup(x, 2) .* 100, x |> collect,
                     colormap = :viridis, colorrange = (0, 1), rasterize = 10)
        if i % 2 == 0
            Colorbar(gs[i][1, 2], p, label = "Relative power (a.u.)", tellheight = true)
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
    save(plotdir("power_spectra", "relative_power.pdf"), fs)
end

begin # * Plot fooof residuals. Bands are 1 S.D.
    # f = Figure()
    Sr_log = map(ustripall.(S), L, meanlayers) do s, l, m
        s = deepcopy(s)
        map(eachslice(s, dims = (:depth, :sessionid)), l) do s, l
            _s = log10.(ustripall(s))
            s .= _s .- (_s |> freqs .|> l .|> log10)
        end
        s = set(s, Dim{:depth} => Dim{:layer}(layernum2name.(parent(m)[:])))
        s = set(s, :layer => DimensionalData.Irregular)
    end

    Sr_log = DimArray(Sr_log |> collect,
                      (Dim{:structure}(lookup(Q, :structure)),))

    for (i, structure) in enumerate(["VISp", "VISl"])
        ax2 = Axis(f[2, i]; xscale = log10,
                   limits = ((3, 300), (-0.1, 3.5)), xtickformat, xlabel = "Frequency (Hz)",
                   ylabel = "Residual spectral density (dB)",
                   title = "Residual spectral density in " * structure) # xticksvisible = false, yaxisposition = :right,
        #    xticklabelsvisible = false,
        vlines!(ax2, [7.0], color = (crimson, 0.42), linestyle = :dash, linewidth = 4)
        vlines!(ax2, [54.4], color = (cornflowerblue, 0.42), linestyle = :dash,
                linewidth = 4)

        for (i, (c, l)) in (reverse ∘ collect ∘ enumerate ∘ zip)(layercolors, layers)
            s = Sr_log[structure = At(structure)][Freq(3 .. 300)]
            s = s[layer = (lookup(s, :layer) .== [l])]
            s = dropdims(mean(s, dims = :layer), dims = :layer)
            d = (length(layers) - i + 1) / 2
            hlines!(ax2, [d]; color = (c, 0.22), linestyle = :dash)
            μ = dropdims(mean(s, dims = :sessionid), dims = :sessionid) .+
                d
            σ = dropdims(std(s, dims = :sessionid), dims = :sessionid)
            band!(ax2, TimeseriesTools.freqs(μ), collect(μ .- σ), collect(μ .+ σ);
                  color = (c, bandalpha))
            lines!(ax2, TimeseriesTools.freqs(μ), μ; color = (c, alpha))
        end
        C = Colorbar(f[2, i][1, 2],
                     colormap = reverse(cgrad(layercolors, categorical = true)),
                     ticks = (range(0, 1, length = (2 * length(layercolors) + 1))[2:2:end],
                              [structure] .* reverse(layers)), ticklabelrotation = π / 2,
                     ticklabelsize = 13)
        # linkxaxes!(ax, ax2)
    end
    f
end

begin # * Residual power supplement
    sf = SixPanel()
    gs = subdivide(sf, 3, 2)

    for (i, structure) in enumerate(structures)
        ax2 = Axis(gs[i]; xscale = log10,
                   limits = ((3, 300), (-0.1, 3.5)), xtickformat, xlabel = "Frequency (Hz)",
                   ylabel = "Residual power (dB)", title = structure) # xticksvisible = false, yaxisposition = :right,
        #    xticklabelsvisible = false,
        vlines!(ax2, [7.0], color = (crimson, 0.42), linestyle = :dash, linewidth = 4)
        vlines!(ax2, [54.4], color = (cornflowerblue, 0.42), linestyle = :dash,
                linewidth = 4)

        for (i, (c, l)) in (reverse ∘ collect ∘ enumerate ∘ zip)(layercolors, layers)
            s = Sr_log[structure = At(structure)][Freq(3 .. 300)]
            s = s[layer = (lookup(s, :layer) .== [l])]
            s = dropdims(mean(s, dims = :layer), dims = :layer)
            d = (length(layers) - i + 1) / 2
            hlines!(ax2, [d]; color = (c, 0.22), linestyle = :dash)
            μ = dropdims(mean(s, dims = :sessionid), dims = :sessionid) .+
                d
            σ = dropdims(std(s, dims = :sessionid), dims = :sessionid)
            band!(ax2, TimeseriesTools.freqs(μ), collect(μ .- σ), collect(μ .+ σ);
                  color = (c, bandalpha))
            lines!(ax2, TimeseriesTools.freqs(μ), μ; color = (c, alpha))
        end
        C = Colorbar(gs[i][1, 2],
                     colormap = reverse(cgrad(layercolors, categorical = true)),
                     ticks = (range(0, 1, length = (2 * length(layercolors) + 1))[2:2:end],
                              ["L"] .* reverse(layers)),
                     ticklabelsize = 14)
        # linkxaxes!(ax, ax2)
    end
    addlabels!(sf)
    wsave(plotdir("power_spectra", "residual_power_supplement.pdf"), sf)
    sf
end

begin # * Calculate the residual power in each band
    Sr = deepcopy(ustripall.(S))
    map(Sr, L) do s, l
        for i in CartesianIndices(l)
            s[:, i] .= (s[:, i]) .- (l[i].(freqs(s[:, i])))
        end
    end
end

begin # * Plot the total residual theta power across channels
    # f = Figure()
    ax = Axis(f[3, 1]; xlabel = "Cortical depth (%)", yticks = WilkinsonTicks(4),
              xtickformat = depthticks,
              ylabel = "Residual 𝜽 ($(theta.left) – $(theta.right) Hz) [$(unit(eltype(S[1][1])))]",
              yticklabelrotation = π / 2,
              title = "Residual θ power")

    for (i, s) in reverse(collect(enumerate(structures)))
        ss = Sr[structure = At(s)][Freq(theta)]
        # ss = upsample(ss, 5, 2)
        x = sum(ss, dims = Freq) ./ step(lookup(ss, 1))
        x = dropdims(x, dims = Freq)
        μ = dropdims(mean(x, dims = :sessionid), dims = :sessionid)
        # σl = dropdims(quantile(x, 0.25, dims = :sessionid), dims = :sessionid)
        # σh = dropdims(quantile(x, 0.75, dims = :sessionid), dims = :sessionid)
        # σ = dropdims(std(x, dims = :sessionid), dims = :sessionid) ./ 2
        # σl = μ .- σ
        # σh = μ .+ σ
        μ, (σl, σh) = bootstrapmedian(x, dims = :sessionid)
        μ, σl, σh = upsample.((μ, σl, σh), 5)

        band!(ax, lookup(μ, 1), collect(σl), collect(σh);
              color = (structurecolors[i], 0.32), label = structures[i])
        lines!(ax, lookup(μ, 1), collect(μ);
               color = (structurecolors[i], alpha), label = structures[i])
    end
    leg = axislegend(ax, position = :lt, nbanks = 2, labelsize = 12, merge = true)
    reverselegend!(leg)
    plotlayerints!(ax, layerints; axis = :x, newticks = false, flipside = false)
    ax.limits = ((0, 1), (nothing, nothing))
    display(f)
end

begin # * Residual gamma power across channels
    # f = Figure()
    ax = Axis(f[3, 2]; xlabel = "Cortical depth (%)", yticks = WilkinsonTicks(4),
              xtickformat = depthticks,
              ylabel = "Residual 𝜸 ($(gamma.left) – $(gamma.right) Hz) [$(unit(eltype(S[1][1])))]",
              yticklabelrotation = π / 2,
              title = "Residual γ power")

    for (i, s) in structures |> enumerate |> collect |> reverse
        ss = Sr[structure = At(s)][Freq(gamma)]
        ss = sum(ss, dims = Freq) ./ step(lookup(ss, Freq))
        ss = dropdims(ss, dims = Freq)
        μ, (σl, σh) = bootstrapmedian(ss, dims = :sessionid)
        μ, σl, σh = upsample.((μ, σl, σh), 5)
        band!(ax, lookup(μ, 1), collect(σl), collect(σh);
              color = (structurecolors[i], bandalpha), label = structures[i])

        lines!(ax, lookup(μ, 1), μ;
               color = (structurecolors[i], alpha), label = structures[i])
    end

    l = axislegend(ax, position = :rt, nbanks = 2, labelsize = 12, merge = true)
    reverselegend!(l)
    plotlayerints!(ax, layerints; axis = :x, newticks = false, flipside = false)
    ax.limits = ((0, 1), (nothing, nothing))
    display(f)
end

# begin # * Plot the spectral width of the gamma band
#     # f = Figure()
#     ax = Axis(f[2, 2]; xlabel = "Cortical depth (%)", yticks = WilkinsonTicks(4),
#               xtickformat = depthticks,
#               ylabel = "𝜸 spectral width ($(gamma.left) – $(gamma.right) Hz) [Hz]",
#               yticklabelrotation = π / 2)

#     for (i, s) in (reverse ∘ collect ∘ enumerate)(structures)
#         ss = Sr[structure = At(s)][Freq(gamma)]
#         ss[ss .< 0] .-= minimum(ss) # ! Ok?
#         df = ustrip(step(lookup(Sr[1], Freq)))
#         N = ss ./ (sum(ss, dims = 1) ./ df) # A density
#         fs = collect(lookup(ss, Freq))
#         μ = sum(fs .* N, dims = Freq) ./ df # Center of mass of gamma band
#         σ = sqrt.(sum((fs .- μ) .^ 2 .* N, dims = 1) ./ df)
#         σ = dropdims(σ, dims = Freq)
#         μ, (σl, σh) = bootstrapmedian(σ, dims = :sessionid)
#         μ, σl, σh = upsample.((μ, σl, σh), 5)
#         band!(ax, lookup(μ, 1), collect(σl), collect(σh);
#               color = (structurecolors[i], bandalpha), label = structures[i])
#         lines!(ax, lookup(μ, 1), μ;
#                color = (structurecolors[i], alpha), label = structures[i])
#     end

#     l = axislegend(ax, position = :lb, nbanks = 2, labelsize = 12, merge = true)
#     reverselegend!(l)
#     plotlayerints!(ax, layerints; axis = :x, newticks = false, flipside = true)
#     ax.limits = ((0, 1), (nothing, nothing))
#     display(f)
# end

addlabels!(f)
f |> display
save(plotdir("power_spectra", "power_spectra.pdf"), f)
save(plotdir("power_spectra", "power_spectra.svg"), f)