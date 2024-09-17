#! /bin/bash
#=
exec julia -t auto "${BASH_SOURCE[0]}" "$@"
=#
# ? Expected execution time: 30 mins
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
Random.seed!(32)

vars = [:csd, :k, :ω]
INTERVAL = SpatiotemporalMotifs.INTERVAL
mainstructure = "VISl"
maincolorrange = [-2.0, 2.0]
session_table = load(datadir("posthoc_session_table.jld2"), "session_table")
oursessions = session_table.ecephys_session_id
path = datadir("calculations")
statsfile = plotdir("theta_wavenumbers", "$mainstructure.txt")
close(open(statsfile, "w"))
ylims = (0.05, 0.95)

begin # * Set up main figure
    mf = TwoPanel()
    mfs = subdivide(mf, 1, 4)
end

stimulus = r"Natural_Images"
uni = load_uni(; stimulus, vars)

# begin # * Current source density
#     csd = getindex.(uni, :csd)[2][:, :, :]
#     csd = dropdims(mean(csd, dims = Trial); dims = Trial)
#     csd = upsample(ustripall(csd), 5, Depth) # For VISl
#     # csd = set(csd, csd[:, end:-1:1])
#     f = Figure()
#     ax = Axis(f[1, 1]; yreversed = true)
#     plotlayermap!(ax, csd)
#     plotlayerints!(ax, uni[2][:layerints])
#     # heatmap(decompose(csd)...; colormap = binarysunset, axis = (; yreversed = true))
#     f
# end

begin # * Supplemental material: average wavenumbers in each region
    f = Figure(size = (720, 1440))

    for i in eachindex(uni)
        k = uni[i][:k]
        ω = uni[i][:ω]
        k = uconvert.(u"mm^-1", k)
        k[ustripall(ω) .< 0] .= NaN * unit(eltype(k)) # Mask out negative frequency periods
        @info "$(round(100*sum(ω .< 0u"Hz")/length(ω), sigdigits=2))% of data has a negative frequency"

        # * Hit
        ax = Axis(f[i, 1], yreversed = true)
        ax.limits = (nothing, ylims)
        structure = metadata(k)[:structure]
        ax.title = structure * ": hit"
        m = nansafe(median, dims = Trial)(k[:, :, lookup(k, Trial) .== true])
        m = dropdims(m, dims = Trial)[𝑡(SpatiotemporalMotifs.INTERVAL)]
        colorrange = maximum(abs.(ustripall(m))) * [-1, 1]
        ints = uni[i][:layerints]
        p = plotlayermap!(ax, m, ints; arrows = true, colorrange) |> first
        if i > 5
            ax.xlabel = "Time (s)"
        end

        if structure == mainstructure # * Plot into main figure
            ax = Axis(mfs[1], yreversed = true)
            ax.xlabel = "Time (s)"
            ax.limits = (nothing, ylims)
            ax.title = structure * ": hit"
            m = nansafe(median, dims = Trial)(k[:, :, lookup(k, Trial) .== true])
            m = dropdims(m, dims = Trial)[𝑡(SpatiotemporalMotifs.INTERVAL)]
            ints = uni[i][:layerints]
            p = plotlayermap!(ax, m, ints; arrows = true, colorrange) |> first
            c = Colorbar(mfs[end]; colorrange = maincolorrange, colormap = defaultcolormap,
                         highclip = defaultcolormap[end], lowclip = defaultcolormap[1])
            c.label = "θ wavenumber ($(unit(eltype(k))))"

            # * Calculate some stats of this figure
            open(statsfile, "a+") do file
                write(file, "## $mainstructure average wavenumbers\n")
                write(file,
                      "Average wavenumber (median ± IQR) = $(only(nansafe(median)(k[:]))) ± $(nansafe(iqr)(k[:])|>only)\n")
                write(file,
                      "Average wavenumber magnitude (median ± IQR) = $(only(nansafe(median)(abs.(k[:])))) ± $(nansafe(iqr)(abs.(k[:]))|>only)\n")
            end
        end

        # * Miss
        ax = Axis(f[i, 2], yreversed = true)
        ax.limits = (nothing, ylims)
        ax.title = structure * ": miss"
        m = nansafe(median, dims = Trial)(k[:, :, lookup(k, Trial) .== false])
        m = dropdims(m, dims = Trial)[𝑡(SpatiotemporalMotifs.INTERVAL)]
        ints = uni[i][:layerints]
        p = plotlayermap!(ax, m, ints; arrows = true, colorrange) |> first
        c = Colorbar(f[i, 3], p)
        c.label = "θ wavenumber ($(unit(eltype(k))))"

        if structure == mainstructure # * Plot into main figure
            ax = Axis(mfs[2], yreversed = true)
            ax.xlabel = "Time (s)"
            ax.limits = (nothing, ylims)
            ax.title = structure * ": miss"
            m = nansafe(median; dims = Trial)(k[:, :, lookup(k, Trial) .== false])
            m = dropdims(m, dims = Trial)[𝑡(SpatiotemporalMotifs.INTERVAL)]
            ints = uni[i][:layerints]
            p = plotlayermap!(ax, m, ints; arrows = true, colorrange = maincolorrange) |>
                first
        end
    end
    addlabels!(f)
    display(f)
    wsave(plotdir("theta_wavenumbers", "supplemental_wavenumber.pdf"), f)
end

begin # * Supplemental material: csd in each region
    f = Figure(size = (720, 1440))

    for i in eachindex(uni)
        k = uni[i][:csd]

        # * Hit
        ax = Axis(f[i, 1], yreversed = true)
        ax.limits = (nothing, ylims)
        structure = metadata(k)[:structure]
        ax.title = structure * ": hit"
        m = nansafe(median, dims = Trial)(k[:, :, lookup(k, Trial) .== true])
        m = dropdims(m, dims = Trial)[𝑡(SpatiotemporalMotifs.INTERVAL)]
        colorrange = maximum(abs.(ustripall(m))) * [-1, 1]
        ints = uni[i][:layerints]
        p = plotlayermap!(ax, m, ints; arrows = false, colorrange,
                          colormap = binarysunset) |> first
        if i > 5
            ax.xlabel = "Time (s)"
        end

        # if structure == mainstructure # * Plot into main figure
        #     ax = Axis(mfs[1], yreversed = true)
        #     ax.xlabel = "Time (s)"
        #     ax.limits = (nothing, ylims)
        #     ax.title = structure * ": hit"
        #     m = nansafe(median, dims = Trial)(k[:, :, lookup(k, Trial) .== true])
        #     m = dropdims(m, dims = Trial)[𝑡(SpatiotemporalMotifs.INTERVAL)]
        #     ints = uni[i][:layerints]
        #     p = plotlayermap!(ax, m, ints; arrows = true, colorrange) |> first
        #     c = Colorbar(mfs[end]; colorrange = maincolorrange, colormap = defaultcolormap,
        #                  highclip = defaultcolormap[end], lowclip = defaultcolormap[1])
        #     c.label = "CSD ($(unit(eltype(k))))"
        # end

        # * Miss
        ax = Axis(f[i, 2], yreversed = true)
        ax.limits = (nothing, ylims)
        ax.title = structure * ": miss"
        m = nansafe(median, dims = Trial)(k[:, :, lookup(k, Trial) .== false])
        m = dropdims(m, dims = Trial)[𝑡(SpatiotemporalMotifs.INTERVAL)]
        ints = uni[i][:layerints]
        p = plotlayermap!(ax, m, ints; arrows = false, colorrange,
                          colormap = binarysunset) |> first
        c = Colorbar(f[i, 3], p)
        c.label = "CSD ($(u"V"*unit(eltype(k))))" # Correct for dropped units by multiplying with V

        # if structure == mainstructure # * Plot into main figure
        #     ax = Axis(mfs[2], yreversed = true)
        #     ax.xlabel = "Time (s)"
        #     ax.limits = (nothing, ylims)
        #     ax.title = structure * ": miss"
        #     m = nansafe(median; dims = Trial)(k[:, :, lookup(k, Trial) .== false])
        #     m = dropdims(m, dims = Trial)[𝑡(SpatiotemporalMotifs.INTERVAL)]
        #     ints = uni[i][:layerints]
        #     p = plotlayermap!(ax, m, ints; arrows = true, colorrange = maincolorrange) |>
        #         first
        # end
    end
    addlabels!(f)
    display(f)
    wsave(plotdir("theta_wavenumbers", "supplemental_csd.pdf"), f)
end

# ? Flashes
stimulus = "flash_250ms"
uni = load_uni(; stimulus, vars)

begin # * Supplemental material: average wavenumbers in each region
    f = SixPanel()
    gs = subdivide(f, 3, 2)

    for i in eachindex(uni)
        k = uni[i][:k]
        ω = uni[i][:ω]
        k = uconvert.(u"mm^-1", k)
        k[ustripall(ω) .< 0] .= NaN * unit(eltype(k)) # Mask out negative frequency periods
        structure = metadata(k)[:structure]

        ax = Axis(gs[i][1, 1], yreversed = true, xlabel = "Time (s)")
        ax.limits = (nothing, ylims)
        ax.title = structure
        m = nansafe(median, dims = (Trial, SessionID))(k)
        m = dropdims(m, dims = (Trial, SessionID))[𝑡(SpatiotemporalMotifs.INTERVAL)]
        ints = uni[i][:layerints]
        p = plotlayermap!(ax, m, ints; arrows = true) |> first
        c = Colorbar(gs[i][1, 2], p)
        c.label = "θ wavenumber ($(unit(eltype(k))))"
        if structure == mainstructure
            colorrange = maincolorrange
            ax = Axis(mfs[3], yreversed = true, xlabel = "Time (s)")
            ax.limits = (nothing, ylims)
            ax.title = structure * ": flashes"
            m = nansafe(median, dims = (Trial, SessionID))(k)
            m = dropdims(m, dims = (Trial, SessionID))[𝑡(SpatiotemporalMotifs.INTERVAL)]
            ints = uni[i][:layerints]
            p = plotlayermap!(ax, m, ints; arrows = true, colorrange) |> first
        end
    end
    addlabels!(f)
    display(f)
    wsave(plotdir("theta_wavenumbers", "supplemental_wavenumber_flash.pdf"), f)
end

begin # * Save main figure
    addlabels!(mf)
    wsave(plotdir("theta_wavenumbers", "theta_wavenumber.pdf"), mf)
    display(mf)
end
