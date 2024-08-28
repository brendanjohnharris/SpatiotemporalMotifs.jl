#! /bin/bash
#=
exec julia -t auto "${BASH_SOURCE[0]}" "$@"
=#
# ? Expected execution time: 15 mins
using DrWatson
@quickactivate "SpatiotemporalMotifs"

import TimeseriesTools: freqs
using FileIO
using Unitful
import DimensionalData: metadata
using MultivariateStats
using SpatiotemporalMotifs
structurecolors, structures, commondepths, parselayernum,
load_calculations, unify_calculations,
plotlayerints!, plotlayermap!,
@preamble
set_theme!(foresight(:physics))
Random.seed!(32)

stimulus = r"Natural_Images"
datafile = datadir("theta_waves_task.jld2")
INT = 𝑡(SpatiotemporalMotifs.INTERVAL)

session_table = load(datadir("session_table.jld2"), "session_table")
oursessions = session_table.ecephys_session_id

path = datadir("calculations")
Q = calcquality(path)[structure = At(structures)]
quality = mean(Q[stimulus = At(stimulus)])

# * Negative frequencies
vars = [:ω]
config = @strdict stimulus vars
data, file = produce_or_load(produce_uni, config, datadir(); filename = savepath)
uni = data["uni"]

begin # * Supplemental material: likelihood of negative frequencies
    F1 = FourPanel()
    F2 = FourPanel()
    g1 = subdivide(F1, 3, 1)
    g2 = subdivide(F2, 3, 1)

    for i in eachindex(uni)
        f = i > 3 ? g2[i - 3] : g1[i]
        k = ustripall(uni[i][:ω][INT]) .< 0
        # k = uconvert.(u"mm^-1", k)

        # * Hit
        ax = Axis(f[1, 1], yreversed = true, limits = (nothing, (0.05, 0.95)))
        structure = metadata(k)[:structure]
        ax.title = structure * ": hit"
        m = mean(k[:, :, lookup(k, Trial) .== true], dims = Trial)
        m = dropdims(m, dims = Trial)
        colorrange = maximum(abs.(ustripall(m))) * [0, 1]
        ints = uni[i][:layerints]
        p = plotlayermap!(ax, m, ints; arrows = (200, 3), colorrange) |> first
        if i > 5
            ax.xlabel = "Time (s)"
        end

        # * Miss
        ax = Axis(f[1, 2], yreversed = true)
        ax.title = structure * ": miss"
        m = mean(k[:, :, lookup(k, Trial) .== false], dims = Trial)
        m = dropdims(m, dims = Trial)
        ints = uni[i][:layerints]
        p = plotlayermap!(ax, m, ints; colorrange) |> first
        c = Colorbar(f[1, 3], p)
        c.label = "P(ω < 0)"
    end
    display(F1)
    wsave(plotdir("theta_waves_task", "supplemental_negative_frequencies_a.pdf"), F1)
    wsave(plotdir("theta_waves_task", "supplemental_negative_frequencies_b.pdf"), F2)
end

# * Theta wavenumber
vars = [:k]
config = @strdict stimulus vars
data, file = produce_or_load(produce_uni, config, datadir(); filename = savepath)
uni = data["uni"]

begin # * Supplemental material: average phase velocity maps in each region
    F1 = FourPanel()
    F2 = FourPanel()
    g1 = subdivide(F1, 3, 1)
    g2 = subdivide(F2, 3, 1)

    for i in eachindex(uni)
        f = i > 3 ? g2[i - 3] : g1[i]
        k = uni[i][:k][:, 2:end, :][INT]
        k = uconvert.(u"mm^-1", k)

        # * Hit
        ax = Axis(f[1, 1], yreversed = true)
        structure = metadata(k)[:structure]
        ax.title = structure * ": hit"
        m = median(k[:, :, lookup(k, Trial) .== true], dims = Trial)
        m = dropdims(m, dims = Trial)
        colorrange = maximum(abs.(ustripall(m))) * [-1, 1]
        ints = uni[i][:layerints]
        p = plotlayermap!(ax, m, ints; arrows = (200, 3), colorrange) |> first
        if i > 5
            ax.xlabel = "Time (s)"
        end

        # * Miss
        ax = Axis(f[1, 2], yreversed = true)
        ax.title = structure * ": miss"
        m = median(k[:, :, lookup(k, Trial) .== false], dims = Trial)
        m = dropdims(m, dims = Trial)
        ints = uni[i][:layerints]
        p = plotlayermap!(ax, m, ints; arrows = (200, 3), colorrange) |> first
        c = Colorbar(f[1, 3], p)
        c.label = "k̃ ($(unit(eltype(k))))"
    end
    display(F1)
    wsave(plotdir("theta_waves_task", "supplemental_wavenumber_a.pdf"), F1)
    wsave(plotdir("theta_waves_task", "supplemental_wavenumber_b.pdf"), F2)
end

# * Theta LFP
vars = [:x]
config = @strdict stimulus vars
data, file = produce_or_load(produce_uni, config, datadir(); filename = savepath)
uni = data["uni"]

begin # * Supplemental material: average phase velocity maps in each region
    F1 = FourPanel()
    F2 = FourPanel()
    g1 = subdivide(F1, 3, 1)
    g2 = subdivide(F2, 3, 1)

    for i in eachindex(uni)
        f = i > 3 ? g2[i - 3] : g1[i]
        k = uni[i][:x][:, 2:end, :][INT]
        k = k .* 1000 # V to mv

        # * Hit
        ax = Axis(f[1, 1], yreversed = true)
        structure = metadata(k)[:structure]
        ax.title = structure * ": hit"
        m = median(k[:, :, lookup(k, Trial) .== true], dims = Trial)
        m = dropdims(m, dims = Trial)
        colorrange = maximum(abs.(ustripall(m))) * [-1, 1]
        ints = uni[i][:layerints]
        p = plotlayermap!(ax, m, ints; colorrange) |> first
        if i > 5
            ax.xlabel = "Time (s)"
        end

        # * Miss
        ax = Axis(f[1, 2], yreversed = true)
        ax.title = structure * ": miss"
        m = median(k[:, :, lookup(k, Trial) .== false], dims = Trial)
        m = dropdims(m, dims = Trial)
        ints = uni[i][:layerints]
        p = plotlayermap!(ax, m, ints; colorrange) |> first
        c = Colorbar(f[1, 3], p)
        c.label = "θ̃ (mV)"
    end
    display(F1)
    wsave(plotdir("supplemental_heatmaps", "supplemental_theta_a.pdf"), F1)
    wsave(plotdir("supplemental_heatmaps", "supplemental_theta_b.pdf"), F2)
end

# * Gamma amplitude
vars = [:r]
config = @strdict stimulus vars
data, file = produce_or_load(produce_uni, config, datadir(); filename = savepath)
uni = data["uni"]

begin # * Supplemental material: average phase velocity maps in each region
    F1 = FourPanel()
    F2 = FourPanel()
    g1 = subdivide(F1, 3, 1)
    g2 = subdivide(F2, 3, 1)

    for i in eachindex(uni)
        f = i > 3 ? g2[i - 3] : g1[i]
        k = uni[i][:r][:, 3:end, :][INT]
        k = k .* 1000 # V to mv

        # * Hit
        ax = Axis(f[1, 1], yreversed = true)
        structure = metadata(k)[:structure]
        ax.title = structure * ": hit"
        m = median(k[:, :, lookup(k, Trial) .== true], dims = Trial)
        m = dropdims(m, dims = Trial)
        colorrange = extrema(ustripall(m))
        ints = uni[i][:layerints]
        p = plotlayermap!(ax, m, ints; colorrange) |> first
        if i > 5
            ax.xlabel = "Time (s)"
        end

        # * Miss
        ax = Axis(f[1, 2], yreversed = true)
        ax.title = structure * ": miss"
        m = median(k[:, :, lookup(k, Trial) .== false], dims = Trial)
        m = dropdims(m, dims = Trial)
        ints = uni[i][:layerints]
        p = plotlayermap!(ax, m, ints; colorrange) |> first
        c = Colorbar(f[1, 3], p)
        c.label = "r̃ (mV)"
    end
    display(F1)
    wsave(plotdir("supplemental_heatmaps", "supplemental_gamma_a.pdf"), F1)
    wsave(plotdir("supplemental_heatmaps", "supplemental_gamma_b.pdf"), F2)
end
