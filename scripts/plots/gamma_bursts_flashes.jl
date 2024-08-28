#! /bin/bash
#=
exec julia -t auto "${BASH_SOURCE[0]}" "$@"
=#
using DrWatson
@quickactivate "SpatiotemporalMotifs"

import TimeseriesTools: freqs
using FileIO
using Unitful
import SpatiotemporalMotifs: plotdir, calcquality, layernum2name, savepath,
                             structurecolors, structures, commondepths, parselayernum,
                             load_calculations, unify_calculations, plotlayerints!,
                             @preamble
@preamble
set_theme!(foresight(:physics))

stimulus = "flash_250ms"
vars = [:y, :aᵧ]

session_table = load(datadir("session_table.jld2"), "session_table")
oursessions = session_table.ecephys_session_id

path = datadir("calculations")
Q = calcquality(path)[structure = At(structures)]
quality = mean(Q[stimulus = At(stimulus)])

config = @strdict stimulus vars
data, file = produce_or_load(config, datadir(); filename = savepath) do config
    out = load_calculations(Q; path, stimulus, vars)
    uni = unify_calculations(out; vars)
    @strdict uni
end
uni = data["uni"]

begin # * Supplemental material: Average gamma oscillation
    f = Figure(size = (720, 800))

    idxs = [[i, j] for i in 1:3, j in 1:2]
    idxs = permutedims(idxs, (2, 1))
    for i in eachindex(uni)
        ax = Axis(f[idxs[i]...], yreversed = true)
        m = median(uni[i][:y], dims = (Trial, SessionID))
        structure = DimensionalData.metadata(m)[:structure]
        ax.title = structure
        m = m .* u"V"
        m = uconvert.(u"μV", m)
        m = dropdims(m, dims = (Trial, SessionID))
        colorrange = maximum(abs.(ustripall(m))) * [-1, 1]

        p = heatmap!(ax, upsample(ustripall(m), 5, 2); colormap = :bone, colorrange)
        c = Colorbar(f[idxs[i]...][1, 2], p)

        vlines!(ax, [0, 0.25]; color = (:black, 0.5), linestyle = :dash, linewidth = 3)

        ints = uni[i][:layerints]
        plotlayerints!(ax, ints)
        if mod(i, 2) == 0
            c.label = "Mean gamma LFP ($(unit(eltype(m))))"
        end
        if i > 4
            ax.xlabel = "Time (s)"
        end
    end
    display(f)
    wsave(plotdir("gamma_bursts_flashes", "supplemental_LFP.pdf"), f)
end

begin # * Supplemental material: Mean gamma amplitude maps in each region
    f = Figure(size = (720, 800))

    idxs = [[i, j] for i in 1:3, j in 1:2]
    idxs = permutedims(idxs, (2, 1))
    for i in eachindex(uni)
        m = median(abs.(uni[i][:aᵧ]), dims = (Trial, SessionID))
        m = m .* u"V"
        m = uconvert.(u"μV", m)
        m = dropdims(m, dims = (Trial, SessionID))

        ax = Axis(f[idxs[i]...], yreversed = true)
        structure = DimensionalData.metadata(m)[:structure]
        ax.title = structure
        colorrange = maximum(ustripall(m)) * [0, 1]

        p = heatmap!(ax, upsample(ustripall(m), 5, 2); colormap = :bone, colorrange)
        c = Colorbar(f[idxs[i]...][1, 2], p)

        vlines!(ax, [0, 0.25]; color = (:white, 0.5), linestyle = :dash, linewidth = 3)

        ints = uni[i][:layerints]
        plotlayerints!(ax, ints)
        if mod(i, 2) == 0
            c.label = "Average 𝑟 ($(unit(eltype(m))))"
        end
        if i > 4
            ax.xlabel = "Time (s)"
        end
    end
    display(f)
    wsave(plotdir("gamma_bursts_flashes", "supplemental_amplitude.pdf"), f)
end
