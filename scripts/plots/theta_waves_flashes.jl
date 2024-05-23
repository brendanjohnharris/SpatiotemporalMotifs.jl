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
                             load_calculations, unify_calculations,
                             plotlayerints!,
                             @preamble
@preamble
set_theme!(foresight(:physics))

stimulus = "flash_250ms"

session_table = load(datadir("session_table.jld2"), "session_table")
oursessions = session_table.ecephys_session_id

path = datadir("calculations")
Q = calcquality(path)[structure = At(structures)]
quality = mean(Q[stimulus = At(stimulus)])

out = load_calculations(Q; path, stimulus);
uni = unify_calculations(out);

begin # * Calculate a global order parameter at each time point
    Og = map(out) do o
        O = map(o) do p
            k = p[:k]
            trials = p[:trials]
            dropdims(mean(sign.(k), dims = Dim{:depth}), dims = Dim{:depth})
        end
    end
end

begin # * Supplemental material: Average theta oscillation
    f = Figure(size = (720, 800))

    idxs = [[i, j] for i in 1:3, j in 1:2]
    idxs = permutedims(idxs, (2, 1))
    for i in eachindex(uni)
        x = uni[i][:x][:, 2:end, :, :] .* u"V"
        x = uconvert.(u"μV", x)
        ax = Axis(f[idxs[i]...], yreversed = true)
        structure = DimensionalData.metadata(x)[:structure]
        ax.title = structure
        m = median(x, dims = (:trial, :sessionid))
        m = dropdims(m, dims = (:trial, :sessionid))
        colorrange = maximum(abs.(ustripall(m))) * [-1, 1]

        p = heatmap!(ax, upsample(ustripall(m), 5, 2); colormap = :bone, colorrange)
        c = Colorbar(f[idxs[i]...][1, 2], p)

        vlines!(ax, [0, 0.25]; color = (:black, 0.5), linestyle = :dash, linewidth = 3)

        ints = uni[i][:layerints]
        plotlayerints!(ax, ints)
        if mod(i, 2) == 0
            c.label = "Mean theta LFP ($(unit(eltype(x))))"
        end
        if i > 4
            ax.xlabel = "Time (s)"
        end
    end
    display(f)
    wsave(plotdir("theta_waves_flashes", "supplemental_LFP.pdf"), f)
end

begin # * Supplemental material: single-trial phase velocity maps in each region
    f = Figure(size = (720, 800))

    idxs = [[i, j] for i in 1:3, j in 1:2]
    idxs = permutedims(idxs, (2, 1))
    for i in eachindex(uni)
        k = uni[i][:k][:, 2:end, :, :]
        k = uconvert.(u"mm^-1", k)
        ax = Axis(f[idxs[i]...], yreversed = true)
        structure = DimensionalData.metadata(k)[:structure]
        ax.title = structure
        m = median(k, dims = (:trial, :sessionid))
        m = dropdims(m, dims = (:trial, :sessionid))
        colorrange = maximum(abs.(ustripall(m))) * [-1, 1]

        p = heatmap!(ax, upsample(ustripall(m), 5, 2); colormap = binarysunset, colorrange)
        c = Colorbar(f[idxs[i]...][1, 2], p)

        q = ustripall((m))[1:100:end, 1:3:end]
        q ./= maximum(abs.(q))
        arrows!(ax, lookup(q, 1), lookup(q, 2), zeros(size(q)), parent(q);
                lengthscale = 0.07,
                normalize = false, color = (:black, 0.4))

        vlines!(ax, [0, 0.25]; color = (:white, 0.5), linestyle = :dash, linewidth = 3)

        ints = uni[i][:layerints]
        plotlayerints!(ax, ints)
        if mod(i, 2) == 0
            c.label = "Average wavenumber ($(unit(eltype(k))))"
        end
        if i > 4
            ax.xlabel = "Time (s)"
        end
    end
    display(f)
    wsave(plotdir("theta_waves_flashes", "supplemental_wavenumber.pdf"), f)
end

begin # * Plot the mean order parameter across time in each layer, and for a surrogate
    Ō = map(Og) do O
        o = dropdims.(mean.(O; dims = Dim{:changetime}), dims = Dim{:changetime})
        sessionids = [DimensionalData.metadata(_o)[:sessionid] for _o in O]
        stack(Dim{:sessionid}(sessionids), o)
    end

    f = Figure()
    ax = Axis(f[1, 1]; xlabel = "Time (s)", ylabel = "Order parameter",
              xautolimitmargin = (0, 0))
    hlines!(ax, [0]; color = (:black, 0.3), linestyle = :dash, linewidth = 1)
    for (i, O) in enumerate(Ō)
        structure = DimensionalData.metadata(O)[:structure]
        μ = dropdims(mean(O, dims = Dim{:sessionid}), dims = Dim{:sessionid})
        σ = dropdims(std(O, dims = Dim{:sessionid}), dims = Dim{:sessionid})
        σ = σ ./ sqrt(size(O, 2)) # SEM
        bargs = [times(μ), μ .- σ, μ .+ σ] .|> ustripall .|> collect
        band!(ax, bargs..., color = (structurecolors[i], 0.3))
        lines!(ax, μ |> ustripall, color = structurecolors[i], label = structure)
    end
    axislegend(ax, position = :lt)
    display(f)
    wsave(plotdir("theta_waves_flashes", "regional_orderparameter.pdf"), f)
end
