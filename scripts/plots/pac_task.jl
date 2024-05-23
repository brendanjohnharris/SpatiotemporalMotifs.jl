#! /bin/bash
#=
exec julia -t auto "${BASH_SOURCE[0]}" "$@"
=#
using DrWatson
@quickactivate "SpatiotemporalMotifs"

import TimeseriesTools: freqs
using FileIO
using Unitful
using Normalization
import SpatiotemporalMotifs: plotdir, calcquality, layernum2name, savepath,
                             structurecolors, structures, commondepths, parselayernum,
                             load_calculations, unify_calculations, plotlayerints!,
                             plotlayermap!,
                             @preamble
@preamble
set_theme!(foresight(:physics))

stimulus = r"Natural_Images"
vars = [:ϕ, :aᵧ]
nmin = 25

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
unidepths = getindex.(uni, :unidepths)
layerints = getindex.(uni, :layerints)
layernames = getindex.(uni, :layernames)
layernums = getindex.(uni, :layernums)

begin # * Normalize amplitudes and generate a burst mask
    r = [abs.(uni[i][:aᵧ]) for i in eachindex(uni)]

    r_h = [r[i][:, :, lookup(r[i], :trial) .== true] for i in eachindex(r)]
    r̂_h = [HalfZScore(_r, dims = [1, 3])(_r) for _r in r_h] # * Ok to normalize over subjects?

    r_m = [r[i][:, :, lookup(r[i], :trial) .== false] for i in eachindex(r)]
    r̂_m = [HalfZScore(_r, dims = [1, 3])(_r) for _r in r_m] # * Ok to normalize over subjects?

    r = []
    GC.gc()
end
begin
    mask_h = [r̂_h[i] .> 2.0 for i in eachindex(r̂_h)]
    mask_m = [r̂_m[i] .> 2.0 for i in eachindex(r̂_m)]

    ϕ = [abs.(uni[i][:ϕ]) for i in eachindex(uni)]

    ϕ_h = [ϕ[i][:, :, lookup(ϕ[i], :trial) .== true] for i in eachindex(ϕ)]
    ϕ̂_h = deepcopy(ϕ_h)
    [ϕ̂_h[i][.!mask_h[i]] .= NaN for i in eachindex(mask_h)]

    ϕ_m = [ϕ[i][:, :, lookup(ϕ[i], :trial) .== false] for i in eachindex(ϕ)]
    ϕ̂_m = deepcopy(ϕ_m)
    [ϕ̂_m[i][.!mask_m[i]] .= NaN for i in eachindex(mask_m)]

    ϕ = []
    uni = []
    GC.gc()
end

begin # * Heatmap of the number of samples
    pn_h = dropdims.(sum.(mask_h, dims = 3), dims = 3)
    pn_m = dropdims.(sum.(mask_m, dims = 3), dims = 3)

    f = Figure(size = (720, 1440))

    for i in eachindex(pn_h)
        # * Hit
        ax = Axis(f[i, 1], yreversed = true)
        structure = DimensionalData.metadata(pn_h[i])[:structure]
        ax.title = structure * ": hit"
        m = ustripall(pn_h[i])
        colorrange = (nmin, max(maximum(m), maximum(ustripall(pn_m[i]))))
        ints = layerints[i]
        p = plotlayermap!(ax, m, ints; arrows = false, colorrange,
                          colormap = :bone, lowclip = :crimson) |>
            first
        if i > 5
            ax.xlabel = "Time (s)"
        end

        # * Miss
        ax = Axis(f[i, 2], yreversed = true)
        ax.title = structure * ": miss"
        m = ustripall(pn_m[i])
        p = plotlayermap!(ax, m, ints; arrows = false, colorrange,
                          colormap = :bone, lowclip = :crimson) |>
            first
        c = Colorbar(f[i, 3], p)
        c.label = "Num. of bursts"
    end

    wsave(plotdir("pac_task", "supplemental_numsamples.pdf"), f)
end

begin # * Heatmap of the PAC entropy
    function pacentropy(ϕ) # ?! Normalised or unnormalized? Domain is always the same, so...?
        _p = dropdims(nansafe(x -> StateSpaceSet(x); dims = 3)(ϕ), dims = 3)
        function ent(x)
            if isempty(x)
                return 3 # Some entropy higher than the maximum possible, 1
            else
                ComplexityMeasures.information_normalized(Shannon(),
                                                          ValueBinning(RectangularBinning(20)),
                                                          StateSpaceSet(x))
            end
        end
        return ent.(_p)
    end

    PAC_h = pacentropy.(ϕ̂_h)
    PAC_m = pacentropy.(ϕ̂_m)

    f = Figure(size = (720, 1440))

    for i in eachindex(PAC_h)
        # * Hit
        ax = Axis(f[i, 1], yreversed = true)
        structure = DimensionalData.metadata(pn_h[i])[:structure]
        ax.title = structure * ": hit"
        m = ustripall(PAC_h[i])
        colorrange = [extrema(m)..., extrema(ustripall(PAC_m[i]))...] |> extrema
        m[pn_h[i] .< nmin] .= 2
        ints = layerints[i]
        p = plotlayermap!(ax, m, ints; arrows = false, colormap = :bone, colorrange,
                          highclip = :crimson) |>
            first
        if i > 5
            ax.xlabel = "Time (s)"
        end

        # * Miss
        ax = Axis(f[i, 2], yreversed = true)
        ax.title = structure * ": miss"
        m = ustripall(PAC_m[i])
        m[pn_m[i] .< nmin] .= 2
        p = plotlayermap!(ax, m, ints; arrows = false, colorrange, highclip = :crimson,
                          colormap = :bone) |>
            first
        c = Colorbar(f[i, 3], p)
        c.label = "PAC entropy"
    end
    wsave(plotdir("pac_task", "supplemental_PAC.pdf"), f)
end
