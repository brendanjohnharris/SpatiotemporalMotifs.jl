#! /bin/bash
# -*- mode: julia -*-
#=
exec julia +1.10.9 -t auto --color=yes "${BASH_SOURCE[0]}" "$@"
=#
using DrWatson
@quickactivate "SpatiotemporalMotifs"

using SpatiotemporalMotifs
import TimeseriesTools: freqs
using Peaks
using FileIO
using SpatiotemporalMotifs
import SpatiotemporalMotifs.layers
using Random
@preamble
set_theme!(foresight(:physics))

begin # * Load the CSD  flashes
    path = datadir("calculations")
    Q = calcquality(path)[Structure = At(structures)]
    outf = load_calculations(Q; stimulus = "flash_250ms", vars = [:csd, :V])
end
begin # * Load the CSD passive natural images
    path = datadir("calculations")
    Q = calcquality(path)[Structure = At(structures)]
    out = load_calculations(Q; stimulus = "Natural_Images_passive", vars = [:csd, :V])
end

begin
    s = 2
    i = 19
    begin # * Our csd
        csdf = outf[s][i][:csd][𝑡 = 1:1562]
        csdf = mapslices(median, csdf, dims = :changetime)
        csdf = dropdims(csdf, dims = :changetime)
        heatmap(csdf |> ustripall,
                axis = (; yreversed = true, limits = ((0, 0.1), nothing)))
    end
    begin # * Our csd
        csd = out[s][i][:csd][𝑡 = 1:1562]
        csd = mapslices(median, csd, dims = :changetime)
        csd = dropdims(csd, dims = :changetime)
        heatmap(csd |> ustripall,
                axis = (; yreversed = true, limits = ((0, 0.1), nothing)))
    end

    begin
        f = Figure()
        ax = Axis(f[1, 1], yreversed = true)
        x = mean([csd, csdf])
        heatmap!(ax, ustripall(x))
        ax.limits = ((0, 0.1), nothing)
        display(f)
    end
end
function find_L4_center(csd::AbstractMatrix,
                        depths;
                        window::Tuple{Real, Real} = (0, 50),
                        fs::Real = 1_000)
    n_time, n_depth = size(csd)
    @assert length(depths)==n_depth "depths length ≠ csd columns"

    # ----- 1. five-point Hamming smoothing (Ulbert et al., 2001) -----
    csd_sm = csd

    # ----- 2. restrict to analysis window -----
    s1 = max(1, round(Int, window[1] / 1_000 * fs) + 1)
    s2 = min(n_time, round(Int, window[2] / 1_000 * fs) + 1)
    @assert s1<s2 "analysis window lies outside data range"

    # ----- 3-5. find earliest, largest sink -----
    earliest = n_time                # sentinel
    best_idx = 1
    best_amp = Inf                   # most-negative ⇒ smaller numbers

    for d in 1:n_depth
        col = view(csd_sm, s1:s2, d)
        amp, rel_idx = findmin(col)          # amp is negative
        t_idx = s1 + rel_idx - 1

        if t_idx < earliest ||
           (t_idx == earliest && amp < best_amp)
            earliest = t_idx
            best_idx = d
            best_amp = amp
        end
    end

    return depths[best_idx]           # L4 centre (µm from pia)
end

begin
    l4 = find_L4_center(parent(csd), lookup(csd, :Depth); window = (0, 50), fs = 1250)
    f = Figure()
    ax = Axis(f[1, 1], title = "L4 center depth", yreversed = true)
    heatmap!(ax, ustripall(csd))
    hlines!(ax, [l4], color = :red, label = "L4 center")
    ax.limits = ((0, 0.25), nothing)
    display(f)
end
