#! /bin/bash
# -*- mode: julia -*-
#=
exec julia +1.10.10 -t auto --color=yes "${BASH_SOURCE[0]}" "$@"
=#
using DrWatson
@quickactivate "SpatiotemporalMotifs"

using SpatiotemporalMotifs
import TimeseriesTools: freqs
using SpatiotemporalMotifs
using Distributed

using Random
@preamble
set_theme!(foresight(:physics))

# begin # * Load the trial LFP for natural images
#     session_table = load(datadir("posthoc_session_table.jld2"), "session_table")
#     oursessions = session_table.ecephys_session_id
#     path = datadir("calculations")
#     Q = calcquality(path)[Structure = At(structures)]
#     Q = Q[SessionID(At(oursessions))]
#     @assert mean(Q[stimulus = At(r"Natural_Images")]) == 1
#     out = load_calculations(Q; stimulus = r"Natural_Images", vars = [:ϕ, :x])
# end

begin # * extract hit trials and plot theta phase
    o = out[2]
    x = map(o) do x
        y = x[:ϕ][:, :, .!x[:trials].hit]
        y = circularmean(y, dims = :changetime)
        dropdims(y, dims = 3)
    end
end

begin
    a = x[45] |> ustripall
    f = Figure()
    ax = Axis3(f[1, 1])

    xs = lookup(a, 𝑡) |> parent
    ys = lookup(a, Depth)
    zs = eachcol(parent(a))

    colorrange = extrema(ys)

    map(zs, ys) do z, y
        lines!(ax, xs, fill(y, length(xs)), z; color = y, colorrange,
               colormap = binarysunset)
    end
    display(f)
end
begin
    a = x[45] |> ustripall
    a = a[𝑡 = -0.05 .. 0.3]
    f = Figure()
    ax = Axis3(f[1, 1], xreversed = true)
    surface!(ax, ustripall(a))
    ax.azimuth = 1.4
    ax.elevation = 0.6

    display(f)
end

begin #  * Plot peak and trough contours
    a = x[45] |> ustripall
    a = a[𝑡 = -0.05 .. 0.3]
    f = Figure()
    ax = Axis(f[1, 1], yreversed = true)

    heatmap!(ax, ustripall(a), colormap = cyclic)
    display(f)
end
