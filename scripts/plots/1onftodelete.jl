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
using Random
@preamble
set_theme!(foresight(:physics))

plot_data, data_file = produce_or_load(Dict(), calcdir("plots");
                                       filename = savepath,
                                       prefix = "figS2") do _
end

begin # * Load the CSD for the example session
    file = "/taiji1/bhar9988/code/DDC/SpatiotemporalMotifs/data/plots/fig2.jld2"
    f = jldopen(file, "r")
end
begin
    is = (10, 50)
    begin # * Without loglinear sampling
        S = f["spontaneous"]["S"][1][:, is...]
        ff, D = fooof(S |> ustripall)
        χ = D[:χ]
        fig = Figure()
        ax = Axis(fig[1, 1], xscale = log10, yscale = log10)
        lines!(ax, S |> ustripall)

        # * Plot the fooof fit
        xs = lookup(S, 1) |> ustripall
        ys = ff.(xs)
        lines!(ax, xs, ys, color = crimson, label = "fooof fit")
        println("uncorrected exponent: ", χ)
    end
    begin # * With loglinear sampling
        _S = f["spontaneous"]["S"][1][:, is...] |> ustripall

        _idxs = log10.(freqs(_S))
        idxs = Interval(extrema(_idxs)...)
        idxs = range(idxs, step = maximum(diff(_idxs))) .|> exp10

        S = _S[𝑓 = Near(idxs)]

        ff, D = fooof(S |> ustripall)
        χ = D[:χ]
        scatter!(ax, S)

        # * Plot the fooof fit
        xs = lookup(S, 1)
        ys = ff.(xs)
        lines!(ax, xs, ys, color = cucumber, label = "fooof fit")
        fig |> display
        println("corrected exponent: ", χ)
    end
end

begin # * Do the same on a langevin simulation
end
