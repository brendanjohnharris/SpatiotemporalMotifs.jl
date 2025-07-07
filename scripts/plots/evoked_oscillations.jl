#! /bin/bash
#=
exec julia +1.10.10 -t auto "${BASH_SOURCE[0]}" "$@"
=#
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

begin # * Load the evoked oscillations for the flashes stimulus
    stimulus = "flash_250ms"
    uni = load_uni(; stimulus, vars = [:θ, :ω]) # Load theta component
    θ = getindex.(uni, :θ)
    ω = getindex.(uni, :ω)
end

begin # * Plot average evoked frequency
    f = Figure()
    ax = Axis(f[1, 1], xlabel = "Frequency (Hz)", ylabel = "Cortical depth",
              title = "Average evoked frequency for flashes stimulus")
    for i in 1:6
        x = ω[i] ./ (2π)
        x = x[𝑡 = 0.25u"s" .. 0.75u"s"]
        x = mean(x, dims = (1, 3, 4))
        x = dropdims(x, dims = (1, 3, 4))
        lines!(ax, reverse(decompose(ustripall(x)))..., color = structurecolors[i],
               label = structures[i])
    end
    axislegend(ax, position = :lt, nbanks = 2)
    display(f)
end
