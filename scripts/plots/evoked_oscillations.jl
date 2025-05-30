#! /bin/bash
#=
exec julia +1.10.9 -t auto "${BASH_SOURCE[0]}" "$@"
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
    uni = load_uni(; stimulus, vars = [:x, :ω]) # Load theta component
    θ = getindex.(uni, :x)
    ω = getindex.(uni, :ω)
end

begin # * Plot average evoked frequency
    f = Figure()
    ax = Axis(f[1, 1])
    for i in 1:6
        x = ω[i] ./ (2π)
        x = x[𝑡 = 0.25u"s" .. 0.75u"s"]
        x = mean(x, dims = (1, 3, 4))
        x = dropdims(x, dims = (1, 3, 4))
        lines!(ax, decompose(ustripall(x))...)
    end
    display(f)
end
