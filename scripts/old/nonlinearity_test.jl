#! /bin/bash
# -*- mode: julia -*-
#=
exec julia -t auto --startup-file=no --color=yes "${BASH_SOURCE[0]}" "$@"
=#

using DrWatson
@quickactivate "SpatiotemporalMotifs"

using SpatiotemporalMotifs
import TimeseriesTools: freqs
using Peaks
using FileIO
using SpatiotemporalMotifs
import SpatiotemporalMotifs.layers
@preamble
set_theme!(foresight(:physics))

sessionid = 1140102579
session = AN.Session(sessionid)
LFP = AN.formatlfp(session; stimulus = r"Natural_Images", structure = "VISp")
x = LFP[:, 6]

begin
    y = deepcopy(x)[1:10000]
    y = bandpass(y, 0.1 .. 100)
    # y = mod.(1:1:1000, 2π)
    # y = colorednoise(0.1:0.1:100)
    using TimeseriesSurrogates
    using ComplexityMeasures

    method = RandomFourier()
    function F(x)
        # est = MissingDispersionPatterns(Dispersion(m = 3, c = 7))
        est = SampleEntropy(x; m = 3, τ = 1)
        complexity_normalized(est, x)
    end
    # function F(x)
    #     mean((x[2:end] - x[1:(end - 1)]) .^ 3)
    # end
    # Buffer y and repeat over windows of length 10000
    T = SurrogateTest(F, parent(y), method; n = 1000)
    p = pvalue(T, tail = :left)
    ziggurat(T.vals)
    vlines!(F(parent(y)))
    current_figure()
end

begin # * Multiscale entropy
    using TimeseriesSurrogates
    using ComplexityMeasures
    y = deepcopy(x)
    # ts = 0.01:0.01:1000
    # y = TimeSeries(ts, mod.(ts, 2π) .^ 2 .+ randn(length(ts)) .* 2)
    s = surrogate(y, RandomFourier())
    using TimeseriesSurrogates
    using ComplexityMeasures
    # est = SampleEntropy(; m = 3, τ = 1, r = 0.2) # Assumes all subseries are standardized
    est = ApproximateEntropy(; m = 3, τ = 1, r = 0.2) # Assumes all subseries are standardized
    # est = Kraskov(Shannon(); k = 1, w = 0)
    # function F(x)
    #     complexity(est, x)
    # end
    function F(x)
        mean((x[2:end] .- x[1:(end - 1)]) .^ 3)
    end

    H = []
    N = 10000 # Number of samples in each entropy window
    n = 10000 # Maximum number of entropy windows
    λs = range(log10(1), log10((length(y) ÷ N)), length = 64) .|> exp10 .|> round .|>
         Int |> unique
    # λs = 1:10:(length(y) ÷ N)
    hs = progressmap(λs; parallel = true) do λ
        x = mean.(buffer(y, λ, 0))
        x = buffer(x, N, 0)
        x = x[randperm(length(x))[1:min(n, length(x))]]
        x = map(x) do y
            y = (y .- mean(y)) ./ std(y)
        end
        mean(F.(x))
    end
    hss = progressmap(λs; parallel = true) do λ
        s = surrogate(y, RandomFourier())
        x = mean.(buffer(s, λ, 0))
        x = buffer(x, N, 0)
        x = x[randperm(length(x))[1:min(n, length(x))]]
        x = map(x) do y
            y = (y .- mean(y)) ./ std(y)
        end
        mean(F.(x))
    end
    f = TwoPanel()
    ax = Axis(f[1, 1]; xscale = log10, xtickformat = x -> string.(round.(Int, x)))
    lines!(ax, 1 ./ (λs .* step(y)), hs, label = "Real")
    lines!(ax, 1 ./ (λs .* step(y)), hss, label = "Surrogate")
    axislegend(ax)
    ax = Axis(f[1, 2]; xscale = log10, xtickformat = x -> string.(round.(Int, x)))
    lines!(ax, 1 ./ (λs .* step(y)), hss .- hs)
    current_figure()
end
