#! /bin/bash
# -*- mode: julia -*-
#=
exec julia -t auto --color=yes "${BASH_SOURCE[0]}" "$@"
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

sessionid = 1052533639
session = AN.Session(sessionid)
LFP = AN.formatlfp(session; structure = "VISam", stimulus = "spontaneous")
S = powerspectrum(LFP, 0.1)
# spectrumplot(S[:, 5])
# current_axis().limits = ((0.1, 500), nothing)
# current_figure()
# L = map(fooof, eachslice(ustripall(S), dims = (2)))
begin
    cfooof = x -> AN.aperiodicfit(x, [1, 300]; aperiodic_mode = "fixed", max_n_peaks = 10,
                                  peak_threshold = 1, peak_width_limits = [1, 50])
    i = 20
    s = collect(eachslice(ustripall(S), dims = (2)))[i][Freq(1 .. 500)]
    L, D = cfooof(s)
    f = Figure()
    ax = Axis(f[1, 1]; xscale = log10, yscale = log10)
    plot!(ax, freqs(s), (L.(freqs(s))))
    lines!(ax, freqs(s), parent(s))
    f
end
begin
    a = collect(eachslice(ustripall(S), dims = (2)))
    @time cfooof(s[1:5:end])
    χ = map(a) do s
        s = s[Freq(1 .. 500)][1:5:end]
        L, D = cfooof(s)
        χ = D[:χ]
    end
end
begin
    plot(χ)
end
begin
    m = LFP .- mean(LFP, dims = 1)
    m = dropdims(mean(abs.(m), dims = 1), dims = 1)
    plot(m)
end
begin
    f = Figure()
    ax = Axis(f[1, 1]; xscale = log10, yscale = log10)
    map(eachslice(S[2:end, 1:6:end]; dims = 2)) do s
        lines!(ax, decompose(s)...)
    end
    f
end
begin
    x = colorednoise(0.01:0.01:100)
    s = powerspectrum(x)
    L, D = fooof(s)
    f = Figure()
    ax = Axis(f[1, 1]; xscale = log10, yscale = log10)
    plot!(ax, freqs(s), (L.(freqs(s))), color = :red)
    lines!(ax, freqs(s), parent(s))
    f
end
