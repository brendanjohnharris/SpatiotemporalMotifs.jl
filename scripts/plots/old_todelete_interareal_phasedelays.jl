#! /bin/bash
#=
exec julia -t auto "${BASH_SOURCE[0]}" "$@"
=#
using DrWatson
@quickactivate "SpatiotemporalMotifs"

import DimensionalData: metadata
using SpatiotemporalMotifs
@preamble
import TimeseriesTools: freqs
using GraphMakie
using GraphMakie.Graphs
using SimpleWeightedGraphs
using LinearAlgebra
set_theme!(foresight(:physics))
Random.seed!(42)

stimulus = r"Natural_Images"
vars = [:ϕ]

session_table = load(datadir("session_table.jld2"), "session_table")
oursessions = session_table.ecephys_session_id

path = datadir("calculations")
Q = calcquality(path)[structure = At(structures)]

config = @strdict stimulus vars
data, file = produce_or_load(produce_out(Q), config, datadir(); filename = savepath,
                             prefix = "out")
out = data["out"]
data, file = produce_or_load(produce_uni, config, datadir(); filename = savepath,
                             prefix = "uni")
uni = data["uni"]

function cylindricalcor(alpha, x)

    # Compute correlation coefficients for sin and cos independently
    rxs = cor(x, sin.(alpha))
    rxc = cor(x, cos.(alpha))
    rcs = cor(sin.(alpha), cos.(alpha))

    # Compute circular-linear correlation
    rho = sqrt((rxc^2 + rxs^2 - 2 * rxc * rxs * rcs) / (1 - rcs^2))

    # Compute p-value
    # pval = 1 - cdf(Chisq(2), n * rho^2)

    return rho
end

begin # * Uniphi
    uniphi = getindex.(uni, :ϕ)
    @assert [dims(uniphi[1], 1) == dims(uniphi[i], 1) for i in eachindex(uniphi)] |> all
    @assert [dims(uniphi[1], 3) == dims(uniphi[i], 3) for i in eachindex(uniphi)] |> all
    unidepths = range(0.05, 0.95, length = 10)
    uniphi = map(uniphi) do phi
        set(phi[Depth(Near(unidepths))], Depth => Depth(unidepths))
    end
    uniphi = stack(Dim{:structure}(structures), uniphi)
    uniphi = uniphi[1:5:end, :, :, :] # Downsample
end

begin # * Massive
    xs = getindex.([visual_cortex_layout], lookup(uniphi, :structure)) .|> first .|> Float32
    ys = getindex.([visual_cortex_layout], lookup(uniphi, :structure)) .|> last .|> Float32
    hs = getindex.([hierarchy_scores], lookup(uniphi, :structure)) .|> Float32
    xs = MinMax(xs)(xs)
    ys = MinMax(ys)(ys)
    hs = MinMax(hs)(hs)

    Δ = [(a, b) for a in eachslice(uniphi, dims = 4), b in eachslice(uniphi, dims = 4)]
    Δ = Δ[filter(!=(0), triu(LinearIndices(Δ), 1))]
    Δ = map(Δ) do (a, b)
        mod.(b .- a .+ π, 2π) .- π
    end
    Δ = stack(Dim{:pair}(eachindex(Δ)), Δ, dims = 4)
    Δ = permutedims(Δ, (4, 1, 2, 3))

    Δx = [b - a for a in xs, b in xs]
    Δx = Δx[filter(!=(0), triu(LinearIndices(Δx), 1))]
    Δy = [b - a for a in ys, b in ys]
    Δy = Δy[filter(!=(0), triu(LinearIndices(Δy), 1))]
    Δh = [b - a for a in hs, b in hs]
    Δh = Δh[filter(!=(0), triu(LinearIndices(Δh), 1))]

    ∂x = Δ ./ Δx
    ∂y = Δ ./ Δy
    ∂h = Δ ./ Δh

    ∂x = dropdims(mean(∂x, dims = 1), dims = 1) # x component of mean vector
    ∂y = dropdims(mean(∂y, dims = 1), dims = 1) # y component of mean vector
    ∂ = sqrt.(∂x .^ 2 .+ ∂y .^ 2) # Norm of mean vector
    ψ = atan.(∂y ./ ∂x) # Angle of mean vector
    ∂h = dropdims(mean(∂h, dims = 1), dims = 1)
    ψ = dropdims(circularmean(ψ, dims = 1), dims = 1)
end

begin # * Plot mean correlation magnitude over layers and time
    ∂̄ = dropdims(circularmean(∂[:, :, lookup(∂, Trial) .== true], dims = Trial),
                  dims = Trial)
end
begin
    f = Figure()
    ax = Axis(f[1, 1], yreversed = true)
    p = heatmap!(ax, ∂̄[:, 3:end] |> ustripall, colormap = :viridis)
    Colorbar(f[1, 2], p)
    f
end

begin # * Calculate ρ_x, ρ_y, ρ_h for unified depths in each subject
    xs = getindex.([visual_cortex_layout], lookup(uniphi, :structure)) .|> first
    ys = getindex.([visual_cortex_layout], lookup(uniphi, :structure)) .|> last
    hs = getindex.([hierarchy_scores], lookup(uniphi, :structure))
    xs = MinMax(xs)(xs)
    ys = MinMax(ys)(ys)
    hs = MinMax(hs)(hs)

    ρ_x = similar(uniphi[structure = At("VISp")])
    ρ_y = similar(ρ_x)
    ρ_h = similar(ρ_x)
    Threads.@threads for (i, uphi) in collect(enumerate(eachslice(uniphi, dims = 1)))
        r = mapslices(uphi, dims = 3) do phi
            cylindricalcor(phi, xs)
        end
        ρ_x[i, :, :] .= dropdims(r, dims = 3)

        r = mapslices(uphi, dims = 3) do phi
            cylindricalcor(phi, ys)
        end
        ρ_y[i, :, :] .= dropdims(r, dims = 3)

        r = mapslices(uphi, dims = 3) do phi
            cylindricalcor(phi, hs)
        end
        ρ_h[i, :, :] .= dropdims(r, dims = 3)
    end
    ρ = sqrt.(ρ_x .^ 2 .+ ρ_y .^ 2)
    ψ = atan.(ρ_y ./ ρ_x)

    ρ_s = similar(ρ_x)
    Threads.@threads for (i, uphi) in collect(enumerate(eachslice(uniphi, dims = 3))) # Random surrogate for each trial
        rp = randperm(length(xs))
        r = mapslices(uphi, dims = 3) do phi
            cylindricalcor(phi, xs[rp])
        end
        ρ_s[:, :, i] .= dropdims(r, dims = 3)
    end
end

begin
    f = Figure()
    ax = Axis(f[1, 1], yreversed = true)
    p = heatmap!(ax, ρ_y[:, :, 100] |> ustripall, colormap = :turbo)
    Colorbar(f[1, 2], p)
    f
end

begin # * Plot mean correlation magnitude over layers and time
    ρ̄ = dropdims(mean(ρ_h[:, :, lookup(ρ, Trial) .== true], dims = Trial),
                  dims = Trial)
end
begin
    f = Figure()
    ax = Axis(f[1, 1], yreversed = true)
    p = heatmap!(ax, ρ̄ |> ustripall, colormap = darksunset,
                 colorrange = symextrema(ustripall(ρ̄)))
    Colorbar(f[1, 2], p)
    f
end

begin # * Plot mean correlation magnitude over layers and time
    ψ̄ = dropdims(circularmean(ψ[:, :, lookup(ρ, Trial) .== true], dims = Trial),
                  dims = Trial)
end
begin
    f = Figure()
    ax = Axis(f[1, 1], yreversed = true)
    p = heatmap!(ax, ψ̄ |> ustripall, colormap = cyclic)
    Colorbar(f[1, 2], p)
    f
end
