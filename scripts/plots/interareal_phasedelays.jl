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
        set(phi[Dim{:depth}(Near(unidepths))], Dim{:depth} => Dim{:depth}(unidepths))
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
end

begin # * Plots
    # !!! Currently for HIT trials
    f = FourPanel()
    gs = subdivide(f, 2, 2)

    begin # * Schematic of calculation
        ax = Axis(gs[1], yreversed = true, aspect = 1)
        hidedecorations!(ax)
        hidespines!(ax)
        outline, fill, l = plot_visual_cortex!(ax; colormap = structurecolors,
                                               fillalpha = 0.2,
                                               strokealpha = 0.8)
        delete!(l)
        plotstructurecenters!(ax)
        arrows!(ax, [100, 100], [550, 550], [65, 0], [0, -80])
    end
    begin # * Correlation to hierarchy score
        ax = Axis(gs[2][1, 1], yreversed = true,
                  title = "θ propagation (hierarchy)", xlabel = "Time (s)",
                  ylabel = "Cortical depth (%)")
        ∂̄ = dropdims(mean(∂h[:, :, lookup(∂h, :trial) .== true],
                           dims = Dim{:trial}),
                      dims = :trial)
        p, _ = plotlayermap!(ax, ∂̄[:, 3:end] |> ustripall, colormap = :viridis)
        Colorbar(gs[2][1, 2], p; label = "Mean hierarchical ∇ (a.u.)")
        f
    end
    begin # * Correlation to position
        ax = Axis(gs[3][1, 1], yreversed = true, ytickformat = depthticks,
                  title = "θ propagation (position)", xlabel = "Time (s)",
                  ylabel = "Cortical depth (%)")
        ∂̄ = dropdims(mean(∂[:, :, lookup(∂, :trial) .== true],
                           dims = Dim{:trial}),
                      dims = :trial)
        p, _ = plotlayermap!(ax, ∂̄[:, 3:end] |> ustripall, colormap = :viridis)
        Colorbar(gs[3][1, 2], p; label = "Mean positional ∇ (a.u.)")
        f
    end
    begin # * Angle over time
        ax = Axis(gs[4][1, 1], yreversed = true, ytickformat = depthticks,
                  title = "θ propagation direction", xlabel = "Time (s)",
                  ylabel = "Cortical depth (%)")
        ∂̄ = dropdims(circularmean(ψ[:, :, lookup(ψ, :trial) .== true],
                                   dims = Dim{:trial}),
                      dims = :trial)
        p, _ = plotlayermap!(ax, ∂̄[:, 3:end] |> ustripall,
                             colormap = reverse(cgrad(:viridis)),
                             colorrange = [0, maximum(∂̄[:, 3:end])])
        Colorbar(gs[4][1, 2], p; label = "Mean ψ (radians)")
        f
    end
    rowsize!(f.layout, 1, Relative(0.5))
    rowsize!(f.layout, 2, Relative(0.5))
    addlabels!(f)
    wsave(plotdir("interareal_phasedelays", "interareal_phasedelays.pdf"), f)
    f
end
