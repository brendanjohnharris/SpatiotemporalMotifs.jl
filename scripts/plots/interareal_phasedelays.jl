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
    Δ = map(Δ) do (a, b) # ! Check interpretation of phase difference to propagation direction
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
        ∂̄ = dropdims(mean(∂h, dims = Dim{:trial}), dims = :trial)
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

begin # * Thalamic comparison
    function load_thalamic_theta(sessionids = nothing;
                                 path = datadir("thalamus_calculations"),
                                 stimulus = r"Natural_Images")
        Q = calcquality(path)
        idxs = lookup(Q, :sessionid) .∈ [sessionids]
        Q = Q[sessionid = idxs, stimulus = At(stimulus)]
        structure = lookup(Q, :structure) |> only
        out = map(lookup(Q, :sessionid)) do sessionid
            D = @strdict structure sessionid stimulus
            data = load(joinpath(path, savepath(D, "jld2")), "ϕ")
            m = dropdims(circularmean(data, dims = 2), dims = 2)
        end
        sessionids = lookup(Q, :sessionid)
        out = DimArray(out, (Dim{:sessionid}(sessionids),))
        return out
    end
    thth = load_thalamic_theta(oursessions)
end

begin # * Check direction of phase difference
    thout = map(out) do o
        [_o[:ϕ] for _o in o if _o[:sessionid] in lookup(thth, :sessionid)]
    end
    thout = progressmap(thout; parallel = true) do o
        progressmap(o; parallel = true) do x
            y = thth[sessionid = At(metadata(x)[:sessionid])]
            I = intersect(Interval(x), Interval(y))
            x = x[Ti(I)]
            y = y[Ti(I)]
            @assert length(y) > 1500
            y = set(y, Ti => times(x))
            mapslices(x; dims = (1, 3)) do x
                phasegrad(y, x)
            end
        end
    end
end

begin # *
    bins = 0.05:0.1:0.95
    mth = map(thout) do x
        map(x) do y
            out = dropdims(mean(y[Ti(0.25u"s" .. 0.75u"s")], dims = (1, 3)); dims = (1, 3))
            out = mean.(SpatiotemporalMotifs.Bins(lookup(out, :depth); bins)(out))
        end
    end
    seshs = [metadata(x)[:sessionid] for x in thout[1]]
    mth = stack.([Dim{:sessionid}(seshs)], mth; dims = 1)
end

begin
    f = Figure()
    ax = Axis(f[1, 1])
    p = heatmap!(ax, thout[2][1][:, :, 120]; colormap = phasecolormap, colorrange = [-π, π])
    Colorbar(f[1, 2], p)
    f
end

begin # * depth wise differences
    f = Figure()
    ax = Axis(f[1, 1])
    for i in 1:6
        d = lookup(mth[i], :bin)
        μ, (l, u) = bootstrapmedian(mth[i]; dims = 1)
        if sum(μ .> 0) < 5
            μ = -μ
            l = -l
            u = -u
        end
        band!(ax, d, collect(l), collect(u), color = (structurecolors[i], 0.3))
        scatterlines!(ax, d, μ, label = structures[i], color = structurecolors[i])
    end
    axislegend(ax)
    f
end
begin # * Time and depth-wise differences. Should we just repeat the wavenumber/order parameter analysis, but for relative phases?? Ask.
