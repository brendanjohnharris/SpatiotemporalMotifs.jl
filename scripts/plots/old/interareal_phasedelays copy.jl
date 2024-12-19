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
using Distributions
set_theme!(foresight(:physics))
Random.seed!(42)

stimulus = r"Natural_Images"
vars = [:ϕ, :ω]

session_table = load(datadir("posthoc_session_table.jld2"), "session_table")
oursessions = session_table.ecephys_session_id

path = datadir("calculations")
Q = calcquality(path)[Structure = At(structures)]
out = load_calculations(Q; stimulus, vars)
out = map(out) do O
    filter(o -> (o[:sessionid] in oursessions), O)
end
datafile = datadir("interareal_phasedelays.jld2")

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

if !isfile(datafile)
    begin # * Functional hierarchy scores
        # ? See https://github.com/AllenInstitute/neuropixels_platform_paper/blob/master/Figure2/comparison_anatomical_functional_connectivity_final.ipynb
        dir = tempdir()
        baseurl = "https://github.com/AllenInstitute/neuropixels_platform_paper/raw/master/data/processed_data"
        f = "FFscore_grating_10.npy"
        Downloads.download(joinpath(baseurl, f), joinpath(dir, f))
        @assert SpatiotemporalMotifs.structures ==
                ["VISp", "VISl", "VISrl", "VISal", "VISpm", "VISam"] # The order of this FC matrix
        FF_score, FF_score_b, FF_score_std = eachslice(load(joinpath(dir, f)), dims = 1)
        FF_score = ToolsArray(FF_score', # Original matrix has higher-order structures on rows
                              (Dim{:structure_1}(SpatiotemporalMotifs.structures),
                               Dim{:structure_2}(SpatiotemporalMotifs.structures)))
    end
    begin # * Approximate 2D channel coordinates
        channels = AN.VisualBehavior.getchannels()
        idxs = .!ismissing.(channels.structure_acronym) .&
               (channels.structure_acronym .∈ [structures]) # All the probes are in one hemisphere
        channels = channels[idxs, :]
        idxs = .!isnan.(channels.left_right_ccf_coordinate) .&
               .!isnan.(channels.anterior_posterior_ccf_coordinate) .&
               .!isnan.(channels.dorsal_ventral_ccf_coordinate)
        channels = channels[idxs, :]
        # @assert all(channels.left_right_ccf_coordinate .> 5000)# All channels are in right hemisphere

        channels = DataFrames.groupby(channels, :ecephys_session_id)
        channels = [c for c in channels if length(unique(c.structure_acronym)) == 6]
        channels = [c for c in channels if only(unique(c.ecephys_session_id)) ∈ oursessions]

        for c in channels
            if mean(c.left_right_ccf_coordinate) .< 5000 # Left hemisphere, so reflect
                c.left_right_ccf_coordinate = .-c.left_right_ccf_coordinate
            end
            xs = hcat(c.left_right_ccf_coordinate,
                      c.anterior_posterior_ccf_coordinate,
                      c.dorsal_ventral_ccf_coordinate)
            xs = xs .- mean(xs, dims = 1)
            U, S, V = svd(xs)
            S[end] = 0.0
            ys = U * Diagonal(S) * V'

            begin # * Check the projection
                as = hcat(ys[:, 1], ys[:, 2], ones(size(ys, 1)))
                bs = ys[:, 3]
                pl = inv(as' * as) * as' * bs
                @assert all(pl[1] .* ys[:, 1] .+ pl[2] .* ys[:, 2] .+ pl[3] .≈ ys[:, 3])
            end

            P = V[:, 1:2]
            zs = xs * P
            lr_cor = cor(zs[:, 1], c.left_right_ccf_coordinate)
            @assert lr_cor > 0.7

            ap_cor = cor(zs[:, 2], c.anterior_posterior_ccf_coordinate)
            @assert ap_cor < -0.7
            @assert eigvals(cov(zs)) ≈ eigvals(cov(xs))[2:end]
            c.x = zs[:, 1]
            c.y = zs[:, 2]
        end

        # cpositions = map(eachrow(channels)) do c
        #     c.id => [
        #         c.left_right_ccf_coordinate,
        #         c.anterior_posterior_ccf_coordinate,
        #         c.dorsal_ventral_ccf_coordinate
        #     ]
        # end |> Dict
        channels = vcat(channels...)
    end

    begin # * Subject-by-subject phasedelays
        unidepths = range(0.05, 0.95, length = 19)
        pxy = map(out...) do o...
            ϕs = map(o) do _o
                ϕ = _o[:ϕ]
                ω = _o[:ω]
                odepths = lookup(ϕ, Depth)
                ϕ = ϕ[Depth(Near(unidepths))]
                ω = ω[Depth(Near(unidepths))]
                idxs = indexin(lookup(ϕ, Depth), odepths)
                cs = DimensionalData.metadata(ϕ)[:depths]
                sortidxs = sortperm(collect(values(cs))) # Assume perfect monotonic correlation to streamline depths. The depths in 'cs' are depths along the probe, smaller = more superficial
                cs = collect(keys(cs))[sortidxs][idxs]
                idxs = [channels.ecephys_channel_id .== c for c in cs]
                idxs = findfirst.(idxs)
                xs = channels[idxs, :x] .|> Float32
                ys = channels[idxs, :y] .|> Float32
                (ϕ, xs, ys, ω)
            end
            [getindex.(ϕs, i) for i in 1:4]
        end
        xs = getindex.(pxy, 2) # ! Maybe set this to a single mean value??
        xs = map(xs) do x
            map(x) do _x
                _x .= mean(_x)
            end
        end
        ys = getindex.(pxy, 3)
        ys = map(ys) do y
            map(y) do _y
                _y .= mean(_y)
            end
        end
        ϕs = getindex.(pxy, 1)
        ωs = getindex.(pxy, 4)
    end

    begin # * Use extra workers if we can
        if haskey(ENV, "JULIA_DISTRIBUTED") && length(procs()) == 1
            using USydClusters
            USydClusters.Physics.addprocs(8; mem = 16, ncpus = 4,
                                          project = projectdir())
            @everywhere using SpatiotemporalMotifs
            @everywhere SpatiotemporalMotifs.@preamble
        end
    end
    Δϕ = pmap(ϕs, ωs) do ϕ, ω # Takes about 5 mins over 32 cores, 180 Gb
        changetimes = lookup.(ϕ, :changetime)
        latency = [maximum(abs.(a .- b))
                   for (a, b) in zip(changetimes[2:end], changetimes[1:(end - 1)])]
        @assert all(latency .< 0.05u"s")
        changetimes = first(changetimes)
        ts = times.(ϕ)
        @assert abs(minimum(first.(ts)) - maximum(first.(ts))) < 0.016u"s"
        ϕ = map(ϕ, ω) do x, w
            x[w .< 0u"Hz"] .= NaN # ! Mask negative frequency periods
            x[1:1562, :, :] # Some trials are a tiny bit shorter
        end
        ϕ = set.(ϕ, [𝑡 => times(ϕ[1])])
        ϕ = set.(ϕ, [Depth => Depth(unidepths)])
        ϕ = set.(ϕ, [Dim{:changetime} => Dim{:changetime}(changetimes)])
        uniphi = stack(Structure(structures), ϕ)
        Δ = [(a, b) for a in eachslice(uniphi, dims = 4), b in eachslice(uniphi, dims = 4)]
        Δ = Δ[filter(!=(0), triu(LinearIndices(Δ), 1))]
        Δ = map(Δ) do (a, b) # ! Check interpretation of phase difference to propagation direction
            mod.(b .- a .+ π, eltype(b)(2π)) .- π
        end
        Δ = stack(Dim{:pair}(eachindex(Δ)), Δ, dims = 4)
        Δ = permutedims(Δ, (4, 1, 2, 3))
    end
    Δxs = pmap(xs) do X
        X = map(X...) do x...
            Δx = [b - a for a in x, b in x]
            Δx = Δx[filter(!=(0), triu(LinearIndices(Δx), 1))]
            Δx = ToolsArray(Δx, (Dim{:pair}(1:length(Δx)),))
        end
        stack(Depth(unidepths), X)
    end
    Δys = pmap(ys) do Y
        Y = map(Y...) do y...
            Δy = [b - a for a in y, b in y]
            Δy = Δy[filter(!=(0), triu(LinearIndices(Δy), 1))]
            Δy = ToolsArray(Δy, (Dim{:pair}(1:length(Δy)),))
        end
        stack(Depth(unidepths), Y)
    end
    Δxys = pmap(Δys, Δxs) do Δy, Δx
        maxs = maximum(sqrt.(Δy .^ 2 .+ Δx .^ 2); dims = :pair) # Normalize so all these measures are comparable. Note that this is constant across depths, due to our projection above
        Δy ./= maxs
        Δx ./= maxs
        return Δx, Δy
    end
    Δxs = first.(Δxys)
    Δys = last.(Δxys)
    @assert maximum(maximum.(Δxs)) ≤ 1
    @assert maximum(maximum.(Δys)) ≤ 1

    @assert structures == [DimensionalData.metadata(phi)[:structure] for phi in ϕs[1]]
    hs = getindex.([hierarchy_scores], structures) .|> Float32
    Δhs = [b - a for a in hs, b in hs] # Make a distance matrix
    Δhs = Δhs[filter(!=(0), triu(LinearIndices(Δhs), 1))]
    Δhs = Δhs ./ maximum(abs.(Δhs)) # Normalize to that the maximum distance between hierarchichal schores is 1
    Δhs = map(Δxs) do Δx # Copy onto shape of Δxs
        h = deepcopy(Δx)
        for _h in eachslice(h, dims = Depth)
            _h .= Δhs
        end
        return h
    end
    @assert maximum(maximum.(Δhs)) ≤ 1

    Δfs = FF_score[filter(!=(0), triu(LinearIndices(FF_score), 1))] # ? Δhs and Δfs match data used for siegle2021 figure 2f, without same-region FC. FF_score is already a 'distance' matrix.
    Δfs = Δfs ./ maximum(abs.(Δfs))
    Δfs = map(Δxs) do Δx
        f = deepcopy(Δx)
        for _f in eachslice(f, dims = Depth)
            _f .= Δfs
        end
        return f
    end
    @assert maximum(maximum.(Δfs)) ≤ 1

    ∂x = progressmap(Δϕ, Δxs; parallel = true) do Δ, Δx
        mapslices(Δ; dims = (:pair, Depth)) do Δ
            .-Δ ./ Δx # Minus because phase increases over time
        end
    end
    ∂y = progressmap(Δϕ, Δys; parallel = true) do Δ, Δy
        mapslices(Δ; dims = (:pair, Depth)) do Δ
            .-Δ ./ Δy # Minus because phase increases over time
        end
    end
    ∂h = progressmap(Δϕ, Δhs; parallel = true) do Δ, Δh
        mapslices(Δ; dims = (:pair, Depth)) do Δ
            .-sign.(Δ) ./ sign.(Δh) # Minus because phase increases over time
        end
    end
    ∂f = progressmap(Δϕ, Δfs; parallel = true) do Δ, Δf
        mapslices(Δ; dims = (:pair, Depth)) do Δ
            .-sign.(Δ) ./ sign.(Δf) # Minus because phase increases over time
        end
    end

    ∂ = progressmap(∂x, ∂y; parallel = true) do dx, dy # * Quite slow
        N = sqrt.(dx .^ 2 .+ dy .^ 2)
        dx = dx ./ N # Normalize x and y components to give a unit vector
        dy = dy ./ N
        dx = nansafe(mean, dims = :pair)(dx)
        dy = nansafe(mean, dims = :pair)(dy)
        m = sqrt.(dx .^ 2 .+ dy .^ 2)
        m = dropdims(m, dims = :pair)
    end
    ∂x = progressmap(∂x; parallel = true) do ∂x
        dropdims(nansafe(mean, dims = :pair)(∂x), dims = :pair) # x component of mean vector
    end
    ∂y = map(∂y) do ∂y
        dropdims(nansafe(mean, dims = :pair)(∂y), dims = :pair) # y component of mean vector
    end
    # ψ = map(∂x, ∂y) do ∂x, ∂y
    #     atan.(∂y ./ ∂x) # Angle of mean vector. Not so meaningful when the spatial coordinates are normalized
    # end
    ∂h = map(∂h) do ∂h
        dropdims(nansafe(mean, dims = :pair)(∂h), dims = :pair)
    end
    ∂f = map(∂f) do ∂f
        dropdims(nansafe(mean, dims = :pair)(∂f), dims = :pair)
    end

    begin # * Average quantities
        ∂x̄ = dropdims(mean(nansafe(mean, dims = :changetime).(∂x)); dims = :changetime)
        ∂ȳ = dropdims(mean(nansafe(mean, dims = :changetime).(∂y)); dims = :changetime)
        ∂̄ = dropdims(mean(nansafe(mean, dims = :changetime).(∂)); dims = :changetime)
        ∂h̄ = dropdims(mean(nansafe(mean, dims = :changetime).(∂h)); dims = :changetime)
        ∂f̄ = dropdims(mean(nansafe(mean, dims = :changetime).(∂f)); dims = :changetime)
    end

    begin
        tagsave(datafile, @strdict unidepths FF_score ∂x̄ ∂ȳ ∂̄ ∂h̄ ∂f̄)
    end
else
    @unpack unidepths, FF_score, ∂x̄, ∂ȳ, ∂̄, ∂h̄, ∂f̄, = load(datafile)
end

begin # * Plots
    f = FourPanel()
    gs = subdivide(f, 2, 2)
    layerints = load(datadir("grand_unified_layers.jld2"), "layerints")
    begin # * Schematic of calculation
        ax = Axis(gs[1], yreversed = true)#, aspect = DataAspect())
        hidedecorations!(ax)
        hidespines!(ax)
        outline, infill, l = plot_visual_cortex!(ax; colormap = structurecolors,
                                                 fillalpha = 0.2,
                                                 strokealpha = 0.8)
        delete!(l)
        plotstructurecenters!(ax)
        arrows!(ax, [-500, -500], [550, 550], [65, 0], [0, -80])
    end
    begin # * Correlation to hierarchy score
        ax = Axis(gs[2][1, 1], yreversed = true,
                  title = "θ propagation (hierarchy)", xlabel = "Time (s)",
                  ylabel = "Cortical depth (%)")
        # ∂̄ = dropdims(mean(∂h, dims = Trial), dims = Trial)
        p, _ = plotlayermap!(ax, ∂h̄[𝑡(SpatiotemporalMotifs.INTERVAL)] |> ustripall,
                             colormap = darksunset,
                             colorrange = symextrema(∂h̄))
        Colorbar(gs[2][1, 2], p; label = "Mean hierarchical ∇ (a.u.)")
        plotlayerints!(ax, layerints; flipside = true, newticks = false,
                       bgcolor = Makie.RGBA(0, 0, 0, 0))
        ax.limits = (nothing, (0.05, 0.95))
        f
    end
    begin # * Correlation to functional hierarchy score
        ax = Axis(gs[4][1, 1], yreversed = true,
                  title = "θ propagation (hierarchy)", xlabel = "Time (s)",
                  ylabel = "Cortical depth (%)")
        # ∂̄ = dropdims(mean(∂h, dims = Trial), dims = Trial)
        p, _ = plotlayermap!(ax, ∂f̄[𝑡(SpatiotemporalMotifs.INTERVAL)] |> ustripall,
                             colormap = darksunset,
                             colorrange = symextrema(∂f̄))
        Colorbar(gs[4][1, 2], p; label = "Mean hierarchical ∇ (a.u.)")
        plotlayerints!(ax, layerints; flipside = true, newticks = false,
                       bgcolor = Makie.RGBA(0, 0, 0, 0))
        ax.limits = (nothing, (0.05, 0.95))
        f
    end
    begin # * Correlation to position
        ax = Axis(gs[3][1, 1], yreversed = true, ytickformat = depthticks,
                  title = "θ propagation (position)", xlabel = "Time (s)",
                  ylabel = "Cortical depth (%)")
        # ∂̄ = dropdims(mean(∂[:, :, lookup(∂, Trial) .== true],
        #                    dims = Trial),
        #               dims = Trial)
        p, _ = plotlayermap!(ax, ∂ȳ[𝑡(SpatiotemporalMotifs.INTERVAL)] |> ustripall,
                             colormap = :inferno)
        Colorbar(gs[3][1, 2], p; label = "Mean positional ∇ (a.u.)")
        plotlayerints!(ax, layerints; flipside = true, newticks = false,
                       bgcolor = Makie.RGBA(0, 0, 0, 0))
        ax.limits = (nothing, (0.05, 0.95))
        f
    end
    # begin # * Angle over time
    #     ax = Axis(gs[4][1, 1], yreversed = true, ytickformat = depthticks,
    #               title = "θ propagation direction", xlabel = "Time (s)",
    #               ylabel = "Cortical depth (%)")
    #     # ∂̄ = dropdims(circularmean(ψ[:, :, lookup(ψ, Trial) .== true],
    #     #                            dims = Trial),
    #     #               dims = Trial)
    #     p, _ = plotlayermap!(ax, ψ̄[𝑡(SpatiotemporalMotifs.INTERVAL)] |> ustripall,
    #                          colormap = binarysunset, colorrange = symextrema(ψ̄))
    #     Colorbar(gs[4][1, 2], p; label = "Mean ψ (radians)")
    #     plotlayerints!(ax, layerints; flipside = true, newticks = false,
    #                    bgcolor = Makie.RGBA(0, 0, 0, 0))
    #     ax.limits = (nothing, (0.05, 0.95))
    #     f
    # end
    rowsize!(f.layout, 1, Relative(0.5))
    rowsize!(f.layout, 2, Relative(0.5))
    addlabels!(f)
    wsave(plotdir("interareal_phasedelays", "interareal_phasedelays.pdf"), f)
    f
end

if false
    begin # * Animate the phase at any given moment, interpolating the background,
        mphi = progressmap(ϕs; parallel = true) do ϕ
            ϕ = map(ϕ) do p
                set(set(dropdims(circularmean(p[1:1562, :, :]; dims = :changetime);
                                 dims = :changetime), 𝑡 => times(ϕs[1][1])[1:1562]),
                    Depth => Depth(unidepths))
            end
            stack(Structure(structures), ϕ |> collect)
        end
        mphi = cat(mphi...; dims = SessionID(oursessions))
        mphi = dropdims(circularmean(mphi; dims = SessionID); dims = SessionID)
    end
    begin
        set_theme!(foresight(:dark))
        d = 19
        # subphi = mphi[:, d, :]
        # subpar = ∂̄[:, d]
        # subh = ∂h̄[:, d]
        # subpsi = ψ̄[:, d]

        n = SpatiotemporalMotifs.DEFAULT_TRIAL_NUM
        subphi = map(ϕs[end]) do p
            set(set(p[1:1562, :, n], 𝑡 => times(ϕs[1][1])[1:1562]),
                Depth => Depth(unidepths))
        end
        subphi = cat(subphi...; dims = Structure(structures))
        subphi = subphi[:, d, :]
        subpar = ∂[end][:, d, n]
        subh = ∂h[end][:, d, n]
        subpsi = ψ[end][:, d, n]

        t = Observable(first(times(mphi)))
        c = lift(t -> subphi[𝑡(At(t))] |> collect, t)
        psi = lift(t -> subpsi[𝑡(At(t))] |> collect, t)
        layout = Point2f[(-350, 350), # VISp
                         (-170, 310), # VISl
                         (-300, 130), # VISrl
                         (-180, 195), # VISal
                         (-475, 240), # VISpm
                         (-450, 140)] # VISam
        o = mean(layout)
        psi1 = lift(psi -> [100.0 .* sin.(psi)], psi)
        psi2 = lift(psi -> [100.0 .* cos.(psi)], psi)
        f = Figure(; size = (600, 600))
        ax = Axis(f[1, 1]; yreversed = true)
        hidedecorations!(ax)
        hidespines!(ax)

        p1, p2, l = plot_visual_cortex!(ax, colors = c, colorrange = (-pi, pi),
                                        colormap = phasecolormap)
        delete!(l)

        scatter!(ax, layout; color = c, colormap = phasecolormap,
                 markersize = 80,
                 colorrange = (-pi, pi))
        text!(ax, layout; text = structures, align = (:center, :center),
              fontsize = 16, color = :black)
        f
    end
    begin # * Animate
        record(f, plotdir("interareal_phasedelays", "interareal_phasedelays.mp4"),
               times(mphi)[1:2:end]) do _t
            t[] = _t
        end
    end
end
