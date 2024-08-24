#! /bin/bash
#=
exec julia -t auto "${BASH_SOURCE[0]}" "$@"
=#
using DrWatson
@quickactivate "SpatiotemporalMotifs"

import TimeseriesTools: freqs
using FileIO
using Unitful
using ModulationIndices
using Normalization
using SpatiotemporalMotifs
@preamble
set_theme!(foresight(:physics))

stimulus = r"Natural_Images"
vars = [:ϕ, :r]

session_table = load(datadir("posthoc_session_table.jld2"), "session_table")
oursessions = session_table.ecephys_session_id

path = datadir("calculations")
Q = calcquality(path)[structure = At(structures)]
quality = mean(Q[stimulus = At(stimulus)])

begin
    config = @strdict stimulus vars
    data, file = produce_or_load(produce_out, config, datadir(); filename = savepath,
                                 prefix = "out")
    out = data["out"]
    GC.gc()
end

begin
    sessionids = getindex.(out[1], :sessionid)
    trials = [getindex.(o, :trials) for o in out]
    ϕs = [getindex.(o, :ϕ) for o in out]
    rs = [getindex.(o, :r) for o in out]
    out = []
    GC.gc()
    unidepths = range(0.05, 0.95, length = 19)
    pacc = map(ϕs, rs) do ϕ, r
        map(ϕ, r) do p, a
            pac(p, a; dims = (Ti, :changetime))
        end
    end
    GC.gc()
    pacc = map(pacc) do pa
        progressmap(pa; parallel = true) do p
            p = p[depth = Near(unidepths)]
            set(p, :depth => unidepths)
        end
    end
    pacc = stack.([Dim{:sessionid}(sessionids)], pacc; dims = 2)
    GC.gc()
end

begin # * Set up figure
    f = TwoPanel()
    gs = subdivide(f, 1, 2)
end

begin
    ax = Axis(gs[1])
    for (i, p) in enumerate(pacc)
        s = structures[i]
        # p = p ./ sum(p, dims = 1)
        μ, (σl, σh) = bootstrapmedian(p |> ustripall; dims = 2)

        mu = upsample(μ, 5)
        l = upsample(σl, 5) |> parent
        h = upsample(σh, 5) |> parent
        band!(ax, lookup(mu, 1), l, h; color = (structurecolormap[s], 0.2),
              label = s)
        lines!(ax, lookup(mu, 1), mu, color = structurecolormap[s], label = s,
               alpha = 0.8,
               linewidth = 4)
        scatter!(ax, lookup(μ, 1), μ, color = structurecolormap[s], label = s,
                 markersize = 10, alpha = 0.8)
    end
    f
end

begin # * Polar plot of coupling angle peaks (find peaks? what if there is more than 1? )
    function phipeak(r, ϕ; n = 50)
        ϕ = mod2pi.(ϕ .+ pi) .- pi
        ϕ = ModulationIndices.tortbin(ϕ; n)
        h = [mean(r[ϕ .== i]) for i in 1:n]
        angles = range(start = -π + π / n, stop = π - π / n, length = n) |> collect
        _, i = findmax(h)
        phimax = angles[i]
    end
    function phipeaks(r, ϕ; kwargs...)
        peaks = progressmap(eachslice(r, dims = :depth),
                            eachslice(ϕ, dims = :depth); parallel = true) do a, b
            phipeak(a, b; kwargs...)
        end
        out = ϕ[1, :, 1] # Just a depth slice
        out .= peaks
        return out
    end
    peaks = map(rs, ϕs) do r, ϕ
        progressmap(r, ϕ; parallel = true) do a, p
            phipeaks(a, p; n = 20)
        end
    end
    peaks = map(peaks) do pa
        map(pa) do p
            p = p[depth = Near(unidepths)]
            set(p, :depth => unidepths)
        end
    end
    peaks = stack.([Dim{:sessionid}(sessionids)], peaks; dims = 2)
end

begin # * Phase annotations
    μs = map(pacc) do p
        μ, (σl, σh) = bootstrapmedian(p |> ustripall; dims = 2)
    end .|> first
    if true
        idx = Depth(Near(0.25))
        alphamin = 0.2
        struc = 1
        tortinset!(gs[1], peaks[struc][idx],
                   colormap = seethrough(structurecolors[struc], alphamin, 1),
                   halign = 0.3, valign = 0.6, color = structurecolors[struc])
        scatter!(ax, [idx.val.val], [μs[struc][idx]], color = structurecolors[struc],
                 markersize = 15,
                 strokecolor = :white, strokewidth = 2)

        idx = Depth(Near(0.3))
        struc = 2
        tortinset!(gs[1], peaks[struc][idx],
                   colormap = seethrough(structurecolors[struc], alphamin, 1),
                   halign = 0.01, valign = 0.4, color = structurecolors[struc])
        scatter!(ax, [idx.val.val], [μs[struc][idx]], color = structurecolors[struc],
                 markersize = 15,
                 strokecolor = :white, strokewidth = 2)

        idx = Depth(Near(0.9))
        struc = 6
        tortinset!(gs[1], peaks[struc][idx],
                   colormap = seethrough(structurecolors[struc], alphamin, 1),
                   halign = 0.87, valign = 0.95, color = structurecolors[struc])
        scatter!(ax, [idx.val.val], [μs[struc][idx]], color = structurecolors[struc],
                 markersize = 15,
                 strokecolor = :white, strokewidth = 2)

        idx = Depth(Near(0.9))
        struc = 3
        tortinset!(gs[1], peaks[struc][idx],
                   colormap = seethrough(structurecolors[struc], alphamin, 1),
                   halign = 0.9, valign = 0.48, color = structurecolors[struc])
        scatter!(ax, [idx.val.val], [μs[struc][idx]], color = structurecolors[struc],
                 markersize = 15,
                 strokecolor = :white, strokewidth = 2)
    end
end

begin
    ax = PolarAxis(gs[2]; theta_as_x = false, thetalimits = (-0.1pi, 1.2pi),
                   rticks = 0:0.25:1, rtickformat = depthticks)
    for (i, p) in enumerate(peaks)
        s = structures[i]
        # p = p ./ sum(p, dims = 1)
        # μ, (σl, σh) = bootstrapmedian(p |> ustripall; dims = 2)
        μ, (σl, σh) = bootstrapaverage(circularmean, p |> ustripall; dims = 2)
        μc, _ = bootstrapmedian(pacc[i] |> ustripall; dims = 2)

        x = unwrap(μ)
        x = upsample(x, 5)
        mu = SpatiotemporalMotifs.wrap.(x; domain = (-π, π))

        x = unwrap(σl)
        x = upsample(x, 5)
        l = SpatiotemporalMotifs.wrap.(x; domain = (-π, π))

        x = unwrap(σh)
        x = upsample(x, 5)
        h = SpatiotemporalMotifs.wrap.(x; domain = (-π, π))

        muc = upsample(μc, 5)
        c = seethrough(structurecolormap[s])
        # band!(ax, lookup(mu, 1), l, h; color = muc |> collect, colormap = c, label = s)
        lines!(ax, lookup(mu, 1), mu, color = muc |> collect, label = s, colormap = c,
               linewidth = 7)
    end
    f
end

begin # * Single-trial PAC
    trialpac = map(ϕs, rs) do ϕ, r
        progressmap(ϕ, r; parallel = true) do p, a
            pac(p, a; dims = Ti)
        end
    end
    GC.gc()
    trialpac = map(trialpac) do pa
        progressmap(pa; parallel = true) do p
            p = p[depth = Near(unidepths)]
            set(p, :depth => unidepths)
        end
    end
    trialpac = map(trialpac, trials) do pa, tr
        map(pa, tr) do p, t
            @assert all(isapprox.(ustripall(lookup(p, :changetime)),
                                  t.change_time_with_display_delay; atol = 1e-1))
            set(p, :changetime => Dim{:trial}(t.hit))
        end
    end
end

begin # * Classification
    regcoef = 0.01
    folds = 5
    repeats = 10
    # bac = progressmap(pacc) do p
    #     classify_kfold(p; regcoef, k = folds, repeats)
    # end
    # p = vcat(pacc...)
    p = pacc[6]
    # p = ToolsArray(p, (Dim{:d}(1:size(p, 1)), dims(pacc[1], 2)))
    H = [getindex.(trialpac, i) for i in eachindex(trialpac[1])]
    H = map(H) do h
        trials = lookup(h[1], :trial)
        h = vcat(h...)
        h = ToolsArray(h, (Dim{:d}(1:size(h, 1)), Dim{:trial}(trials)))
        h = h[1:5:end, :]
    end
    bac = classify_kfold.(H; regcoef = 0.1, k = folds, repeats)

    # x = dropdims(mean(p[:, lookup(p, :trial)]; dims = 2); dims = 2)
    # x = (x .- dropdims(mean(p[:, .!lookup(p, :trial)]; dims = 2); dims = 2)) ./ x
    # plot(x)
end

begin
    f = Figure()
    ax = Axis(f[1, 1])
    lines!.(ax, [p for p in pacc[1][5:5:end]])
    f
end

begin # * For one subject, plot the layerwise PAC
    i = 5
    ϕ = ϕs[i]
    r = rs[i]
    ϕ = map(ϕ) do p
        d = dims(p, 2)
        p = vcat(eachslice(p, dims = 3)...)
        p = ToolsArray(p, (𝑡(1:size(p, 1)), d))
    end
    r = map(r) do p
        d = dims(p, 2)
        p = vcat(eachslice(p, dims = 3)...)
        p = ToolsArray(p, (𝑡(1:size(p, 1)), d))
    end

    pacc = [pac(ϕ[i], r[i]; dims = 1) for i in eachindex(ϕ)]
    # outs = map(eachindex(pacc)) do i
    #     s = structures[i]
    #     μ, (σl, σh) = bootstrapmedian(pacc[i] |> ustripall; dims = 2)
    # end
    # μs = first.(outs)
    # σs = last.(outs)

end

begin
    f = Figure()
    ax = Axis(f[1, 1])
    for i in eachindex(pacc)
        lines!(ax, lookup(pacc[i], 1), pacc[i], color = structurecolormap[structures[i]],
               label = structures[i],
               alpha = 0.8,
               linewidth = 4)
    end
    f
end
begin
    f = Figure()
    ax = Axis(f[1, 1], limits = ((0, 1), (nothing, nothing)),
              xlabel = "Cortical depth (%)", ylabel = "PAC", backgroundcolor = :transparent,
              xtickformat = depthticks)
    for i in eachindex(pacc)
        s = structures[i]
        μ = μs[i]
        (σl, σh) = σs[i]
        mu = upsample(μ, 5)
        l = upsample(σl, 5) |> parent
        h = upsample(σh, 5) |> parent
        band!(ax, lookup(mu, 1), l, h; color = (structurecolormap[s], 0.2),
              label = s)
        lines!(ax, lookup(mu, 1), mu, color = structurecolormap[s], label = s,
               alpha = 0.8,
               linewidth = 4)
        scatter!(ax, lookup(μ, 1), μ, color = structurecolormap[s], label = s,
                 markersize = 10, alpha = 0.8)
    end
    axislegend(ax; position = :lt, nbanks = 2, merge = true)
    f
end
