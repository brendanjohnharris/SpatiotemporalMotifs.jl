#! /bin/bash
#=
exec julia -t auto "${BASH_SOURCE[0]}" "$@"
=#
using DrWatson
@quickactivate "SpatiotemporalMotifs"

import TimeseriesTools: freqs
using FileIO
using Unitful
using Images
import CairoMakie.Axis
using SpatiotemporalMotifs
import SpatiotemporalMotifs: HistBins
using Peaks
@preamble
set_theme!(foresight(:physics))

stimulus = r"Natural_Images"
vars = [:r]
α = 0.9

session_table = load(datadir("posthoc_session_table.jld2"), "session_table")
oursessions = session_table.ecephys_session_id

begin # * Extract burst mask from each trial. Takes about 15 minutes on 32 cores
    thr = 2.0
    file = datadir("gamma_bursts.jld2")
    layerints = load(datadir("grand_unified_layers.jld2"), "layerints")
    if isfile(file)
        m, Δt, Δx, ints, schemr = load(file, "m", "Δt", "Δx", "ints", "schemr")
    else
        path = datadir("calculations")
        Q = calcquality(path)[Structure = At(structures)]
        quality = mean(Q[stimulus = At(stimulus)])
        config = @strdict stimulus vars
        data, _ = produce_or_load(produce_out(Q), config, datadir(); filename = savepath,
                                  prefix = "out")
        out = data["out"]
        out = map(out) do O # * Filter to posthoc sessions
            filter(o -> (o[:sessionid] in oursessions), O)
        end
        idx = getindex.(out[1], :sessionid) .== SpatiotemporalMotifs.DEFAULT_SESSION_ID
        schemr = only(out[6][idx])[:r][:, :, SpatiotemporalMotifs.DEFAULT_TRIAL_NUM]

        unidata, _ = produce_or_load(produce_uni, config, datadir(); filename = savepath)
        uni = unidata["uni"]
        @assert uni[1][:oursessions] == oursessions
        ints = getindex.(uni, :layerints)
        Bs = map(out) do O
            mtx = map(O) do o
                r = o[:r] .* 1000 # V to mV
                r = HalfZScore(r, dims = [1, 3])(r) # Normalized over time and trials
                pks, proms = findpeaks(r; minprom = thr)
                m = maskpeaks(r; minprom = thr) .> 0

                if false && metadata(r)[:sessionid] == 1048189115
                    ff = TwoPanel()
                    ax = Axis(ff[1, 1]; title = "Gamma amplitude", xlabel = "Time (s)",
                              ylabel = "Depth (μm)")
                    heatmap!(ax, r[:, :, 5] |> ustripall)
                    ax = Axis(ff[1, 2]; title = "Burst mask", xlabel = "Time (s)",
                              ylabel = "Depth (μm)")
                    heatmap!(ax, m[:, :, 5] |> ustripall)
                    display(ff)
                end

                Δt = mapslices(m; dims = 1) do m
                    @assert m isa AbstractVector
                    @assert eltype(m) <: Bool
                    length(m) < 1000 && @warn "Input seems too short to be a time series..."
                    m = Matrix(m') # label_components segfaults with vectors
                    cs = label_components(m)
                    c = component_lengths(cs) |> last
                    c = c * samplingperiod(r)
                end

                # * Replace unified depths with probe depths
                origdepths = sort(metadata(r)[:depths] |> values |> collect)
                _m = set(m, Depth => origdepths)
                _m = rectify(_m, dims = Depth)
                # m = m[1:3:end, :, :] # Downsample for speed

                Δx = mapslices(m; dims = 2) do m
                    @assert m isa AbstractVector
                    @assert eltype(m) <: Bool
                    length(m) > 30 && @warn "Input seems too long to be layer-wise..."
                    m = Matrix(m') # label_components segfaults with vectors
                    cs = label_components(m)
                    c = component_lengths(cs) |> last
                    c = c * step(lookup(_m, Depth)) .* u"μm"
                end
                return m, Δt, Δx
            end
            m, Δt, Δx = [getindex.(mtx, i) for i in eachindex(mtx[1])]
            return m, Δt, Δx
        end
        m, Δt, Δx = [getindex.(Bs, i) for i in eachindex(Bs[1])]
        tagsave(file, @strdict m Δt Δx ints schemr)
    end
end

begin # * Setup plot
    f = FourPanel()
    gs = subdivide(f, 2, 2)
end

begin # * Schematic diagram
    ax = Axis3(gs[1]; perspectiveness = 0.25, viewmode = :stretch, xlabel = "Time (ms)",
               ylabel = "Cortical depth (%)",
               ytickformat = depthticks,
               title = "γ-burst detection")
    ax.zspinesvisible = false
    ax.zgridvisible = false
    hidezdecorations!(ax)
    ax.ygridvisible = false
    ax.xgridvisible = false
    ax.xypanelvisible = true
    ax.xypanelcolor = (:black, 0.05)
    s = schemr[𝑡(0.12u"s" .. 0.25u"s")] |> ustripall
    s = set(s, 𝑡 => (lookup(s, 𝑡) .- minimum(lookup(s, 𝑡))) .* 1000)
    # s = upsample(s, 25, 1)
    # s = reverse(s, dims = 2)
    x = [lookup(s, 𝑡) for i in 1:size(s, 2)]
    y = [repeat([i], size(s, 1)) for i in lookup(s, Depth)]
    z = eachcol(parent(s))
    cmap = cgrad(layercolors)
    surface!(ax, s, alpha = 0.3; colormap = cmap, rasterize = 5)
    map(1.0 .- (eachindex(x) ./ length(x)), x, y, z) do i, x, y, z
        pks, vals = findmaxima(z, 10)

        pks, proms = peakproms(pks, z; minprom = 5e-5)

        lines!(ax, x, y, z; linewidth = 3, color = z, alpha = 0.1, colormap = cmap,
               colorrange = extrema(s), rasterize = 5)

        if !isempty(pks)
            pks, widths, leftedge, rightedge = peakwidths(pks, z, proms)

            scatter!(ax, x[pks], y[pks], z[pks]; markersize = 10, color = z[pks],
                     colormap = cmap, colorrange = extrema(s))

            # l = findmin(abs.(x .- leftedge)) |> last

            # r = findmin(abs.(x .- rightedge)) |> last
            # scatter!(ax, [x[l]], y[[l]], z[[l]]; markersize = 10, color = z[pks],
            #          colormap = cmap, colorrange = extrema(s))
            # scatter!(ax, [x[r]], y[[r]], z[[r]]; markersize = 10, color = z[pks],
            #          colormap = cmap, colorrange = extrema(s))
        end
    end

    linesegments!(ax, [(80, 0.5, 4e-5), (80, 1, 4e-5)]; linewidth = 5,
                  color = crimson)
    text!(ax, ((80, 0.5, 4e-5) .+ (80, 1, 4e-5)) ./ 2; text = "Width",
          offset = (20.0, 15.0))
    linesegments!(ax, [(30, 0.8, 4e-5), (60, 0.8, 4e-5)]; linewidth = 5,
                  color = cornflowerblue)
    text!(ax, ((30, 0.8, 4e-5) .+ (60, 0.8, 4e-5)) ./ 2; text = "Duration",
          align = (:right, :center),
          offset = (-40.0, 5.0))
    tightlimits!(ax)
    ax.azimuth = -1.3π / 3
    f
end

begin # * Heatmap of burst likelihood
    function pf!(g, i; compact = false, kwargs...)
        vm = m[i] # * VISp
        structure = metadata(vm[1])[:structure]
        ts = -0.25u"s" .. 0.75u"s"
        unidepths = 0.05:0.1:0.95
        ax = Axis(g[1, 1]; title = "Burst likelihood ($structure)", xlabel = "Time (s)",
                  yreversed = true,
                  limits = (extrema(ts |> ustripall), extrema(unidepths)))

        if compact
            if i[1] < 3
                hidexdecorations!(ax)
            end
        end
        Nb = map(vm) do m
            m = dropdims(mean(m .> 0, dims = 3), dims = 3) # * Probability of bursting, over trials
            m = m[𝑡(ts), Depth(Near(unidepths))]
            m = set(m, Depth => unidepths)
        end
        N = mean(Nb)
        p, ps = plotlayermap!(ax, N; colorrange = (0, 0.35), rasterize = 5)

        if !compact || i in [2, 4, 6]
            Colorbar(g[1, 2], p)
        end
        plotlayerints!(ax, ints[i])
        return p
    end
    pf!(gs[2], 1) # * Plot VISp for main text
    current_figure()

    sf = FourPanel()
    gsf = subdivide(sf, 3, 2)
    pf!.(gsf[:], eachindex(gsf); compact = true)
    addlabels!(sf)
    wsave(plotdir("gamma_bursts_task", "supplemental_burst_likelihood.pdf"), sf)
end

begin # * Distribution of burst durations
    bins = range(0, 1, length = 11)
    tbins = map(Δt) do tt
        ts = map(tt) do t
            B = HistBins(lookup(t, Depth); bins)
            b = B.(eachslice(t, dims = 3))
            b = map(x -> mean.(x), b)
            dropdims(mean(cat(b..., dims = 2), dims = 2), dims = 2)
        end
        hcat(ts...) # * Depth × trials
    end
end
begin # * Plot durations
    ax = Axis(gs[3]; title = "Burst duration", xlabel = "Cortical depth (%)",
              xtickformat = depthticks,
              ytickformat = x -> string.(round.(Int, x .* 1000)),
              ylabel = "Duration (ms)", limits = ((0, 1), (0.025, 0.08)))
    map((reverse ∘ collect ∘ enumerate)(structures)) do (i, s)
        μ, (σl, σh) = bootstrapmedian(tbins[i] |> ustripall; dims = 2)
        mu = upsample(μ, 10)
        l = upsample(σl, 10) |> parent
        h = upsample(σh, 10) |> parent
        band!(ax, lookup(mu, 1), l, h; color = (structurecolormap[s], 0.2), label = s)
        lines!(ax, lookup(mu, 1), mu, color = structurecolormap[s], label = s, alpha = α)
        scatter!(ax, lookup(μ, 1), μ, color = structurecolormap[s], label = s, alpha = α)
    end

    l = axislegend(ax; merge = true, nbanks = 2, position = :rt, framevisible = true)
    reverselegend!(l)
    plotlayerints!(ax, layerints; newticks = false, flipside = false, axis = :x)
    f
end
begin # * Distribution of burst widths
    bins = range(-0.25u"s", 0.75u"s", length = 25)
    xbins = map(Δx) do xx
        xs = map(xx) do x
            B = HistBins(lookup(x, 𝑡); bins)
            b = B.(eachslice(x, dims = 3))
            b = map(x -> mean.(x), b)
            dropdims(mean(cat(b..., dims = 2), dims = 2), dims = 2)
        end
        hcat(xs...) # * Time × trials
    end
end
begin # * Plot widths
    ax = Axis(gs[4]; title = "Burst width", xlabel = "Time (s)",
              ylabel = "Width (μm)", limits = ((-0.25, 0.75), (nothing, nothing)))
    vlines!(ax, [0.0, 0.25]; color = (:black, 0.5), linestyle = :dash, linewidth = 3)
    map((reverse ∘ collect ∘ enumerate)(structures)) do (i, s)
        μ, (σl, σh) = bootstrapmedian(xbins[i] |> ustripall; dims = 2)
        mu = upsample(μ, 10)
        l = upsample(σl, 10) |> parent
        h = upsample(σh, 10) |> parent
        band!(ax, lookup(mu, 1), l, h; color = (structurecolormap[s], 0.2), label = s)
        lines!(ax, lookup(mu, 1), mu, color = structurecolormap[s], label = s, alpha = α,
               linewidth = 4)
        scatter!(ax, lookup(μ, 1), μ, color = structurecolormap[s], label = s,
                 markersize = 10, alpha = α)
    end

    l = axislegend(ax; merge = true, nbanks = 2, position = :lt, framevisible = true)
    reverselegend!(l)
    f
end

begin # * Final adjustments
    addlabels!(f)
    wsave(plotdir("gamma_bursts_task", "gamma_bursts_task.pdf"), f)
end

# begin # * Supplemental material: average phase velocity maps in each region
#     f = Figure(size = (720, 1440))

#     for i in eachindex(uni)
#         k = uni[i][:aᵧ] .|> abs
#         k = k .* u"V"
#         k = uconvert.(u"μV", k)

#         # * Hit
#         ax = Axis(f[i, 1], yreversed = true)
#         structure = DimensionalData.metadata(k)[:structure]
#         ax.title = structure * ": hit"
#         m = median(k[:, :, lookup(k, Trial) .== true], dims = Trial)
#         m = dropdims(m, dims = Trial)
#         colorrange = maximum(abs.(ustripall(m))) * [0, 1]
#         ints = uni[i][:layerints]
#         p = plotlayermap!(ax, m, ints; arrows = false, colorrange, colormap = :bone) |>
#             first
#         if i > 5
#             ax.xlabel = "Time (s)"
#         end

#         # * Miss
#         ax = Axis(f[i, 2], yreversed = true)
#         ax.title = structure * ": miss"
#         m = median(k[:, :, lookup(k, Trial) .== false], dims = Trial)
#         m = dropdims(m, dims = Trial)
#         ints = uni[i][:layerints]
#         p = plotlayermap!(ax, m, ints; arrows = false, colorrange, colormap = :bone) |>
#             first
#         c = Colorbar(f[i, 3], p)
#         c.label = "Average 𝑟 ($(unit(eltype(k))))"
#     end
#     display(f)
#     wsave(plotdir("gamma_bursts_task", "supplemental_amplitudes.pdf"), f)
# end
