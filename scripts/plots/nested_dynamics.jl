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
using ModulationIndices
import CairoMakie.Axis
using SpatiotemporalMotifs
import SpatiotemporalMotifs: HistBins, structures
using Peaks
@preamble
set_theme!(foresight(:physics))

session_table = load(datadir("posthoc_session_table.jld2"), "session_table")
oursessions = session_table.ecephys_session_id

begin # * Set up master plot
    mf = SixPanel()
    mgs = subdivide(mf, 3, 2)
end

begin # * Burst masks and schematic
    stimulus = r"Natural_Images"
    vars = [:r]
    α = 0.9

    begin # * Setup plot
        f = TwoPanel()
        gs = subdivide(f, 1, 2)
    end

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
            out = load_calculations(Q; stimulus, vars)
            out = map(out) do O # * Filter to posthoc sessions
                filter(o -> (o[:sessionid] in oursessions), O)
            end
            idx = getindex.(out[1], :sessionid) .== SpatiotemporalMotifs.DEFAULT_SESSION_ID
            schemr = only(out[6][idx])[:r][:, :, 26] # Select some trial, 26

            uni = load_uni(; stimulus, vars)
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
                        length(m) < 1000 &&
                            @warn "Input seems too short to be a time series..."
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

    begin # * Schematic diagram
        ax = Axis3(mgs[1]; perspectiveness = 0.25, viewmode = :stretch,
                   xlabel = "Time (ms)",
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
        cmap = cgrad(:inferno)
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
        # pf!(gs[2], 1) # * Plot VISp for main text
        # current_figure()

        sf = SixPanel()
        gsf = subdivide(sf, 3, 2)
        pf!.(gsf[:], eachindex(gsf); compact = true)
        addlabels!(sf, labelformat)
        wsave(plotdir("nested_dynamics", "supplemental_burst_likelihood.pdf"), sf)
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
        ax = Axis(gs[1]; title = "Burst duration", xlabel = "Cortical depth (%)",
                  xtickformat = depthticks,
                  ytickformat = x -> string.(round.(Int, x .* 1000)),
                  ylabel = "Duration (ms)", limits = ((0, 1), (0.025, 0.08)))
        map((reverse ∘ collect ∘ enumerate)(structures)) do (i, s)
            μ, (σl, σh) = bootstrapmedian(tbins[i] |> ustripall; dims = 2)
            mu = upsample(μ, 10)
            l = upsample(σl, 10) |> parent
            h = upsample(σh, 10) |> parent
            band!(ax, lookup(mu, 1), l, h; color = (structurecolormap[s], 0.2), label = s)
            lines!(ax, lookup(mu, 1), collect(mu), color = structurecolormap[s], label = s,
                   alpha = α)
            scatter!(ax, lookup(μ, 1), collect(μ), color = structurecolormap[s], label = s,
                     alpha = α)
        end

        l = axislegend(ax; merge = true, nbanks = 2, position = :rt, framevisible = true,
                       labelsize = 10)
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
        ax = Axis(mgs[2]; title = "Burst width", xlabel = "Time (s)",
                  ylabel = "Width (μm)", limits = ((-0.25, 0.75), (nothing, nothing)))
        vlines!(ax, [0.0, 0.25]; color = (:black, 0.5), linestyle = :dash, linewidth = 3)
        map((reverse ∘ collect ∘ enumerate)(structures)) do (i, s)
            μ, (σl, σh) = bootstrapmedian(xbins[i] |> ustripall; dims = 2)
            mu = upsample(μ, 10)
            l = upsample(σl, 10) |> parent
            h = upsample(σh, 10) |> parent
            band!(ax, lookup(mu, 1), l, h; color = (structurecolormap[s], 0.2), label = s)
            lines!(ax, lookup(mu, 1), collect(mu), color = structurecolormap[s], label = s,
                   alpha = α,
                   linewidth = 4)
            scatter!(ax, lookup(μ, 1), collect(μ), color = structurecolormap[s], label = s,
                     markersize = 10, alpha = α)
        end

        l = axislegend(ax; merge = true, nbanks = 2, position = :lt, framevisible = true,
                       labelsize = 10)
        reverselegend!(l)
    end

    begin # * Final adjustments
        addlabels!(f, labelformat)
        wsave(plotdir("nested_dynamics", "supplemental_durations.pdf"), f)
    end # !! Also compare to surrogate durations?

    begin # * Write out stats
        statsfile = plotdir("nested_dynamics", "nested_dynamics.txt")
        close(open(statsfile, "w")) # Create the file or clear it
        open(statsfile, "a+") do file
            write(file, "\n# Burst durations\n")
            write(file, "Mean duration (s) = $(mean(vcat(parent.(tbins)...)))\n")
            write(file, "Duration std (s) = $(std(vcat(parent.(tbins)...)))\n")
        end
        open(statsfile, "a+") do file
            write(file, "\n# Burst widths\n")
            write(file, "Mean width (μm) = $(mean(vcat(parent.(xbins)...)))\n")
            write(file, "Width std (s) = $(std(vcat(parent.(xbins)...)))\n")
        end
        for T in [0u"s" .. 0.25u"s", 0.25u"s" .. 0.5u"s"]
            N = 1e6
            open(statsfile, "a+") do file
                write(file, "\n# Burst widths $T\n")
                subxbins = getindex.(xbins, [T], [:])
                x = getindex.([SpatiotemporalMotifs.hierarchy_scores], structures)
                _y = stack(Structure(structures), subxbins)
                ds = dims(_y) |> collect
                ds[2] = SessionID(oursessions)
                _y = rebuild(_y; dims = Tuple(ds))
                _y = dropdims(mean(_y, dims = :bin); dims = :bin)
                y = stack(Depth([0]), [_y])

                write(file, "Mean width (μm) = $(mean(vcat(subxbins...)))\n")
                write(file, "Width std (s) = $(std(vcat(subxbins...)))\n")

                # * Group level
                μ, σ, 𝑝 = hierarchicalkendall(x, ustrip.(y), :group; N) .|> first
                write(file, "\n## Group level")
                write(file, "\nmedian τ = $μ")
                write(file, "\n95\\% conf. = $σ")
                write(file, "\n𝑝 = $𝑝")
                write(file, "\nN = $N")
                write(file, "\n")

                # * Individual level
                μ, σ, 𝑝 = hierarchicalkendall(x, ustrip.(y), :individual; N) .|> first
                write(file, "\n## Individual level")
                write(file, "\nmedian τ = $μ")
                write(file, "\nIQR = $(σ[2] - σ[1])")
                write(file, "\n𝑝 = $𝑝")
                write(file, "\n")
            end
        end
    end
end

function phipeak(r, ϕ; n = 20)
    ϕ = mod2pi.(ϕ .+ pi) .- pi
    ϕ = ModulationIndices.tortbin(ϕ; n)
    h = [mean(r[ϕ .== i]) for i in 1:n]
    angles = range(start = -π + π / n, stop = π - π / n, length = n) |>
             collect
    _, i = findmax(h)
    phimax = angles[i]
end

begin # * Global and spatiotemporal PAC
    vars = [:ϕ, :r]

    path = datadir("calculations")
    Q = calcquality(path)[Structure = At(structures)]
    Q = Q[SessionID = At(oursessions)]
    quality = mean(Q[stimulus = At(stimulus)])

    stimuli = [r"Natural_Images", "spontaneous", "flash_250ms"]
    pQ = calcquality(datadir("power_spectra"))
    for stimulus in stimuli # * Average comodulograms
        _Q = pQ[stimulus = At(stimulus), Structure = At(structures)]
        subsessions = intersect(oursessions, lookup(_Q, SessionID))
        if length(subsessions) < length(oursessions)
            @warn "Power spectra calculations are incomplete, proceeding regardless"
        end
        _Q = _Q[SessionID(At(subsessions))]
        filebase = stimulus == "spontaneous" ? "" : "_$stimulus"
        f = Figure(size = (900, 1080))

        begin # * Load data
            S = map(lookup(_Q, Structure)) do structure
                out = map(lookup(_Q, SessionID)) do sessionid
                    if _Q[SessionID = At(sessionid), Structure = At(structure)] == 0
                        return nothing
                    end
                    filename = savepath((@strdict sessionid structure stimulus), "jld2",
                                        datadir("power_spectra"))
                    C = load(filename, "C")
                    # S = load(filename, "sC")
                    # return (C .- S) ./ median(S)
                end
                idxs = .!isnothing.(out)
                out = out[idxs]

                m = DimensionalData.metadata.(out)
                out = map(out) do o
                    dropdims(mean(o, dims = Chan), dims = Chan)
                end
                out = stack(SessionID(lookup(_Q, SessionID)[idxs]), out, dims = 3)
                return out
            end
        end

        begin # * Supplemental average comodulograms
            f = SixPanel()
            gs = subdivide(f, 3, 2)
            map(gs, structures, S) do g, structure, s
                ax = Axis(g[1, 1]; title = structure, xlabel = "Phase frequency (Hz)",
                          ylabel = "Amplitude frequency (Hz)")
                s = dropdims(mean(s, dims = SessionID); dims = SessionID)
                s = upsample(s, 5, 1)
                s = upsample(s, 5, 2)
                p = heatmap!(ax, s; colormap = seethrough(reverse(sunrise)), rasterize = 5)
                Colorbar(g[1, 2], p; label = "Modulation index")
            end
            addlabels!(f, labelformat)
            display(f)
            wsave(plotdir("nested_dynamics", "comodulograms_$stimulus.pdf"), f)
        end

        if stimulus == "spontaneous" # * Plot into main figure
            structure = "VISl"
            ax = Axis(mgs[3][1, 1]; title = "$structure comodulogram",
                      xlabel = "Phase frequency (Hz)",
                      ylabel = "Amplitude frequency (Hz)")
            s = S[lookup(_Q, Structure) .== structure] |> only
            display(s)
            s = dropdims(mean(s, dims = SessionID); dims = SessionID)
            s = upsample(s, 5, 1)
            s = upsample(s, 5, 2)
            p = heatmap!(ax, s .* 10^4; colormap = seethrough(reverse(sunrise)),
                         rasterize = 5)
            Colorbar(mgs[3][1, 2], p; label = "Mean PAC (×10⁴)")
        end
    end

    uni = load_uni(; stimulus, vars)

    unidepths = getindex.(uni, :unidepths)
    layerints = getindex.(uni, :layerints)
    layernames = getindex.(uni, :layernames)
    layernums = getindex.(uni, :layernums)

    begin # * Normalize amplitudes and generate a burst mask
        r = [abs.(uni[i][:r]) for i in eachindex(uni)]

        ϕ = [uni[i][:ϕ] for i in eachindex(uni)]
        ϕ = [mod2pi.(x .+ pi) .- pi for x in ϕ]

        uni = []
        PAC = progressmap(ϕ, r) do ϕ, r
            pac(ϕ, r; dims = Trial)
        end
        # ϕ_pref = map(ϕ, r) do ϕ, r
        #     peaks = map(eachslice(ϕ; dims = (𝑡, Depth)),
        #                 eachslice(r; dims = (𝑡, Depth))) do ϕ, r
        #         phipeak(r, ϕ; n = 20)
        #     end
        # end
        # begin
        #     fax = heatmap(decompose(PAC[2])...; axis = (; yreversed = true))
        #     f = fax.figure
        #     selection = (𝑡(0.3u"s"..0.45u"s"), Depth(0.15..0.4))
        #     tortinset!(f[1, 1], peaks[idx],
        #                colormap = seethrough(structurecolors[struc], alphamin, 1),
        #                halign = 0.3, valign = 0.6, color = structurecolors[struc])
        # end
        GC.gc()
    end

    begin # * Supplemental figure: spatiotemporal PAC over all regions
        f = SixPanel()
        gs = subdivide(f, 3, 2)
        for (g, l, P) in zip(gs, layerints, PAC)
            s = metadata(P)[:structure]
            ax = Axis(g[1, 1]; title = s, yreversed = true,
                      limits = (nothing, (extrema(lookup(P, Depth)))), xlabel = "Time (s)")
            p = plotlayermap!(ax, ustripall(P[𝑡 = SpatiotemporalMotifs.INTERVAL]), l;
                              rasterize = 5) |>
                first
            Colorbar(g[1, 2], p, label = "PAC")

            if s == "VISl"
                ax = Axis(mgs[6][1, 1]; title = "$s spatiotemporal PAC", yreversed = true,
                          limits = (nothing, (extrema(lookup(P, Depth)))),
                          xlabel = "Time (s)")
                p = plotlayermap!(ax,
                                  ustripall(P[𝑡 = SpatiotemporalMotifs.INTERVAL]) .* 10^3,
                                  l; rasterize = 5) |>
                    first
                Colorbar(mgs[6][1, 2], p, label = "PAC (×10³)")
            end
        end
        addlabels!(f, labelformat)
        wsave(plotdir("nested_dynamics", "supplemental_spatiotemporal_pac.pdf"), f)
        f
    end

    r = []
    ϕ = []
end

begin # * Layer-wise PAC
    begin
        lfile = datadir("layerwise_pac.jld2")
        Q = Q[stimulus = At(stimulus)]
        if isfile(lfile)
            pacc, peaks = load(lfile, "pacc", "peaks")
        else
            begin
                out = load_calculations(Q; stimulus, vars)
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
                        pac(p, a; dims = (𝑡, :changetime))
                    end
                end
                GC.gc()
                pacc = map(pacc) do pa
                    map(pa) do p
                        p = p[Depth = Near(unidepths)]
                        set(p, Depth => unidepths)
                    end
                end
                # pacc = stack.([SessionID(sessionids)], pacc; dims = 2)
                pacc = [cat(p...; dims = SessionID(sessionids)) for p in pacc]
                GC.gc()
            end
            begin # * Polar plot of coupling angle peaks (find peaks? what if there is more than 1? )
                function phipeak(r, ϕ; n = 50)
                    ϕ = mod2pi.(ϕ .+ pi) .- pi
                    ϕ = ModulationIndices.tortbin(ϕ; n)
                    h = [mean(r[ϕ .== i]) for i in 1:n]
                    angles = range(start = -π + π / n, stop = π - π / n, length = n) |>
                             collect
                    _, i = findmax(h)
                    phimax = angles[i]
                end
                function phipeaks(r, ϕ; kwargs...)
                    peaks = map(eachslice(r, dims = Depth),
                                eachslice(ϕ, dims = Depth)) do a, b
                        phipeak(a, b; kwargs...)
                    end
                    out = ϕ[1, :, 1] # Just a depth slice
                    out .= peaks
                    return out
                end
                peaks = progressmap(rs, ϕs) do r, ϕ
                    map(r, ϕ) do a, p
                        phipeaks(a, p; n = 20)
                    end
                end
                peaks = map(peaks) do pa
                    map(pa) do p
                        p = p[Depth = Near(unidepths)]
                        set(p, Depth => unidepths)
                    end
                end
                peaks = stack.([SessionID(sessionids)], peaks; dims = 2)
            end
            tagsave(lfile, @strdict pacc peaks)
        end
    end

    begin # * Plot layer-wise PAC
        ax = Axis(mgs[4]; xlabel = "Cortical depth (%)", ylabel = "Median PAC",
                  xtickformat = depthticks, limits = ((0, 1), nothing),
                  title = "Layerwise θ-γ PAC")
        for (i, p) in enumerate(pacc)
            s = structures[i]
            # p = p ./ sum(p, dims = 1)
            μ, (σl, σh) = bootstrapmedian(p |> ustripall; dims = 2)

            mu = upsample(μ, 5)
            l = upsample(σl, 5) |> parent
            h = upsample(σh, 5) |> parent
            band!(ax, lookup(mu, 1), l, h; color = (structurecolormap[s], 0.2),
                  label = s)
            lines!(ax, lookup(mu, 1), collect(mu), color = structurecolormap[s], label = s,
                   alpha = 0.8,
                   linewidth = 4)
            scatter!(ax, lookup(μ, 1), collect(μ), color = structurecolormap[s], label = s,
                     markersize = 10, alpha = 0.8)
        end
        l = axislegend(ax; merge = true, nbanks = 2, position = :lt, framevisible = true,
                       labelsize = 10)
        layerints = load(datadir("grand_unified_layers.jld2"), "layerints")
        plotlayerints!(ax, layerints; axis = :x, newticks = false, flipside = false)
        f
    end

    begin # * Phase annotations
        μs = map(pacc) do p
            μ, (σl, σh) = bootstrapmedian(p |> ustripall; dims = 2)
        end .|> first
        if true
            gg = mgs[4]
            idx = Depth(Near(0.25))
            alphamin = 0.2
            struc = 1
            tortinset!(gg, peaks[struc][idx],
                       colormap = seethrough(structurecolors[struc], alphamin, 1),
                       halign = 0.3, valign = 0.6, color = structurecolors[struc])
            scatter!(ax, [idx.val.val], [μs[struc][idx]], color = structurecolors[struc],
                     markersize = 15,
                     strokecolor = :white, strokewidth = 2)

            idx = Depth(Near(0.3))
            struc = 2
            tortinset!(gg, peaks[struc][idx],
                       colormap = seethrough(structurecolors[struc], alphamin, 1),
                       halign = 0.01, valign = 0.4, color = structurecolors[struc])
            scatter!(ax, [idx.val.val], [μs[struc][idx]], color = structurecolors[struc],
                     markersize = 15,
                     strokecolor = :white, strokewidth = 2)

            idx = Depth(Near(0.9))
            struc = 6
            tortinset!(gg, peaks[struc][idx],
                       colormap = seethrough(structurecolors[struc], alphamin, 1),
                       halign = 0.87, valign = 0.95, color = structurecolors[struc])
            scatter!(ax, [idx.val.val], [μs[struc][idx]], color = structurecolors[struc],
                     markersize = 15,
                     strokecolor = :white, strokewidth = 2)

            idx = Depth(Near(0.9))
            struc = 3
            tortinset!(gg, peaks[struc][idx],
                       colormap = seethrough(structurecolors[struc], alphamin, 1),
                       halign = 0.9, valign = 0.48, color = structurecolors[struc])
            scatter!(ax, [idx.val.val], [μs[struc][idx]], color = structurecolors[struc],
                     markersize = 15,
                     strokecolor = :white, strokewidth = 2)
        end
    end

    begin
        ax = PolarAxis(mgs[5]; theta_as_x = false, thetalimits = (-0.1pi, 1.2pi),
                       rticks = 0:0.25:1, rtickformat = depthticks,
                       title = "Layerwise PAC angle")
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
            lines!(ax, lookup(mu, 1), collect(mu), color = muc |> collect, label = s,
                   colormap = c,
                   linewidth = 7)
        end
    end
end

begin
    addlabels!(mf, labelformat)
    wsave(plotdir("nested_dynamics", "nested_dynamics.pdf"), mf)
end
