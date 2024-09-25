#! /bin/bash
#=
exec julia -t auto "${BASH_SOURCE[0]}" "$@"
=#
using DrWatson
@quickactivate "SpatiotemporalMotifs"
using SpatiotemporalMotifs
using Peaks
@preamble
set_theme!(foresight(:physics))
set_theme!(; fontsize = 15)
stimulus = r"Natural_Images"
sessionid = SpatiotemporalMotifs.DEFAULT_SESSION_ID
trial = 14 #10# SpatiotemporalMotifs.DEFAULT_TRIAL_NUM
ΔT = 𝑡(-0.25u"s" .. 0.75u"s")
δT = 𝑡(0u"s" .. 0.25u"s")
Q = calcquality(datadir("power_spectra"))

for structure in structures
    begin
        file = savepath((; sessionid, stimulus, structure), "jld2")
        file = datadir("calculations", file)
        file = jldopen(file, "r")

        layernames = ToolsArray(file["layerinfo"][1], (Depth(file["streamlinedepths"])))

        V = file["V"][ΔT][:, :, trial] .* 1000 # V to mV
        V = set(V, Depth(file["streamlinedepths"]))

        x = file["x"][ΔT][:, :, trial] .* 1000 # V to mV
        x = set(x, Depth(file["streamlinedepths"]))
        y = file["y"][ΔT][:, :, trial] .* 1000 # V to mV
        y = set(y, Depth(file["streamlinedepths"]))

        ϕ = file["ϕ"][ΔT][:, :, trial]
        ϕ = set(ϕ, Depth(file["streamlinedepths"]))
        r = file["r"][ΔT][:, :, trial] .* 1000 # V to mV
        r = set(r, Depth(file["streamlinedepths"]))
        _r = file["r"][ΔT] .* 1000 # V to mV
        r̂ = HalfZScore(_r, dims = [1, 3])(_r) # Normalized over time and trials
        _r = set(_r[:, :, trial], Depth(file["streamlinedepths"]))
        r̂ = set(r̂[:, :, trial], Depth(file["streamlinedepths"]))

        k = file["k"][ΔT][:, :, trial]
        k = set(k, Depth(file["streamlinedepths"]))
        k = uconvert.(u"mm^-1", k)

        ω = file["ω"][ΔT][:, :, trial]
        ω = set(ω, Depth(file["streamlinedepths"]))
        v = file["v"][ΔT][:, :, trial]
        v = set(v, Depth(file["streamlinedepths"]))

        spikes = file["spiketimes"]

        close(file)

        unitdepths = load_unitdepths(Q[SessionID = (lookup(Q, SessionID) .== sessionid),
                                       Structure = At(structures),
                                       stimulus = (lookup(Q, :stimulus) .==
                                                   string(stimulus))])

        # spikes = AN.getspiketimes(sessionid, structure)
        spikes = map(collect(spikes)) do (u, sp)
            t = ustripall(only(refdims(x, :changetime)))
            sp = sp[sp .∈ [t - 0.25 .. t + 0.75]] .- t
            u => sp
        end |> Dict
        spikes = Dict(k => v for (k, v) in spikes if 1 ≤ length(v))
        unitdepths = [unitdepths[unitdepths.id .== u, :streamlinedepth]
                      for u in keys(spikes)] .|>
                     only .|> Float64
        spikes = values(spikes) |> collect
        spikes = ToolsArray(spikes, (Depth(unitdepths),))
        # sort!(spikes)
    end

    begin # * Set up figure layout
        f = Figure(size = (720, 360) .* 1.25)

        # g2 = subdivide(f[1, 1:3], 1, 2) # * Middle row
        g3 = GridLayout(f[1:2, 1]) # * Bottom left
        g4 = subdivide(f[1:2, 2], 2, 1) # * Bottom center
        g5 = subdivide(f[1:2, 3], 2, 1) # * Bottom right

        rowsize!(g3, 1, Relative(0.4))
        # rowsize!(f.layout, 1, Relative(0.35))
        # rowsize!(f.layout, 2, Relative(0.35))
    end

    begin # * Heatmaps

        # * LFP
        vax = Axis(g3[1, 1]; title = "$structure LFP (mV)", yreversed = true,
                   xticks = [-0.25, 0.0, 0.25, 0.5], xlabel = "Time (s)",
                   limits = ((-0.25, 0.75), (0.01, 0.95)), yticksvisible = false)
        vp = plotlayermap!(vax, V, layernames; colormap = lfpcolormap,
                           colorrange = extrema(V)) |> first
        Colorbar(g3[1, 2], vp)

        # vax2 = Axis(gls[2, 1][1, 1], yreversed = true, title = "$structure LFP (mV)",
        #             xticks = [0.0, 0.125, 0.25])
        # vp2 = plotlayermap!(vax2, V[δT], layernames; colormap = :bone,
        #                     colorrange = extrema(V)) |> first
        # Colorbar(gls[2, 1][1, 2], vp2)

        # * Theta and Gamma
        # xax = Axis(gls[1, 2][1, 1], title = "θ LFP (mV)", yreversed = true,
        #            xticks = [-0.25, 0.0, 0.25, 0.5])
        # xp = plotlayermap!(xax, x, layernames; colormap = lfpcolormap,
        #                    colorrange = symextrema(x)) |> first
        # Colorbar(gls[1, 2][1, 2], xp)

        # yax = Axis(gls[2, 2][1, 1], title = "γ LFP (mV)", yreversed = true,
        #            xticks = [0.0, 0.125, 0.25])
        # yp = plotlayermap!(yax, y[δT], layernames; colormap = lfpcolormap,
        #                    colorrange = symextrema(y)) |> first
        # Colorbar(gls[2, 2][1, 2], yp)

        # * Phase and amplitude maps
        phax = Axis(g4[1][1, 1], title = "θ phase (radians)", yreversed = true,
                    xticks = [-0.25, 0.0, 0.25, 0.5], xlabel = "Time (s)",
                    limits = ((-0.25, 0.75), (0.01, 0.95)), yticksvisible = false)
        php = plotlayermap!(phax, ustripall(ϕ), layernames; colormap = phasecolormap,
                            colorrange = [-pi, pi], domain = -π .. π) |> first
        Colorbar(g4[1][1, 2], php; ticks = ([-pi, 0, pi], ["-𝜋", "0", "𝜋"]))

        rax = Axis(g4[2][1, 1], title = "γ amplitude (a.u.)", yreversed = true,
                   xticks = [0.0, 0.125, 0.25], xlabel = "Time (s)",
                   limits = ((0, 0.25), (0.01, 0.95)), yticksvisible = false)
        rp = plotlayermap!(rax, ustripall(_r[δT]) .* 100, layernames;
                           colormap = binarysunset,
                           colorrange = [0, maximum(ustripall(_r[δT]) .* 100)]) |> first
        Colorbar(g4[2][1, 2], rp)

        # * Wavenumber and Theta-gamma coupling
        kax = Axis(g5[1][1, 1], title = "θ wavenumber ($(unit(eltype(k))))",
                   yreversed = true, xticks = [-0.25, 0.0, 0.25, 0.5], xlabel = "Time (s)",
                   limits = ((-0.25, 0.75), (0.01, 0.95)), yticksvisible = false)
        colorrange = (-10, 10)
        K = ustripall(k)
        K[K .> colorrange[2] .+ 0.1] .= colorrange[2] .+ 0.1
        K[K .≤ colorrange[1] .- 0.1] .= colorrange[1] .- 0.1
        kp = plotlayermap!(kax, K, layernames; colormap = lfpcolormap,
                           colorrange, arrows = (150, 4), arrowsize = 10,
                           arrowcolor = (:white, 0.4), lengthscale = 0.1,
                           highclip = :cornflowerblue, lowclip = :crimson) |>
             first
        plotlayermap!(kax, ustripall(ω) .< 0; colormap = cgrad([:transparent, :white]))
        # plotlayermap!(kax, ustripall(ω) .< 0; colormap = cgrad([:transparent, :white]))
        Colorbar(g5[1][1, 2], kp)

        # hideydecorations!(xax)
        # hideydecorations!(vax2)
        # hideydecorations!(yax)
        # hideydecorations!(phax)
        # hideydecorations!(rax)
        # hideydecorations!(kax)

        vax.yticksvisible = phax.yticksvisible = kax.yticksvisible = false

        f
    end

    begin # * Plot masked peaks over phi
        pks, proms = findpeaks(r̂; minprom = 2.5)
        mask = maskpeaks(r̂; minprom = 2.5) .> 0

        ax = Axis(g5[2][1, 1], title = "θ-γ PAC", xticks = [0.0, 0.125, 0.25],
                  yreversed = true, xlabel = "Time (s)",
                  limits = ((0, 0.25), (0.01, 0.95)),
                  yticksvisible = false)

        pp = plotlayermap!(ax, ustripall(ϕ), layernames;
                           colormap = phasecolormap,
                           colorrange = (-pi, pi), domain = -π .. π) |> first
        # for i in eachindex(spikes)
        #     depth = lookup(spikes, Depth)[i]
        #     depth = fill(depth, length(spikes[i]))
        #     scatter!(ax, spikes[i], depth, color = :white,
        #              markersize = 5)
        # end
        contour!(ax, ustripall(mask[𝑡 = δT][Depth = 0.05 .. 0.95]); levels = [0.5],
                 color = (cucumber, 1),
                 linewidth = 2)
        pks = filter(!isempty, pks)
        pks = [p[𝑡 = δT] for p in pks]
        pks = filter(!isempty, pks)
        for d in lookup(pks, Depth)
            p = times(pks[Depth = At(d)])
            d = fill(d, length(p))
            scatter!(ax, collect(ustripall(p)), d, color = cucumber, markersize = 10,
                     strokewidth = 1, strokecolor = :white)
        end
        Colorbar(g5[2][1, 2], pp; ticks = ([-pi, 0, pi], ["-𝜋", "0", "𝜋"]))
        # ax.limits = (extrema(ustripall(δT.val)), extrema(lookup(r̂, 2)))
        # hideydecorations!(ax)
        f
    end

    if false # * Spiky? Spikiness is important. Gamma power is doubled for visualization
        ax = Axis(g2[2][1, 1], xlabel = "Time (s)",
                  title = "Nested oscillations")
        d = 7 # Which channel?
        lines!(ax, lookup(V[ΔT][:, d], 𝑡) |> ustripall, V[ΔT][:, d] |> ustripall |> collect,
               color = (cucumber, 0.8),
               label = "LFP",
               linewidth = 3)
        lines!(ax, decompose(x[ΔT][:, d] .- minimum(x[ΔT][:, d]) .- 0.5 |> ustripall)...,
               color = (crimson, 0.8),
               label = L"\theta", linewidth = 3)
        # lines!(ax, bandpass(V |> ustripall, 15 .. 30)[ustripall(T), d ] .- 0.1) # harmonic
        lines!(ax,
               decompose(r[ΔT][:, d] .* 2 .- minimum(V[ΔT][:, d]) .+ minimum(x[ΔT][:, d]) .-
                         0.8 |>
                         ustripall)...,
               color = (cornflowerblue, 0.8), label = L"\gamma", linewidth = 3)

        ax2 = Axis(g2[2][1, 1]; limits = ((0, 0.25), (0.01, 0.95)),
                   yreversed = true)
        for i in eachindex(spikes)
            depth = lookup(spikes, Depth)[i]
            depth = fill(depth, length(spikes[i]))
            scatter!(ax2, spikes[i], depth, color = (:black, 0.4),
                     markersize = 5)
        end
        axislegend(ax, orientation = :horizontal, position = :lt, framevisible = true,
                   padding = fill(4, 4))
        vlines!(ax, [0, 0.25], color = (:black, 0.5), linestyle = :dash, linewidth = 3)
        hideydecorations!(ax)
        hidedecorations!(ax2)
        ax.limits = ((-0.25, 0.75), (-0.6, 0.3))
        f
    end

    begin # * Plot power spectrum during spontaneous
        ff = Figure()
        D = @strdict sessionid structure
        push!(D, "stimulus" => "spontaneous")
        sfilename = savepath(D, "jld2", datadir("power_spectra"))
        S = load(sfilename, "S")[:, 4:end]
        depths = lookup(S, Depth)

        s = mean(S, dims = Depth)[:, 1]
        pks, vals = findmaxima(s, 10)
        pks, proms = peakproms(pks, s)
        promidxs = (proms ./ vals .> 0.25) |> collect
        maxs = maximum(S, dims = Depth)[:, 1]
        pks = pks[promidxs]
        pks = TimeseriesTools.freqs(s)[pks]
        vals = s[Freq(At(pks))]

        colorrange = extrema(depths)
        ax = Axis(ff[1, 1]; xscale = log10, yscale = log10, xtickformat = "{:.0f}",
                  limits = ((2, 100), (1e-12, 5e-10)), xgridvisible = true,
                  ygridvisible = true, topspinevisible = true,
                  xminorticksvisible = true, yminorticksvisible = true,
                  xminorgridvisible = true, yminorgridvisible = true,
                  xminorgridstyle = :dash,
                  title = "VISp LFP spectrum")
        p = traces!(ax, S[2:end, :]; colormap = cgrad(sunset, alpha = 0.4),
                    linewidth = 3, colorrange)
        scatter!(ax, ustrip.(pks), collect(ustrip.(vals)), color = :black,
                 markersize = 10, marker = :dtriangle)
        text!(ax, ustrip.(pks), collect(ustrip.(vals));
              text = string.(round.(eltype(pks), pks, digits = 1)),
              align = (:center, :bottom), color = :black, rotation = 0,
              fontsize = 12,
              offset = (0, 3))

        c = Colorbar(ff[1, 2]; label = "Channel depth (μm)", colorrange,
                     colormap = sunset)
        f
    end

    begin # * Save figure
        # for i in 1:(length(gls) - 1)
        #     try
        #         colgap!(gls[i], 1, Relative(0.02))
        #     catch e
        #     end
        # end
        # rowgap!(f.layout, 1, Relative(0.1))
        # rowsize!(f.layout, 1, Relative(0.3))
        # colgap!(f.layout, 1, Relative(0.1))
        addlabels!(f)

        save(plotdir("schematic", "single_trial_schematic_$structure.pdf"), f)
        f
    end

    if false # * Spike MUA spectrum
        S = AN.Session(sessionid)
        spikes = AN.getspikes(S)
        ts = 3655:0.005:4567
        B = x -> HistBins(x; bins = ts)

        bs = progressmap(collect(values(spikes))[1:100]; parallel = true) do x
            x = x[minimum(ts) .< x .< maximum(ts)]
            b = B(x)(x)
            b = length.(b)
        end
        b = set(bs[1], sum(bs))
        b = set(b, Dim{:bin} => 𝑡(ts[1:(end - 1)]))

        s = log10.(spectrum(b, 0.5)[𝑓 = 1 .. 50])
        lines(freqs(s), s)
        ax = current_axis()
        ax.xlabel = "Hz"
        ax.ylabel = "power"
        ax.title = "MUA spectrum"
        current_figure()
    end
    f
end
