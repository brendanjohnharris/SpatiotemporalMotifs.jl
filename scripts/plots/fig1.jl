#! /bin/bash
#=
exec julia +1.10.10 -t auto "${BASH_SOURCE[0]}" "$@"
=#
using DrWatson
@quickactivate "SpatiotemporalMotifs"
using SpatiotemporalMotifs
import AllenNeuropixelsBase as ANB
using PythonCall
using Peaks
@preamble
set_theme!(foresight(:physics))
set_theme!(; fontsize = 15)

begin # * Parameters
    stimulus = r"Natural_Images"
    sessionid = SpatiotemporalMotifs.DEFAULT_SESSION_ID
    trial = SpatiotemporalMotifs.DEFAULT_TRIAL_NUM

    config = @strdict stimulus sessionid trial
end

plot_data, data_file = produce_or_load(copy(config), calcdir("plots");
                                       filename = savepath("fig1")) do config
    @unpack stimulus, sessionid, trial = config
    Q = calcquality(calcdir("power_spectra"))

    plot_data = map(structures) do structure
        file = savepath((; sessionid, stimulus, structure), "jld2")
        file = calcdir("calculations", file)
        out = jldopen(file, "r") do file
            layernames = ToolsArray(file["layerinfo"][1],
                                    (Depth(file["streamlinedepths"]),))
            datadepths = file["streamlinedepths"]
            spikes = file["spiketimes"]

            V = file["V"][:, :, trial] .|> Float32
            x = file["θ"][:, :, trial] .|> Float32
            y = file["γ"][:, :, trial] .|> Float32
            ϕ = file["ϕ"][:, :, trial] .|> Float32
            k = file["k"][:, :, trial] .|> Float32
            ω = file["ω"][:, :, trial] .|> Float32
            v = file["v"][:, :, trial] .|> Float32
            θ = deepcopy(ϕ)
            γ = deepcopy(y)
            r = file["r"] .|> Float32 # Don't select trial, because we need to normalize over time later on

            out = @strdict V x y ϕ k ω v θ γ r layernames datadepths spikes
            return out
        end
        unitdepths = load_unitdepths(Q[SessionID = (lookup(Q, SessionID) .== sessionid),
                                       Structure = At([structure]),
                                       stimulus = (lookup(Q, :stimulus) .==
                                                   stimulus)])
        out["unitdepths"] = unitdepths

        begin # * Stimulus examples
            sessionid = SpatiotemporalMotifs.DEFAULT_SESSION_ID
            session = ANB.Session(sessionid)
            df = session.pyObject.stimulus_templates |> ANB.py2df

            imgs = map(enumerate(eachrow(df))) do (i, d)
                _img = d.unwarped
                img = fill(Makie.RGBA(0, 0, 0, 0), size(_img))
                is = _img[.!isnan.(_img)] ./ 255
                img[.!isnan.(_img)] .= Makie.RGBA.(is, is, is, 1.0)
                return img
            end

            _img = df[1, :unwarped]
            img = fill(Makie.RGBA(0, 0, 0, 0), size(_img))
            img[.!isnan.(_img)] .= Makie.RGBA.(0.5, 0.5, 0.5, 1.0)
            prepend!(imgs, [img])

            out["stimulus_examples"] = imgs |> unique
        end
        return out
    end
    return Dict(structures .=> plot_data)
end

begin # ? Figure 1A
    @info "Plotting wave schematic"
    structure = "VISl"
    ΔT = (-0.0u"s" .. 0.52u"s") |> 𝑡
    ΔD = Depth(0.15 .. 0.8)
    depth_colormap = SpatiotemporalMotifs.layercolormap
    begin
        data = plot_data[structure]
        @unpack θ, γ, r, datadepths, spikes, unitdepths = data
        begin
            # Normalized theta by taking phase
            θ = θ[ΔT]
            θ = set(θ, Depth(datadepths))[ΔD] |> ustripall
            θ = -1im * θ .|> exp .|> real
            θ[1, :] .= θ[2, :]

            γ = γ[ΔT]
            γ = set(γ, Depth(datadepths))[ΔD] |> ustripall

            _r = deepcopy(r)[ΔT]
            r = r[ΔT]
            r = set(r, Depth(datadepths))[ΔD] |> ustripall
            r̂ = HalfZScore(_r, dims = [1, 3])(_r) # Normalized over time and trials
            _r = set(_r[:, :, trial], Depth(datadepths))[ΔD]
            r = set(r̂[:, :, trial], Depth(datadepths))[ΔD] |> ustripall

            spikes = map(collect(spikes)) do (u, sp)
                t = ustripall(only(refdims(θ, :changetime)))
                sp = sp[sp .∈ [t - 0.25 .. t + 0.75]] .- t
                u => sp
            end |> Dict
            spikes = Dict(k => v for (k, v) in spikes if 1 ≤ length(v))
            unitdepths = [unitdepths[unitdepths.id .== u, :streamlinedepth]
                          for u in keys(spikes)] .|>
                         only .|> Float64
            spikes = values(spikes) |> collect
            spikes = ToolsArray(spikes, (Depth(unitdepths),))
            # * Jitter spike depths for visualization
            spikes = set(spikes,
                         Depth => lookup(spikes, Depth) .+ 0.01 .* randn(length(spikes)))
            sidxs = sortperm(lookup(spikes, Depth))
            spikes = set(spikes[sidxs],
                         Depth => DimensionalData.Lookups.ForwardOrdered)
            spikes = spikes[ΔD]
            spikes = [s[s .∈ [ustripall((ΔT.val.left + 0.02u"s") .. ΔT.val.right)]]
                      for s in spikes]
            # sort!(spikes)
        end

        begin # * Normalization and gamma patches
            θ = .-θ # for visualization
            γ = set(γ, ZScore(γ |> parent; dims = 1)(γ))
            r = set(r, MinMax(r |> parent; dims = 1)(r))

            x = upsample(θ, 8, 2) # Smoother signal
        end

        begin
            # * Generate peaks, weighted by depths
            f = Figure(; size = (800, 500), backgroundcolor = :transparent)

            peaks = similar(x)
            peakgrid = Iterators.product(lookup(x)...)
            g = (x, y; σ) -> exp.(.-(x .^ 2 + y .^ 2 / σ^2) ./ 0.00025)

            t1 = 0.39
            t2 = 0.14
            dt = 0.05
            is = [(t1, 0.7) # * Lower 2
                  (t1 + dt, 0.25) # * Upper 2
                  (t2, 0.65) # * Lower 1
                  (t2 + dt, 0.3) # * Upper 1
                  (0.1, 0.54)]
            σs = [2, 4, 6, 7, 6]
            amps = [0.75, 0.75, 1.0, 1.0, 0] .* 1.5

            peaks .= sum([[a .* g(x .- i[1], y .- i[2]; σ) for (x, y) in peakgrid]
                          for (i, σ, a) in zip(is, σs, amps)])

            G = ((xy) -> sin(500xy[1] + 10xy[2])).(peakgrid)
            G = G .* peaks .^ 2
            S = x .+ G

            ax = Axis3(f[1, 1]; backgroundcolor = :transparent, viewmode = :stretch,
                       perspectiveness = 0.0)
            ax.xypanelvisible = ax.yzpanelvisible = ax.xzpanelvisible = false
            hidedecorations!(ax)
            ax.elevation = 0.9
            ax.azimuth = -π / 2 + 0.1
            # * Plot the theta wave
            L = tukey(size(x, :𝑡), 0.4)
            L[(length(L) ÷ 2):end] .= 1 # Opaque right edge
            l = repeat(L, 1, size(x, 2))
            l = MinMax(l)(l)
            cmat = set(x, repeat(lookup(x, Depth)', size(x, 1), 1))
            cmat = MinMax(cmat)(cmat)
            cmat = [cgrad(depth_colormap)[i] for i in cmat]
            cmat = [Makie.RGBA(c.r, c.g, c.b, _l) for (c, _l) in zip(cmat, l)]
            surface!(ax, decompose(S)..., color = collect(cmat), alpha = 0.9,
                     specular = 0.1,
                     rasterize = 5)

            cl = Makie.Contours.contours(decompose(peaks)..., [0.25])
            cl = Makie.Contours.levels(cl)[1]
            for l in Makie.Contours.lines(cl)
                local ux, uy = Makie.Contours.coordinates(l)
                xs = [findfirst(lookup(S, 1) .> _x) for _x in ux]
                ys = [findfirst(lookup(S, 2) .> _y) for _y in uy]
                ixs = .!isnothing.(xs) .& .!isnothing.(ys)
                xs, ys = xs[ixs], ys[ixs]
                xx = Point3f.(zip(ux, uy, x[xs, ys]))
                lines!(ax, xx; linewidth = 4, alpha = 0.2,
                       color = getindex.([cmat], xs, ys))#, color = :crimson)
            end
            thelines = range(start = 1, stop = size(S, 2), length = 8) .|> round .|> Int
            for i in axes(S, 2)[thelines]
                sy = lookup(S, 2)[i]
                z = S[:, i]
                xs = lookup(z, 1)
                color = cgrad(depth_colormap)[i ./ (size(S, 2))]
                color = fill(color, length(xs))
                color = [Makie.RGBA(c.r, c.g, c.b, l) for (c, l) in zip(color, L)]
                lines!(ax, xs, fill(sy, length(xs)), z |> collect;
                       color,
                       linewidth = 4, linecap = :round, alpha = 0.7)
            end

            begin # * Add spike dots
                spikes = set(spikes, reverse(collect(spikes))) # Looks better
                for (i, sp) in enumerate(spikes)
                    ys = fill(lookup(spikes, Depth)[i], length(sp))
                    zs = [S[𝑡(Near(x)), Depth(Near(y))] for (x, y) in zip(sp, ys)]
                    scatter!(ax, zip(sp, ys, zs) .|> Point3f; markersize = 3,
                             color = :black,
                             alpha = 0.3)
                end
            end

            ax.yreversed = true
            hidespines!(ax)
            f
        end
    end

    wsave(plotdir("fig1", "fig1A.pdf"), f;
          px_per_unit = 10)
    f

    begin # * Save fig1A csv files
        using CSV, DataFrames

        outdir = datadir("source_data", "fig1A")
        mkpath(outdir)

        # 1. Surface data (S) - the main 3D surface with theta + gamma overlay
        surface_df = DataFrame(S)
        rename_cols!(surface_df)
        CSV.write(joinpath(outdir, "surface.csv"), surface_df)

        # 2. Base theta signal (x) before gamma overlay
        theta_df = DataFrame(x)
        rename_cols!(theta_df)
        CSV.write(joinpath(outdir, "theta.csv"), theta_df)

        # 3. Gamma peaks envelope (peaks)
        peaks_df = DataFrame(peaks)
        rename_cols!(peaks_df)
        CSV.write(joinpath(outdir, "gamma_peaks.csv"), peaks_df)

        # 4. Gamma oscillation (G) - G is a plain array, need to reconstruct with dims
        G_da = set(x, G)  # Use x's dimensions for G
        gamma_df = DataFrame(G_da)
        rename_cols!(gamma_df)
        CSV.write(joinpath(outdir, "gamma_oscillation.csv"), gamma_df)

        # 5. Spike data
        spike_df = DataFrame("time (s)" => Float64[], "cortical depth (%)" => Float64[])
        for (i, sp) in enumerate(spikes)
            d = lookup(spikes, Depth)[i]
            for t in sp
                push!(spike_df, (t, d))
            end
        end
        CSV.write(joinpath(outdir, "spikes.csv"), spike_df)

        @info "Saved Figure 1A data to $outdir"
    end

    begin # * Save a few representative stimulus images
        map(enumerate(plot_data["VISp"]["stimulus_examples"])) do (i, img)
            wsave(plotdir("fig1", "stimulus_examples", "natural_images_$i.png"), img)
        end
    end
    f
end

begin # ? Figure 1: C--G
    @info "Plotting single-trial schematic"
    ΔT = 𝑡(-0.25u"s" .. 0.75u"s")
    δT = 𝑡(0u"s" .. 0.25u"s")
    for structure in structures
        data = plot_data[structure]
        begin
            @unpack V, x, y, ϕ, k, ω, v, θ, γ, r, layernames, datadepths, spikes, unitdepths = data

            V = V[ΔT] .* 1000 # V to mV
            V = set(V, Depth(datadepths))

            x = x[ΔT] .* 1000 # V to mV
            x = set(x, Depth(datadepths))
            y = y[ΔT] .* 1000 # V to mV
            y = set(y, Depth(datadepths))

            ϕ = ϕ[ΔT]
            ϕ = set(ϕ, Depth(datadepths))

            _r = deepcopy(r[ΔT]) .* 1000 # V to mV
            r = r[ΔT] .* 1000 # V to mV
            r = set(r, Depth(datadepths))
            r̂ = HalfZScore(_r, dims = [1, 3])(_r) # Normalized over time and trials
            _r = set(_r[:, :, trial], Depth(datadepths))
            r̂ = set(r̂[:, :, trial], Depth(datadepths))

            k = k[ΔT]
            k = set(k, Depth(datadepths))
            k = uconvert.(u"mm^-1", k)

            ω = ω[ΔT]
            ω = set(ω, Depth(datadepths))
            v = v[ΔT]
            v = set(v, Depth(datadepths))

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

            g3 = GridLayout(f[1:2, 1]) # * Bottom left
            g4 = subdivide(f[1:2, 2], 2, 1) # * Bottom center
            g5 = subdivide(f[1:2, 3], 2, 1) # * Bottom right

            rowsize!(g3, 1, Relative(0.4))
        end

        begin # * Heatmaps

            # * LFP
            vax = Axis(g3[1, 1]; title = "$structure LFP (mV)", yreversed = true,
                       xticks = [-0.25, 0.0, 0.25, 0.5], xlabel = "Time (s)",
                       limits = ((-0.25, 0.75), (0.01, 0.95)), yticksvisible = false)
            vp = plotlayermap!(vax, V, layernames; colormap = lfpcolormap,
                               colorrange = extrema(V)) |> first
            Colorbar(g3[1, 2], vp)

            # * Phase and amplitude maps
            phax = Axis(g4[1][1, 1], title = "θ phase (radians)", yreversed = true,
                        xticks = [-0.25, 0.0, 0.25, 0.5], xlabel = "Time (s)",
                        limits = ((-0.25, 0.75), (0.01, 0.95)), yticksvisible = false)
            php = plotlayermap!(phax, ustripall(ϕ), layernames; colormap = phasecolormap,
                                colorrange = [-pi, pi], domain = -π .. π) |> first
            Colorbar(g4[1][1, 2], php; ticks = ([-pi, 0, pi], ["-𝜋", "0", "𝜋"]))

            rax = Axis(g4[2][1, 1], title = "γ amplitude (arb. units)", yreversed = true,
                       xticks = [0.0, 0.125, 0.25], xlabel = "Time (s)",
                       limits = ((0, 0.25), (0.01, 0.95)), yticksvisible = false)
            rp = plotlayermap!(rax, ustripall(_r[δT]) .* 100, layernames;
                               colormap = binarysunset,
                               colorrange = [0, maximum(ustripall(_r[δT]) .* 100)]) |> first
            Colorbar(g4[2][1, 2], rp)

            # * Wavenumber and Theta-gamma coupling
            kax = Axis(g5[1][1, 1], title = "θ wavenumber ($(unit(eltype(k))))",
                       yreversed = true, xticks = [-0.25, 0.0, 0.25, 0.5],
                       xlabel = "Time (s)",
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
            Colorbar(g5[1][1, 2], kp)

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
            f
        end

        if false # * Spiky? Spikiness is important. Gamma power is doubled for visualization
            ax = Axis(g2[2][1, 1], xlabel = "Time (s)",
                      title = "Nested dynamics")
            d = 7 # Which channel?
            lines!(ax, lookup(V[ΔT][:, d], 𝑡) |> ustripall,
                   V[ΔT][:, d] |> ustripall |> collect,
                   color = (cucumber, 0.8),
                   label = "LFP",
                   linewidth = 3)
            lines!(ax,
                   decompose(x[ΔT][:, d] .- minimum(x[ΔT][:, d]) .- 0.5 |> ustripall)...,
                   color = (crimson, 0.8),
                   label = L"\theta", linewidth = 3)
            lines!(ax,
                   decompose(r[ΔT][:, d] .* 2 .- minimum(V[ΔT][:, d]) .+
                             minimum(x[ΔT][:, d]) .-
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

        begin # * Save figure
            addlabels!(f, labelformat)

            wsave(plotdir("fig1", "single_trial_schematic_$structure.pdf"), f)
            f
        end

        if structure == "VISl" # * Save csv source data
            outdir = datadir("source_data", "fig1CtoG")
            mkpath(outdir)

            # Panel C: LFP heatmap
            lfp_df = DataFrame(ustripall(V))
            rename_cols!(lfp_df)
            rename!(lfp_df, names(lfp_df)[end] => "lfp (mV)")
            CSV.write(joinpath(outdir, "lfp.csv"), lfp_df)

            # Panel D: Theta phase heatmap
            phase_df = DataFrame(ustripall(ϕ))
            rename_cols!(phase_df)
            rename!(phase_df, names(phase_df)[end] => "Phase (radians)")
            CSV.write(joinpath(outdir, "theta_phase.csv"), phase_df)

            # Panel E: Gamma amplitude heatmap
            gamma_amp_df = DataFrame(ustripall(_r[δT]) .* 100)
            rename_cols!(gamma_amp_df)
            rename!(gamma_amp_df, names(gamma_amp_df)[end] => "Amplitude (arb. units)")
            CSV.write(joinpath(outdir, "gamma_amplitude.csv"), gamma_amp_df)

            # Panel F: Wavenumber heatmap
            wavenumber_df = DataFrame(ustripall(k))
            rename_cols!(wavenumber_df)
            rename!(wavenumber_df, names(wavenumber_df)[end] => "Wavenumber (mm^-1)")
            CSV.write(joinpath(outdir, "wavenumber.csv"), wavenumber_df)

            # Panel G: PAC mask and peaks
            pac_mask_df = DataFrame(mask[𝑡 = δT])
            rename_cols!(pac_mask_df)
            rename!(pac_mask_df, names(pac_mask_df)[end] => "Mask")
            CSV.write(joinpath(outdir, "pac_mask.csv"), pac_mask_df)

            # PAC peaks (gamma burst times and depths)
            pac_peaks_df = DataFrame("time (s)" => Float64[],
                                     "cortical depth (%)" => Float64[])
            for d in lookup(pks, Depth)
                p = times(pks[Depth = At(d)])
                for t in ustripall(p)
                    push!(pac_peaks_df, (t, d))
                end
            end
            CSV.write(joinpath(outdir, "pac_peaks.csv"), pac_peaks_df)

            @info "Saved Figure 1 C-G data to $outdir"
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
        if false # * Plot power spectrum during spontaneous
            ff = Figure()
            D = @strdict sessionid structure
            push!(D, "stimulus" => "spontaneous")
            sfilename = savepath(D, "jld2", calcdir("power_spectra"))
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
            ff
        end
        f
    end
end
