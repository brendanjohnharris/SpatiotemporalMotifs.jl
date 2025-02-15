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

begin # * Parameters
    stimulus = r"Natural_Images"
    sessionid = SpatiotemporalMotifs.DEFAULT_SESSION_ID
    trial = 14 # SpatiotemporalMotifs.DEFAULT_TRIAL_NUM
    # ΔT = SpatiotemporalMotifs.INTERVAL |> 𝑡
    ΔT = (-0.02u"s" .. 0.5u"s") |> 𝑡
    structure = "VISl"
    ΔD = Depth(0.15 .. 0.8)
    depth_colormap = SpatiotemporalMotifs.layercolormap

    config = Dict{String, Any}()
    @pack! config = stimulus, sessionid, trial, structure
end

plot_data, data_file = produce_or_load(copy(config), datadir("plots"); filename = savepath,
                                       prefix = "fig1A") do config # * Extract plot data from full calculations
    @unpack stimulus, sessionid, trial, structure = config
    Q = calcquality(datadir("power_spectra"))
    file = savepath((; sessionid, stimulus, structure), "jld2")
    file = datadir("calculations", file)
    out = jldopen(file, "r") do file
        θ = file["ϕ"][:, :, trial]
        γ = file["y"][:, :, trial]
        r = file["r"]
        datadepths = file["streamlinedepths"]
        spikes = file["spiketimes"]
        out = Dict{String, Any}()
        @pack! out = θ, γ, r, datadepths, spikes
        return out
    end
    unitdepths = load_unitdepths(Q[SessionID = (lookup(Q, SessionID) .== sessionid),
                                   Structure = At([structure]),
                                   stimulus = (lookup(Q, :stimulus) .==
                                               stimulus)])
    out["unitdepths"] = unitdepths
    return out
end

begin
    @unpack θ, γ, r, datadepths, spikes, unitdepths = plot_data
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
        # θ = θ # for visualization
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
        amps = [0.75, 0.75, 1.0, 1.0, 0]

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
                 rasterize = 6)

        cl = Makie.Contours.contours(decompose(peaks)..., [0.25])
        cl = Makie.Contours.levels(cl)[1]
        for l in Makie.Contours.lines(cl)
            local ux, uy = Makie.Contours.coordinates(l)
            xs = [findfirst(lookup(S, 1) .> _x) for _x in ux]
            ys = [findfirst(lookup(S, 2) .> _y) for _y in uy]
            ixs = .!isnothing.(xs) .& .!isnothing.(ys)
            xs, ys = xs[ixs], ys[ixs]
            xx = Point3f.(zip(ux, uy, x[xs, ys]))
            # yy = [pop!(xx)]
            # while !isempty(xx)
            #     idx = findmin([norm(yy[end] .- x) for x in xx]), 1
            #     push!(yy, popat!(xx, idx))
            # end
            lines!(ax, xx; linewidth = 4, alpha = 0.2, color = getindex.([cmat], xs, ys))#, color = :crimson)
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
                scatter!(ax, zip(sp, ys, zs) .|> Point3f; markersize = 3, color = :black,
                         alpha = 0.3)
            end
        end

        ax.yreversed = true
        f
    end
end

save(plotdir("schematic", "fig1A.pdf"), f;
     px_per_unit = 10)
f

if false # * Save a few representative stimulus images
    import AllenNeuropixelsBase as ANB
    using PythonCall
    sessionid = SpatiotemporalMotifs.DEFAULT_SESSION_ID
    session = ANB.Session(sessionid)
    df = session.pyObject.stimulus_templates |> ANB.py2df

    map(enumerate(eachrow(df))) do (i, d)
        _img = d.unwarped
        img = fill(Makie.RGBA(0, 0, 0, 0), size(_img))
        is = _img[.!isnan.(_img)] ./ 255
        img[.!isnan.(_img)] .= Makie.RGBA.(is, is, is, 1.0)
        save(plotdir("schematic", "natural_images", "natural_images_$i.png"), img)
    end

    _img = df[1, :unwarped]
    img = fill(Makie.RGBA(0, 0, 0, 0), size(_img))
    img[.!isnan.(_img)] .= Makie.RGBA.(0.5, 0.5, 0.5, 1.0)
    save(plotdir("schematic", "natural_images", "natural_images_0.png"), img)
end
f
