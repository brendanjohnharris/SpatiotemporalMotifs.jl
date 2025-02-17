#! /bin/bash
#=
exec julia -t auto "${BASH_SOURCE[0]}" "$@"
=#
using DrWatson
@quickactivate "SpatiotemporalMotifs"

import DimensionalData: metadata
using SpatiotemporalMotifs
@preamble
set_theme!(foresight(:physics))
Random.seed!(42)

plot_data, data_file = produce_or_load(Dict(), datadir("plots");
                                       filename = savepath,
                                       prefix = "fig6") do _
    stimulus = r"Natural_Images"
    vars = [:ϕ, :r]

    session_table = load(datadir("posthoc_session_table.jld2"), "session_table")
    oursessions = session_table.ecephys_session_id

    path = datadir("calculations")
    Q = calcquality(path)[Structure = At(structures)][SessionID = At(oursessions)]
    out = load_calculations(Q; stimulus, vars)
    Qs = calcquality(datadir("power_spectra"))[Structure = At(structures)][SessionID = At(oursessions)]
    unitdepths = load_unitdepths(Qs)

    begin # * Format spikes to depth dataframe. These have already been rectified
        spikes = map(out) do o
            regionalspikes = map(o) do _o
                DataFrame(collect.([keys(_o[:spiketimes]), values(_o[:spiketimes])]),
                          [:ecephys_unit_id, :spiketimes])
            end
            vcat(regionalspikes...)
        end
        D = vcat(spikes...)
        spikes = innerjoin(D, unitdepths, on = :ecephys_unit_id)
        filter!(:streamlinedepth => ∈(0 .. 1), spikes) # Cortex only
    end

    begin # * Calculate spike--phase and spike--amplitude coupling across layers. Takes about 30 minutes over 64 cores, 125 GB
        # The idea here is that we have this larger spike_lfp file, and subset it to the
        # relevant data for this plot
        if isfile(datadir("spike_lfp.jld2"))
            pspikes = load(datadir("spike_lfp.jld2"), "pspikes") # Delete this file to recalculate
        else
            pspikes = deepcopy(spikes)
        end
        requiredcols = [:pairwise_phase_consistency, :trial_pairwise_phase_consistency,
            :pairwise_phase_consistency_pvalue, :trial_pairwise_phase_consistency_pvalue,
            :pairwise_phase_consistency_angle, :pairwise_phase_consistency_angle,
            :spike_amplitude_coupling, :trial_spike_amplitude_coupling]
        if !all(hasproperty.([pspikes], requiredcols))
            for s in eachindex(structures)
                ϕ = getindex.(out[s], :ϕ)
                r = getindex.(out[s], :r)
                spc!(pspikes, ustripall.(ϕ)) # * PPC spike--phase coupling
                sac!(pspikes, ustripall.(r)) # * Mean normalized amplitude spike--amplitude coupling
            end
            save(datadir("spike_lfp.jld2"), "pspikes", pspikes)
        end
    end
    pspikes = subset(pspikes, :pairwise_phase_consistency => ByRow(!isnan))
    pspikes = subset(pspikes, :pairwise_phase_consistency => ByRow(>(0)))
    pspikes = subset(pspikes, :spike_amplitude_coupling => ByRow(!isnan))
    select!(pspikes, Not([:spiketimes])) # We don't need the heavy spiketimes for this plot
    select!(pspikes, Not([:filtering]))
    select!(pspikes, Not([:trial_pairwise_phase_consistency]))
    select!(pspikes, Not([:trial_pairwise_phase_consistency_pvalue]))
    select!(pspikes, Not([:trial_spike_amplitude_coupling]))

    unitdepths = subset(unitdepths, :stimulus => ByRow(==("spontaneous")))
    unitdepths = subset(unitdepths, :spc => ByRow(!isnan))
    unitdepths = subset(unitdepths, :spc => ByRow(>(0)))
    unitdepths = subset(unitdepths, :sac => ByRow(!isnan))

    return (@strdict pspikes unitdepths)
end

begin
    @unpack pspikes, unitdepths = plot_data
    begin # * Set up figure
        f = SixPanel()
        gs = subdivide(f, 3, 2)
        layerints = load(datadir("plots", "grand_unified_layers.jld2"), "layerints")
        bins = range(0, 1, length = 11)
    end

    begin # * Plot, for each structure, spontaneous SPC
        @info "Plotting spontaneous spike--phase coupling"
        ax = Axis(gs[1], xlabel = "Cortical depth (%)", ylabel = "PPC",
                  title = "Spike-phase coupling (θ)",
                  limits = ((0.05, 0.95), (0, nothing)),
                  xtickformat = depthticks)
        for structure in reverse(structures)
            idxs = unitdepths.structure_acronym .== structure
            allsesh_pspikes = @views unitdepths[idxs, :]
            mss = map(unique(allsesh_pspikes.ecephys_session_id)) do sesh
                _pspikes = @views allsesh_pspikes[allsesh_pspikes.ecephys_session_id .== sesh,
                                                  :]
                ys = _pspikes.spc
                xs = _pspikes.streamlinedepth .|> Float32
                ys = ys .|> Float32
                # hexbin(xs, ys)
                B = HistBins(xs; bins)
                _ms = rectify(B(ys), dims = :bin)
            end

            X = cat(mss...; dims = 2)
            X = map(eachrow(X)) do x
                x = vcat(x...)
            end
            _ms = pmap(bootstrapmean, X)
            # _ms = bootstrapmedian.(X)
            σs = last.(_ms)
            _ms = first.(_ms)
            _σl = first.(σs)
            _σh = last.(σs)
            # _ms, (_σl, _σh) = bootstrapmedian(X; dims = 2)
            idxs = .!isnan.(_ms) .& .!isnan.(_σl) .& .!isnan.(_σh)
            _ms = _ms[idxs]
            _σl = _σl[idxs]
            _σh = _σh[idxs]
            # _σs = map((args...) -> nansafe(confidence)(collect(args)) |> first, mss...)
            # _ms = map((args...) -> nansafe(mean)(collect(args)) |> first, mss...)

            _ls = lookup(_ms, 1)
            ms = upsample(_ms, 10)
            σl = upsample(_σl, 10)
            σh = upsample(_σh, 10)
            ls = lookup(ms, 1)
            band!(ax, ls, collect(σl), collect(σh);
                  color = (structurecolormap[structure], 0.3), label = structure)
            lines!(ax, ls, collect(ms), color = (structurecolormap[structure], 0.7),
                   label = structure)
            scatter!(ax, _ls, collect(_ms); color = structurecolormap[structure],
                     label = structure)
        end
        l = axislegend(ax, merge = true, nbanks = 2)
        plotlayerints!(ax, layerints; axis = :x, flipside = false, newticks = false,
                       bgcolor = Makie.RGBA(0, 0, 0, 0))
        reverselegend!(l)
    end

    begin # * Plot, for each structure, spontaneous SAC
        @info "Plotting spontaneous spike--amplitude coupling"
        ax = Axis(gs[3], xlabel = "Cortical depth (%)", ylabel = "SAC",
                  title = "Spike-amplitude coupling (γ)",
                  limits = ((0.05, 0.95), (1.1, 1.5)),
                  xtickformat = depthticks)
        for structure in reverse(structures)
            idxs = unitdepths.structure_acronym .== structure
            allsesh_pspikes = @views unitdepths[idxs, :]
            mss = map(unique(allsesh_pspikes.ecephys_session_id)) do sesh
                _pspikes = @views allsesh_pspikes[allsesh_pspikes.ecephys_session_id .== sesh,
                                                  :]
                ys = _pspikes.sac
                xs = _pspikes.streamlinedepth .|> Float32
                ys = ys .|> Float32
                # hexbin(xs, ys)
                B = HistBins(xs; bins)
                _ms = rectify(B(ys), dims = :bin)
            end

            X = cat(mss...; dims = 2)
            X = map(eachrow(X)) do x
                x = vcat(x...)
            end
            _ms = pmap(bootstrapmean, X)
            # _ms = bootstrapmedian.(X)
            σs = last.(_ms)
            _ms = first.(_ms)
            _σl = first.(σs)
            _σh = last.(σs)
            # _ms, (_σl, _σh) = bootstrapmedian(X; dims = 2)
            idxs = .!isnan.(_ms) .& .!isnan.(_σl) .& .!isnan.(_σh)
            _ms = _ms[idxs]
            _σl = _σl[idxs]
            _σh = _σh[idxs]
            # _σs = map((args...) -> nansafe(confidence)(collect(args)) |> first, mss...)
            # _ms = map((args...) -> nansafe(mean)(collect(args)) |> first, mss...)

            _ls = lookup(_ms, 1)
            ms = upsample(_ms, 10)
            σl = upsample(_σl, 10)
            σh = upsample(_σh, 10)
            ls = lookup(ms, 1)
            band!(ax, ls, collect(σl), collect(σh);
                  color = (structurecolormap[structure], 0.3), label = structure)
            lines!(ax, ls, collect(ms), color = (structurecolormap[structure], 0.7),
                   label = structure)
            scatter!(ax, _ls, collect(_ms); color = structurecolormap[structure],
                     label = structure)
        end
        l = axislegend(ax, merge = true, nbanks = 2)
        plotlayerints!(ax, layerints; axis = :x, flipside = false, newticks = false,
                       bgcolor = Makie.RGBA(0, 0, 0, 0))
        reverselegend!(l)
    end

    begin # * Plot, for each structure, SPC
        @info "Plotting structure-wise spike--phase coupling"
        ax = Axis(gs[2], xlabel = "Cortical depth (%)", ylabel = "PPC",
                  title = "Spike-phase coupling (θ)",
                  limits = ((0.05, 0.95), (0, nothing)),
                  xtickformat = depthticks)
        for structure in reverse(structures)
            idxs = pspikes.structure_acronym .== structure
            allsesh_pspikes = @views pspikes[idxs, :]
            mss = map(unique(allsesh_pspikes.ecephys_session_id)) do sesh
                _pspikes = @views allsesh_pspikes[allsesh_pspikes.ecephys_session_id .== sesh,
                                                  :]
                ys = _pspikes.pairwise_phase_consistency
                @assert sum(isnan.(ys)) / length(ys) < 0.3
                idxs = .!isnan.(ys) .& .!ismissing.(ys) .& (ys .> 0) # Negative values come from low sample sizes
                xs = _pspikes.streamlinedepth[idxs] .|> Float32
                ys = ys[idxs] .|> Float32
                # hexbin(xs, ys)
                B = HistBins(xs; bins)
                _ms = rectify(B(ys), dims = :bin)
            end

            X = cat(mss...; dims = 2)
            X = map(eachrow(X)) do x
                x = vcat(x...)
            end
            _ms = pmap(bootstrapmean, X)
            # _ms = bootstrapmedian.(X)
            σs = last.(_ms)
            _ms = first.(_ms)
            _σl = first.(σs)
            _σh = last.(σs)
            # _ms, (_σl, _σh) = bootstrapmedian(X; dims = 2)
            idxs = .!isnan.(_ms) .& .!isnan.(_σl) .& .!isnan.(_σh)
            _ms = _ms[idxs]
            _σl = _σl[idxs]
            _σh = _σh[idxs]
            # _σs = map((args...) -> nansafe(confidence)(collect(args)) |> first, mss...)
            # _ms = map((args...) -> nansafe(mean)(collect(args)) |> first, mss...)

            _ls = lookup(_ms, 1)
            ms = upsample(_ms, 10)
            σl = upsample(_σl, 10)
            σh = upsample(_σh, 10)
            ls = lookup(ms, 1)
            band!(ax, ls, collect(σl), collect(σh);
                  color = (structurecolormap[structure], 0.3), label = structure)
            lines!(ax, ls, collect(ms), color = (structurecolormap[structure], 0.7),
                   label = structure)
            scatter!(ax, _ls, collect(_ms); color = structurecolormap[structure],
                     label = structure)
        end
        l = axislegend(ax, merge = true, nbanks = 2)
        plotlayerints!(ax, layerints; axis = :x, flipside = false, newticks = false,
                       bgcolor = Makie.RGBA(0, 0, 0, 0))
        reverselegend!(l)
    end
    begin # * Plot, for each structure, SAC
        @info "Plotting structure-wise spike--amplitude coupling"
        ax = Axis(gs[4], xlabel = "Cortical depth (%)", ylabel = "SAC",
                  title = "Spike-amplitude coupling (γ)",
                  limits = ((0.05, 0.95), (1.1, 1.5)),
                  xtickformat = depthticks)
        for structure in reverse(structures)
            idxs = pspikes.structure_acronym .== structure
            allsesh_pspikes = @views pspikes[idxs, :]
            mss = map(unique(allsesh_pspikes.ecephys_session_id)) do sesh
                _pspikes = @views allsesh_pspikes[allsesh_pspikes.ecephys_session_id .== sesh,
                                                  :]
                ys = _pspikes.spike_amplitude_coupling
                xs = _pspikes.streamlinedepth .|> Float32
                ys = ys .|> Float32
                B = HistBins(xs; bins)
                _ms = rectify(B(ys), dims = :bin)
            end

            X = cat(mss...; dims = 2)
            X = map(eachrow(X)) do x
                x = vcat(x...)
            end
            _ms = pmap(bootstrapmean, X)
            σs = last.(_ms)
            _ms = first.(_ms)
            _σl = first.(σs)
            _σh = last.(σs)
            idxs = .!isnan.(_ms) .& .!isnan.(_σl) .& .!isnan.(_σh)
            _ms = _ms[idxs]
            _σl = _σl[idxs]
            _σh = _σh[idxs]

            _ls = lookup(_ms, 1)
            ms = upsample(_ms, 10)
            σl = upsample(_σl, 10)
            σh = upsample(_σh, 10)
            ls = lookup(ms, 1)
            band!(ax, ls, collect(σl), collect(σh);
                  color = (structurecolormap[structure], 0.3), label = structure)
            lines!(ax, ls, collect(ms), color = (structurecolormap[structure], 0.7),
                   label = structure)
            scatter!(ax, _ls, collect(_ms); color = structurecolormap[structure],
                     label = structure)
        end
        l = axislegend(ax, merge = true, nbanks = 2)
        reverselegend!(l)
        plotlayerints!(ax, layerints; axis = :x, flipside = false, newticks = false,
                       bgcolor = Makie.RGBA(0, 0, 0, 0))
    end
    f
end

begin # * Preferred phases for spontaneous
    @info "Plotting preferred spike phases for spontaneous"
    ax = PolarAxis(gs[5]; theta_as_x = false, thetalimits = (0, 1.2pi),
                   rticks = 0:0.25:1, rtickformat = depthticks,
                   title = "Layerwise PPC angle")

    for structure in reverse(structures)
        idxs = unitdepths.structure_acronym .== structure
        allsesh_pspikes = @views unitdepths[idxs, :]
        mss = map(unique(allsesh_pspikes.ecephys_session_id)) do sesh
            _pspikes = @views allsesh_pspikes[allsesh_pspikes.ecephys_session_id .== sesh,
                                              :]
            _ys = _pspikes.spc
            ys = _pspikes.spc_angle
            xs = _pspikes.streamlinedepth .|> Float32
            ys = ys .|> Float32
            B = HistBins(xs; bins)
            _ms = rectify(B(ys), dims = :bin)
            c = rectify(B(_ys), dims = :bin)
            return _ms, c
        end
        css = last.(mss)
        mss = first.(mss)

        X = cat(mss...; dims = 2)
        X = map(eachrow(X)) do x
            x = vcat(x...)
        end
        _ms = pmap(Base.Fix1(bootstrapaverage, circularmean), X)
        σs = last.(_ms)
        _ms = first.(_ms)
        _σl = first.(σs)
        _σh = last.(σs)

        C = cat(mss...; dims = 2)
        C = map(eachrow(C)) do x
            x = vcat(x...)
        end
        _cs = pmap(bootstrapmean, C)
        _cs = first.(_cs)

        # _ms, (_σl, _σh) = bootstrapaverage(circularmean, X; dims = 2)
        idxs = .!isnan.(_ms) .& .!isnan.(_σl) .& .!isnan.(_σh)
        _ms = _ms[idxs]

        # _cs, _ = bootstrapmedian(C; dims = 2)
        _cs = _cs[idxs]

        _ls = lookup(_ms, 1)

        x = unwrap(_ms)
        x = upsample(x, 5)
        ms = SpatiotemporalMotifs.wrap.(x; domain = (-π, π))

        cs = upsample(_cs, 10)
        cs = MinMax(cs)(cs)

        c = seethrough(structurecolormap[structure])
        # band!(ax, lookup(mu, 1), l, h; color = muc |> collect, colormap = c, label = s)
        lines!(ax, lookup(ms, 1), collect(ms), color = cs |> collect, label = structure,
               colormap = c,
               linewidth = 7)
    end
end

begin # * Preferred phases
    @info "Plotting preferred spike phases"
    ax = PolarAxis(gs[6]; theta_as_x = false, thetalimits = (0, 1.2pi),
                   rticks = 0:0.25:1, rtickformat = depthticks,
                   title = "Layerwise PPC angle")

    for structure in reverse(structures)
        idxs = pspikes.structure_acronym .== structure
        allsesh_pspikes = @views pspikes[idxs, :]
        mss = map(unique(allsesh_pspikes.ecephys_session_id)) do sesh
            _pspikes = @views allsesh_pspikes[allsesh_pspikes.ecephys_session_id .== sesh,
                                              :]
            _ys = _pspikes.pairwise_phase_consistency
            ys = _pspikes.pairwise_phase_consistency_angle
            xs = _pspikes.streamlinedepth .|> Float32
            ys = ys .|> Float32
            B = HistBins(xs; bins)
            _ms = rectify(B(ys), dims = :bin)
            c = rectify(B(_ys), dims = :bin)
            return _ms, c
        end
        css = last.(mss)
        mss = first.(mss)

        X = cat(mss...; dims = 2)
        X = map(eachrow(X)) do x
            x = vcat(x...)
        end
        _ms = pmap(Base.Fix1(bootstrapaverage, circularmean), X)
        σs = last.(_ms)
        _ms = first.(_ms)
        _σl = first.(σs)
        _σh = last.(σs)

        C = cat(mss...; dims = 2)
        C = map(eachrow(C)) do x
            x = vcat(x...)
        end
        _cs = pmap(bootstrapmean, C)
        _cs = first.(_cs)

        # _ms, (_σl, _σh) = bootstrapaverage(circularmean, X; dims = 2)
        idxs = .!isnan.(_ms) .& .!isnan.(_σl) .& .!isnan.(_σh)
        _ms = _ms[idxs]

        # _cs, _ = bootstrapmedian(C; dims = 2)
        _cs = _cs[idxs]

        _ls = lookup(_ms, 1)

        x = unwrap(_ms)
        x = upsample(x, 5)
        ms = SpatiotemporalMotifs.wrap.(x; domain = (-π, π))

        cs = upsample(_cs, 10)
        cs = MinMax(cs)(cs)

        c = seethrough(structurecolormap[structure])
        # band!(ax, lookup(mu, 1), l, h; color = muc |> collect, colormap = c, label = s)
        lines!(ax, lookup(ms, 1), collect(ms), color = cs |> collect, label = structure,
               colormap = c,
               linewidth = 7)
    end
end

begin
    addlabels!(f, labelformat)
    linkyaxes!(gs[1] |> contents |> first, gs[2] |> contents |> first)
    wsave(plotdir("fig6", "spike_lfp.pdf"), f)
    f
end
