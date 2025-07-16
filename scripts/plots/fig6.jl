#! /bin/bash
#=
exec julia +1.10.10 -t auto "${BASH_SOURCE[0]}" "$@"
=#
using DrWatson
@quickactivate "SpatiotemporalMotifs"

import DimensionalData: metadata
using SpatiotemporalMotifs
@preamble
set_theme!(foresight(:physics))
Random.seed!(42)

plot_data, data_file = produce_or_load(Dict(), calcdir("plots");
                                       filename = savepath,
                                       prefix = "fig6") do _
    stimulus = r"Natural_Images"
    outfile = calcdir("out&stimulus=$(stimulus).jld2")
    vars = [:ϕ, :r]

    session_table = load(calcdir("posthoc_session_table.jld2"), "session_table")
    oursessions = session_table.ecephys_session_id

    path = calcdir("calculations")
    Q = calcquality(path)[Structure = At(structures)][SessionID = At(oursessions)]
    out = load_calculations(Q; stimulus, vars)
    Qs = calcquality(calcdir("power_spectra"))[Structure = At(structures)][SessionID = At(oursessions)]
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

    begin # * Add trial info to dataframe
        sessionids = unique(spikes.ecephys_session_id)
        outfile = calcdir("out&stimulus=Natural_Images.jld2")
        trials = jldopen(outfile, "r") do f
            map(sessionids) do sessionid
                trials = f["VISp/$(sessionid)/trials"]
            end
        end
        trials = ToolsArray(trials, (SessionID(sessionids),))
        spikes.hitmiss = map(eachrow(spikes)) do row
            # isempty(row[:trial_pairwise_phase_consistency]) && return
            _trials = trials[SessionID = At(row[:ecephys_session_id])]
            return _trials.hit
        end
    end

    begin # * Calculate spike--phase and spike--amplitude coupling across layers. Takes about 30 minutes over 64 cores, 125 GB
        # The idea here is that we have this larger spike_lfp file, and subset it to the
        # relevant data for this plot
        if isfile(calcdir("spike_lfp.jld2"))
            pspikes = load(calcdir("spike_lfp.jld2"), "pspikes") # Delete this file to recalculate
        else
            pspikes = deepcopy(spikes)
            idxs = pspikes.stimulus .== [r"Natural_Images"]
            pspikes = pspikes[idxs, :]
        end

        requiredcols = [:pairwise_phase_consistency,
            :trial_pairwise_phase_consistency,
            :pairwise_phase_consistency_pvalue, :trial_pairwise_phase_consistency_pvalue,
            :pairwise_phase_consistency_angle,
            :pairwise_phase_consistency_angle,
            :spike_amplitude_coupling,
            :trial_spike_amplitude_coupling,
            :onset_pairwise_phase_consistency,
            :offset_pairwise_phase_consistency,
            :onset_pairwise_phase_consistency_pvalue,
            :offset_pairwise_phase_consistency_pvalue,
            :onset_pairwise_phase_consistency_angle,
            :offset_pairwise_phase_consistency_angle,
            :hit_onset_pairwise_phase_consistency,
            :hit_offset_pairwise_phase_consistency,
            :hit_onset_pairwise_phase_consistency_pvalue,
            :hit_offset_pairwise_phase_consistency_pvalue,
            :hit_onset_pairwise_phase_consistency_angle,
            :hit_offset_pairwise_phase_consistency_angle,
            :miss_onset_pairwise_phase_consistency,
            :miss_offset_pairwise_phase_consistency,
            :miss_onset_pairwise_phase_consistency_pvalue,
            :miss_offset_pairwise_phase_consistency_pvalue,
            :miss_onset_pairwise_phase_consistency_angle,
            :miss_offset_pairwise_phase_consistency_angle,
            :hitmiss]

        if !all(hasproperty.([pspikes], requiredcols))
            for s in eachindex(structures)
                @info "Calculating spike coupling for $(structures[s])"
                ϕ = map(out[s]) do o
                    o[:ϕ][𝑡 = SpatiotemporalMotifs.INTERVAL] # Phase during trial interval
                end
                r = map(out[s]) do o
                    o[:r][𝑡 = SpatiotemporalMotifs.INTERVAL] # Amplitude during trial interval
                end
                spc!(pspikes, ustripall.(ϕ)) # * PPC spike--phase coupling
                sac!(pspikes, ustripall.(r)) # * Mean normalized amplitude spike--amplitude coupling
            end

            save(calcdir("spike_lfp.jld2"), "pspikes", pspikes)
        end
    end
    begin # * add trial info
        sessionids = unique(pspikes.ecephys_session_id)
        trials = jldopen(outfile, "r") do f
            map(sessionids) do sessionid
                trials = f["VISp/$(sessionid)/trials"]
            end
        end
        trials = ToolsArray(trials, (SessionID(sessionids),))
        pspikes.hitmiss = map(eachrow(pspikes)) do row
            isempty(row[:trial_pairwise_phase_consistency]) && return
            _trials = trials[SessionID = At(row[:ecephys_session_id])]
            return _trials.hit
        end

        begin # * Add unit layers and rectified change times
            pspikes.rectified_change_times = Vector{Vector{<:Quantity}}(undef,
                                                                        nrow(pspikes))
            pspikes.layer = Vector{String}(undef, nrow(pspikes))

            sessionids = unique(pspikes.ecephys_session_id)
            structs = unique(pspikes.structure_acronym)
            jldopen(outfile, "r") do f
                map(sessionids) do sessionid
                    map(structs) do structure
                        lfp = f["$(structure)/$(sessionid)/V"]
                        rectified_change_times = lookup(lfp, :changetime) |> collect
                        idxs = pspikes.ecephys_session_id .== sessionid
                        idxs = idxs .& (pspikes.structure_acronym .== structure)
                        pspikes[idxs, :rectified_change_times] .= [rectified_change_times]

                        layermap = f["$(structure)/$(sessionid)/layernames"]
                        map(findall(idxs)) do idx
                            depth = pspikes[idx, :probedepth]
                            layer = layermap[Depth = Near(depth)]
                            if layer ∈ ["or", "scwm", "cing"]
                                i = findlast(lookup(layermap, 1) .< depth) |> last
                                layer = layermap[i] # As in Unify Calculations.
                            end
                            pspikes[idx, :layer] = "L" *
                                                   SpatiotemporalMotifs.layers[parselayernum(layer)]
                        end
                    end
                end
            end

            pspikes.trial_change_times = map(eachrow(pspikes)) do row
                isempty(row[:trial_pairwise_phase_consistency]) && return
                _trials = trials[SessionID = At(row[:ecephys_session_id])]
                return _trials.rectified_change_times
            end
        end
        begin # * How does firing rate vary with depth? Make a psth for each theta-sensitive neuron
            binwidth = 15u"ms"
            timebins = range(SpatiotemporalMotifs.INTERVAL, step = binwidth) |>
                       ustripall |>
                       intervals
            dt0 = ustrip.(extrema(SpatiotemporalMotifs.INTERVAL))

            @info "Calculating trial PSTH rates for each spike"
            trial_psth_rates = map(eachrow(pspikes)) do row # Takes about 10 minutes. Breaks if distributed...
                isempty(row[:spiketimes]) && return
                spiketimes = row[:spiketimes]
                spiketimes = spiketrain(spiketimes)

                trial_spikes = map(row.rectified_change_times) do t
                    # spike aligned to stimulus onset time
                    dt = IntervalSets.Interval((dt0 .+ ustrip.(t))...)
                    ts = spiketimes[𝑡 = dt] # Extract
                    ts = set(ts, 𝑡 => times(ts) .- ustrip.(t)) # Align to t
                end
                trial_psth_rates = map(trial_spikes) do ts
                    if isempty(ts)
                        return nothing
                    else
                        psth = groupby(ts, 𝑡 => Bins(timebins)) .|> sum
                        psth_rates = psth ./ ustrip(uconvert(u"s", binwidth))
                        psth_rates = set(psth_rates, 𝑡 => mean.(times(psth_rates)))
                    end
                end

                # if isempty(trial_spikes)
                # psth_rates = nothing
                # else
                # psth_rates = mean(trial_spikes) ./ ustrip(binwidth)
                # psth_rates = set(psth_rates, 𝑡 => mean.(times(psth_rates)))
                # end
                return trial_psth_rates
            end
            pspikes.trial_psth_rates = trial_psth_rates
        end

        begin # * Do a proper hit/miss comparison by removing pre-stimulus baseline
            pspikes.zscored_trial_psth_rates = map(eachrow(pspikes)) do row
                isempty(row[:trial_psth_rates]) && return
                psths = row[:trial_psth_rates]
                bpsths = filter(!isnothing, psths)

                # * Remove pre-stimulus baseline
                # bpsth = map(bpsths) do psth
                #     if !(psth isa RegularTimeSeries)
                #         psth = rectify(psth, dims = 𝑡)
                #     end
                #     psth[𝑡 = (-0.25 .. 0.0)]
                # end |> mean
                tpsth = rectify.(bpsths, dims = 𝑡) |> mean # Just mean psth for whole trial

                μ = mean(tpsth)
                σ = std(tpsth)
                psths = map(psths) do psth
                    if isnothing(psth)
                        return nothing
                    elseif !(psth isa RegularTimeSeries)
                        psth = rectify(psth, dims = 𝑡)
                    end

                    psth = (psth .- μ) ./ σ

                    return psth
                end
            end
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

begin # * Set up figure
    @unpack pspikes, unitdepths = plot_data
    begin # * Set up figure
        f = SixPanel()
        gs = subdivide(f, 3, 2)

        sf = SixPanel()
        sgs = subdivide(sf[1, 1:2], 1, 3)
        sgss = subdivide(sf[2:3, 1:2], 3, 2)

        layerints = load(calcdir("plots", "grand_unified_layers.jld2"), "layerints")
        bins = range(0, 1, length = 11)
    end
end

begin
    begin # * Plot, for each structure, spontaneous SPC
        @info "Plotting spontaneous spike--phase coupling"
        ax = Axis(sgs[1], xlabel = "Cortical depth (%)", ylabel = "PPC",
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
        ax = Axis(sgs[2], xlabel = "Cortical depth (%)", ylabel = "SAC",
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
        ax = Axis(gs[1], xlabel = "Cortical depth (%)", ylabel = "PPC",
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
        ax = Axis(gs[2], xlabel = "Cortical depth (%)", ylabel = "SAC",
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
    ax = PolarAxis(sgs[3]; theta_as_x = false, thetalimits = (0, 1.2pi),
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

        C = cat(css...; dims = 2)
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
    ax = PolarAxis(gs[3]; theta_as_x = false, thetalimits = (0, 1.2pi),
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

# begin # * Plot mean layer-wise PPC for hit/miss trials
#     sessionid = sessionids[1]
#     x = pspikes[pspikes.ecephys_session_id .== sessionid, :]
#     structure = "VISl"
#     x = x[x.structure_acronym .== structure, :]
#     x = x[0.0 .< x.streamlinedepth .< 0.4, :]

#     ppcs = x.trial_pairwise_phase_consistency
#     hitmiss = x.hitmiss
#     hits = map(ppcs, hitmiss) do p, h
#         filter(!isnan, p[h])
#     end
#     misses = map(ppcs, hitmiss) do p, h
#         filter(!isnan, p[.!h])
#     end
# end
# begin
#     density(filter(!isnan, vcat(hits...)))
#     vlines!(mean(filter(!isnan, vcat(hits...))), color = cornflowerblue, label = "mean hit")
#     density!(filter(!isnan, vcat(misses...)))
#     vlines!(mean(filter(!isnan, vcat(misses...))), color = crimson, label = "mean miss")
#     axislegend(; position = :lt)
#     current_axis().xlabel = "PPC"
#     current_axis().ylabel = "Density"
#     current_figure() |> display
# end

begin # * Filter only for sensitive and significant neurons
    begin # * Only use neurons that have significant PPC's in both hit and miss cases
        ps = pspikes.hit_onset_pairwise_phase_consistency_pvalue
        ps[isnan.(ps)] .= 1.0
        ps[ps .< 0] .= 0.0
        ps = adjust(ps, BenjaminiHochberg())
        significant_idxs = ps .< SpatiotemporalMotifs.PTHR

        ps = pspikes.miss_onset_pairwise_phase_consistency_pvalue
        ps[isnan.(ps)] .= 1.0
        ps[ps .< 0] .= 0.0
        ps = adjust(ps, BenjaminiHochberg())
        significant_idxs = significant_idxs .& (ps .< SpatiotemporalMotifs.PTHR)

        pspikes = pspikes[significant_idxs, :]
    end

    begin # * Only neurons that are sensitive to phase
        cutoff_ppc = 0.1
        sensitive_idxs = pspikes.pairwise_phase_consistency .> cutoff_ppc
        pspikes = pspikes[sensitive_idxs, :]
    end
end

# if false # * Preferred phase, hit vs miss, across layers
#     @info "Plotting preferred spike phases"

#     begin # * Set up figure
#         f = Figure()
#         layerints = load(calcdir("plots", "grand_unified_layers.jld2"), "layerints")
#         bins = range(0, 1, length = 11)
#     end

#     ax = PolarAxis(f[1, 1]; theta_as_x = false, thetalimits = (0, 2pi),
#                    rticks = 0:0.25:1, rtickformat = depthticks,
#                    title = "Layerwise PPC angle", rlimits = (0.0, 0.8))

#     for structure in reverse(structures)
#         idxs = pspikes.structure_acronym .== structure
#         allsesh_pspikes = @views pspikes[idxs, :]
#         allangles = map(unique(allsesh_pspikes.ecephys_session_id)) do sesh
#             _pspikes = @views allsesh_pspikes[allsesh_pspikes.ecephys_session_id .== sesh,
#                                               :]
#             _ys = _pspikes.hit_onset_pairwise_phase_consistency
#             ys = _pspikes.hit_onset_pairwise_phase_consistency_angle
#             xs = _pspikes.streamlinedepth .|> Float32
#             ys = ys .|> Float32
#             B = HistBins(xs; bins)
#             _ms = rectify(B(ys), dims = :bin)
#             c = rectify(B(_ys), dims = :bin)
#             return _ms, c
#         end
#         css = last.(allangles)
#         allangles = first.(allangles)

#         catppc = cat(allangles...; dims = 2)
#         catppc = map(eachrow(catppc)) do x
#             x = vcat(x...)
#         end
#         _ms = pmap(Base.Fix1(bootstrapaverage, circularmean), catppc)
#         σs = last.(_ms)
#         _ms = first.(_ms)
#         _σl = first.(σs)
#         _σh = last.(σs)

#         C = cat(css...; dims = 2)
#         C = map(eachrow(C)) do x
#             x = vcat(x...)
#         end
#         _cs = pmap(bootstrapmean, C)
#         _cs = first.(_cs)

#         # _ms, (_σl, _σh) = bootstrapaverage(circularmean, X; dims = 2)
#         idxs = .!isnan.(_ms) .& .!isnan.(_σl) .& .!isnan.(_σh)
#         _ms = _ms[idxs]

#         # _cs, _ = bootstrapmedian(C; dims = 2)
#         _cs = _cs[idxs]

#         _ls = lookup(_ms, 1)

#         x = unwrap(_ms)
#         x = upsample(x, 5)
#         ms = SpatiotemporalMotifs.wrap.(x; domain = (-π, π))

#         cs = upsample(_cs, 10)
#         cs = MinMax(cs)(cs)

#         c = seethrough(structurecolormap[structure])
#         # band!(ax, lookup(mu, 1), l, h; color = muc |> collect, colormap = c, label = s)
#         lines!(ax, lookup(ms, 1), collect(ms), color = cs |> collect, label = structure,
#                colormap = c,
#                linewidth = 7)
#     end
#     display(f)
# end

begin # * Plot on the same polar axis the hit and miss angles for VISl only
    structure = "VISl"

    bins = range(0, 1, length = 11)

    ax = PolarAxis(gs[4]; theta_as_x = false, thetalimits = (0, 2pi),
                   rticks = 0:0.25:1, rtickformat = depthticks,
                   title = "Layerwise PPC angle", rlimits = (0.0, 0.8))
    cols = [:hit_onset_pairwise_phase_consistency_angle,
        :miss_onset_pairwise_phase_consistency_angle]
    for col in cols
        idxs = pspikes.structure_acronym .== structure
        allsesh_pspikes = @views pspikes[idxs, :]
        allangles = map(unique(allsesh_pspikes.ecephys_session_id)) do sesh
            _pspikes = @views allsesh_pspikes[allsesh_pspikes.ecephys_session_id .== sesh,
                                              :]
            ys = _pspikes[:, col]
            xs = _pspikes.streamlinedepth .|> Float32
            ys = ys .|> Float32
            B = HistBins(xs; bins)
            _ms = rectify(B(ys), dims = :bin)
            return _ms
        end

        catppc = cat(allangles...; dims = 2)
        catppc = map(eachrow(catppc)) do x
            x = vcat(x...)
        end
        _ms = pmap(Base.Fix1(bootstrapaverage, circularmean), catppc)
        σs = last.(_ms)
        _ms = first.(_ms)
        _σl = first.(σs)
        _σh = last.(σs)

        idxs = .!isnan.(_ms) .& .!isnan.(_σl) .& .!isnan.(_σh)
        _ms = _ms[idxs]

        _ls = lookup(_ms, 1)

        x = unwrap(_ms)
        x = upsample(x, 5)
        ms = SpatiotemporalMotifs.wrap.(x; domain = (-π, π))

        lines!(ax, lookup(ms, 1), collect(ms), label = structure, linewidth = 7)
    end
    display(f)
end

begin # * Plot the mean firing rate of theta-sensitive neurons. For hit vs miss across layers
    structure = "VISp" # !! Change to VISl

    idxs = pspikes.structure_acronym .== [structure]
    subspikes = pspikes[idxs, :]

    # * Average psths
    subspikes.hit_psth = map(eachrow(subspikes)) do row
        psths = row.zscored_trial_psth_rates
        psths = psths[row.hitmiss]
        psths = filter(!isnothing, psths)
        if isempty(psths)
            return nothing
        else
            return convert.(Float32, mean(psths))
        end
    end
    subspikes.miss_psth = map(eachrow(subspikes)) do row
        psths = row.zscored_trial_psth_rates
        psths = psths[.!row.hitmiss]
        psths = filter(!isnothing, psths)
        if isempty(psths)
            return nothing
        else
            return convert.(Float32, mean(psths))
        end
    end
end

begin # * Plot hit & miss firing rate in different layers (onset)
    intt = -0.02 .. 0.2
    f = Figure(size = (800, 300))
    ax = Axis(f[1, 1], title = "Hit trials", xlabel = "Time (s)",
              ylabel = "Normalized firing rate (Hz)")
    ax2 = Axis(f[1, 2], title = "Miss trials", xlabel = "Time (s)",
               ylabel = "Normalized firing rate (Hz)")

    ics = zip(["Supragranular", "Granular", "Infragranular"],
              [["L1", "L2/3"], ["L4"], ["L5", "L6"]])

    for (compartment, clayers) in ics
        lidxs = indexin(clayers, "L" .* SpatiotemporalMotifs.layers)
        ss = subspikes.layer .∈ [clayers]
        layer_spikes = subspikes[ss, :]
        hit_psth = filter(!isnothing, layer_spikes.hit_psth)
        hit_psth = rectify.(hit_psth, dims = 𝑡)
        hit_psth = ToolsArray(hit_psth, (Unit(1:size(hit_psth, 1)),)) |> stack
        hit_psth = hit_psth[𝑡 = intt]

        _μ, (_σl, _σh) = bootstrapmedian(hit_psth, dims = 2)

        # σ = std(hit_psth, dims = 2)
        # σ = dropdims(σ, dims = 2)
        # σl = μ .- σ / 2
        # σh = μ .+ σ / 2

        μ = upsample(_μ, 5)
        σl = upsample(_σl, 5)
        σh = upsample(_σh, 5)

        ts = collect(lookup(μ, 𝑡))
        band!(ax, ts, parent(σl), parent(σh), label = compartment,
              color = (mean(layercolors[lidxs]), 0.3))
        lines!(ax, μ, label = compartment,
               color = mean(layercolors[lidxs]))
        scatter!(ax, _μ, color = mean(layercolors[lidxs]),
                 markersize = 15, label = compartment)

        miss_psth = filter(!isnothing, layer_spikes.miss_psth)
        miss_psth = filter(x -> !any(isnan, x), miss_psth)
        miss_psth = ToolsArray(miss_psth, (Unit(1:size(miss_psth, 1)),)) |> stack
        miss_psth = miss_psth[𝑡 = intt]
        _μ, (_σl, _σh) = bootstrapmedian(miss_psth, dims = 2)

        μ = upsample(_μ, 5)
        σl = upsample(_σl, 5)
        σh = upsample(_σh, 5)

        ts = collect(lookup(μ, 𝑡))
        band!(ax2, ts, parent(σl), parent(σh), label = compartment,
              color = (mean(layercolors[lidxs]), 0.3))
        lines!(ax2, μ, label = compartment,
               color = mean(layercolors[lidxs]))
        scatter!(ax2, _μ, color = mean(layercolors[lidxs]),
                 markersize = 15, label = compartment)
    end

    axislegend(ax; position = :rt, merge = true)
    axislegend(ax2; position = :rt, merge = true)
    linkyaxes!([ax, ax2])
    display(f)
end
begin # * Plot hit & miss firing rate in different layers (offset
    intt = 0.23 .. 0.45
    f = Figure(size = (800, 300))
    ax = Axis(f[1, 1], title = "Hit trials", xlabel = "Time (s)",
              ylabel = "Normalized firing rate (Hz)")
    ax2 = Axis(f[1, 2], title = "Miss trials", xlabel = "Time (s)",
               ylabel = "Normalized firing rate (Hz)")

    ics = zip(["Supragranular", "Granular", "Infragranular"],
              [["L1", "L2/3"], ["L4"], ["L5", "L6"]])

    for (compartment, clayers) in ics
        lidxs = indexin(clayers, "L" .* SpatiotemporalMotifs.layers)
        ss = subspikes.layer .∈ [clayers]
        layer_spikes = subspikes[ss, :]
        hit_psth = filter(!isnothing, layer_spikes.hit_psth)
        hit_psth = rectify.(hit_psth, dims = 𝑡)
        hit_psth = ToolsArray(hit_psth, (Unit(1:size(hit_psth, 1)),)) |> stack
        hit_psth = hit_psth[𝑡 = intt]

        _μ, (_σl, _σh) = bootstrapmedian(hit_psth, dims = 2)

        # σ = std(hit_psth, dims = 2)
        # σ = dropdims(σ, dims = 2)
        # σl = μ .- σ / 2
        # σh = μ .+ σ / 2

        μ = upsample(_μ, 5)
        σl = upsample(_σl, 5)
        σh = upsample(_σh, 5)

        ts = collect(lookup(μ, 𝑡))
        band!(ax, ts, parent(σl), parent(σh), label = compartment,
              color = (mean(layercolors[lidxs]), 0.3))
        lines!(ax, μ, label = compartment,
               color = mean(layercolors[lidxs]))
        scatter!(ax, _μ, color = mean(layercolors[lidxs]),
                 markersize = 15, label = compartment)

        miss_psth = filter(!isnothing, layer_spikes.miss_psth)
        miss_psth = filter(x -> !any(isnan, x), miss_psth)
        miss_psth = ToolsArray(miss_psth, (Unit(1:size(miss_psth, 1)),)) |> stack
        miss_psth = miss_psth[𝑡 = intt]
        _μ, (_σl, _σh) = bootstrapmedian(miss_psth, dims = 2)

        μ = upsample(_μ, 5)
        σl = upsample(_σl, 5)
        σh = upsample(_σh, 5)

        ts = collect(lookup(μ, 𝑡))
        band!(ax2, ts, parent(σl), parent(σh), label = compartment,
              color = (mean(layercolors[lidxs]), 0.3))
        lines!(ax2, μ, label = compartment,
               color = mean(layercolors[lidxs]))
        scatter!(ax2, _μ, color = mean(layercolors[lidxs]),
                 markersize = 15, label = compartment)
    end

    axislegend(ax; position = :rt, merge = true)
    axislegend(ax2; position = :rt, merge = true)
    linkyaxes!([ax, ax2])
    display(f)
end

begin
    begin # * Average psths + hit-miss difference (z-scored)
        structure = "VISp"

        idxs = pspikes.structure_acronym .== [structure]
        subspikes = pspikes[idxs, :]
        subspikes.hit_psth = map(eachrow(subspikes)) do row
            psths = row.zscored_trial_psth_rates
            psths = psths[row.hitmiss]
            psths = filter(!isnothing, psths)
            isempty(psths) && return nothing
            return convert.(Float32, mean(psths))
        end

        subspikes.miss_psth = map(eachrow(subspikes)) do row
            psths = row.zscored_trial_psth_rates
            psths = psths[.!row.hitmiss]
            psths = filter(!isnothing, psths)
            isempty(psths) && return nothing
            return convert.(Float32, mean(psths))
        end

        subspikes.diff_psth_z = map(eachrow(subspikes)) do row
            h, m = row.hit_psth, row.miss_psth
            (h === nothing || m === nothing) && return nothing
            d = h .- m                                    # hit - miss
            μ, σ = mean(d), std(d)
            σ == 0.0f0 && return nothing                     # guard flat traces
            return (d .- μ) ./ σ                           # z-scored difference
        end
    end

    begin # * Plot hit - miss firing rate (one axis)
        intt = -0.0 .. 0.2
        f = Figure(size = (800, 300))
        ax = Axis(f[1, 1], title = "Hit − Miss (z)", xlabel = "Time (s)",
                  ylabel = "Δ rate (z-score)")

        ics = zip(["Supragranular", "Granular", "Infragranular"],
                  [["L1", "L2/3"], ["L4"], ["L5", "L6"]])

        for (compartment, clayers) in ics
            lidxs = indexin(clayers, "L" .* SpatiotemporalMotifs.layers)
            ss = subspikes.layer .∈ [clayers]
            layer_spikes = subspikes[ss, :]

            diff_psth = filter(!isnothing, layer_spikes.diff_psth_z)
            diff_psth = filter(x -> !any(isnan, x), diff_psth)
            diff_psth = ToolsArray(diff_psth, (Unit(1:size(diff_psth, 1)),)) |> stack
            diff_psth = diff_psth[𝑡 = intt]

            _μ, (_σl, _σh) = bootstrapmedian(diff_psth, dims = 2)

            μ = upsample(_μ, 5)
            σl = upsample(_σl, 5)
            σh = upsample(_σh, 5)

            ts = collect(lookup(μ, 𝑡))
            col = mean(layercolors[lidxs])

            band!(ax, ts, parent(σl), parent(σh), label = compartment,
                  color = (col, 0.3))
            lines!(ax, μ, label = compartment, color = col)
            scatter!(ax, _μ, color = col, markersize = 15, label = compartment)
        end

        axislegend(ax; position = :rt, merge = true)
        display(f)
    end
end

begin
    addlabels!(f, labelformat)
    addlabels!(sf, labelformat)
    wsave(plotdir("fig6", "spike_lfp.pdf"), f)
    wsave(plotdir("fig6", "spike_lfp_spontaneous.pdf"), sf)
    f
end
