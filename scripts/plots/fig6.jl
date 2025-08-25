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

if !isfile(calcdir("plots", savepath("fig6", Dict(), "jld2")))
    if nprocs() == 1
        if SpatiotemporalMotifs.CLUSTER()
            using USydClusters
            ourprocs = USydClusters.Physics.addprocs(6; mem = 33, ncpus = 8,
                                                     project = projectdir(),
                                                     queue = "l40s")
        else
            addprocs(6) # Equal to the number of structures
        end
    end
    @everywhere using SpatiotemporalMotifs
    @everywhere SpatiotemporalMotifs.@preamble
end

plot_data, data_file = produce_or_load(Dict(), calcdir("plots");
                                       filename = savepath("fig6")) do _
    stimulus = r"Natural_Images"
    outfile = calcdir("out&stimulus=$(stimulus).jld2")
    vars = [:ϕ, :r]

    session_table = load(calcdir("plots", "posthoc_session_table.jld2"), "session_table")
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

    begin # * Calculate spike--phase and spike--amplitude coupling across layers. Takes about 60 minutes over 6 workers
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
            structurepspikes = groupby(pspikes, :structure_acronym) |> deepcopy
            subouts = map(keys(structurepspikes)) do k
                view(out, findfirst(structures .== [k.structure_acronym]))
            end
            dframes = pmap(structurepspikes, keys(structurepspikes),
                           subouts) do dframe, k, O
                structure = k.structure_acronym
                @info "Calculating spike coupling for $(structure)"
                ϕ = map(only(O)) do o
                    o[:ϕ][𝑡 = SpatiotemporalMotifs.INTERVAL] # Phase during trial interval
                end
                r = map(only(O)) do o
                    o[:r][𝑡 = SpatiotemporalMotifs.INTERVAL] # Amplitude during trial interval
                end
                spc!(dframe, ustripall.(ϕ)) # * PPC spike--phase coupling
                sac!(dframe, ustripall.(r)) # * Mean normalized amplitude spike--amplitude coupling
                return dframe
            end
            pspikes = vcat(dframes...)
            tagsave(calcdir("spike_lfp.jld2"), Dict("pspikes" => pspikes))
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
        begin # * How does firing rate vary with depth?
            binwidth = 15u"ms"
            timebins = range(SpatiotemporalMotifs.INTERVAL, step = binwidth) |>
                       ustripall |>
                       intervals
            dt0 = ustrip.(extrema(SpatiotemporalMotifs.INTERVAL))

            @info "Calculating trial PSTH rates for each spike" # Takes about 30 minutes
            @withprogress name="PSTH" begin
                threadmax = size(pspikes, 1)
                threadlog = 0
                trial_psth_rates = map(eachrow(pspikes)) do row # ! Takes about 10 minutes. Breaks if distributed...
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
                    threadlog += 1
                    @logprogress threadlog / threadmax
                    return trial_psth_rates
                end
                pspikes.trial_psth_rates = trial_psth_rates
            end
        end

        @info "Calculating normalized PSTH"
        begin # * Do a proper hit/miss comparison
            @withprogress name="Normalized PSTH" begin
                threadmax = size(pspikes, 1)
                threadlog = 0
                pspikes.zscored_trial_psth_rates = map(eachrow(pspikes)) do row
                    if isnothing(row[:trial_psth_rates]) || isempty(row[:trial_psth_rates])
                        return
                    end
                    psths = row[:trial_psth_rates]
                    bpsths = filter(!isnothing, psths)
                    if isempty(bpsths)
                        return
                    end

                    # * Remove pre-stimulus baseline
                    # bpsth = map(bpsths) do psth
                    #     if !(psth isa RegularTimeSeries)
                    #         psth = rectify(psth, dims = 𝑡)
                    #     end
                    #     psth[𝑡 = (-0.25 .. 0.0)]
                    # end |> mean
                    tpsth = rectify.(bpsths, dims = 𝑡) |> mean # Mean rate across all time and trials (both hit and miss)
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
                    threadlog += 1
                    @logprogress threadlog / threadmax
                    return psths
                end
            end
        end

        @info "Calculating hit/miss PSTH"
        pspikes.hit_psth = map(eachrow(pspikes)) do row
            psths = row.zscored_trial_psth_rates
            if isnothing(psths)
                return nothing
            else
                psths = psths[row.hitmiss]
                psths = filter(!isnothing, psths)
                if isempty(psths)
                    return nothing
                else
                    return convert.(Float32, mean(psths))
                end
            end
        end
        pspikes.miss_psth = map(eachrow(pspikes)) do row
            psths = row.zscored_trial_psth_rates
            if isnothing(psths)
                return nothing
            else
                psths = psths[.!row.hitmiss]
                psths = filter(!isnothing, psths)
                if isempty(psths)
                    return nothing
                else
                    return convert.(Float32, mean(psths))
                end
            end
        end
    end

    begin # * Format hit and miss mean psths
        hit_psths = ToolsArray(pspikes.hit_psth, (Unit(pspikes.ecephys_unit_id),))
        hit_psths = filter(!isnothing, hit_psths)
        hit_psths = stack(hit_psths)

        miss_psths = ToolsArray(pspikes.miss_psth, (Unit(pspikes.ecephys_unit_id),))
        miss_psths = filter(!isnothing, miss_psths)
        miss_psths = stack(miss_psths)
    end

    pspikes = subset(pspikes, :pairwise_phase_consistency => ByRow(!isnan))
    pspikes = subset(pspikes, :pairwise_phase_consistency => ByRow(>(0)))
    pspikes = subset(pspikes, :spike_amplitude_coupling => ByRow(!isnan))
    select!(pspikes, Not([:spiketimes])) # We don't need the heavy spiketimes for this plot
    select!(pspikes, Not([:filtering]))
    select!(pspikes, Not([:trial_pairwise_phase_consistency]))
    select!(pspikes, Not([:trial_pairwise_phase_consistency_pvalue]))
    select!(pspikes, Not([:trial_spike_amplitude_coupling]))
    select!(pspikes, Not([:trial_psth_rates]))
    select!(pspikes, Not([:zscored_trial_psth_rates]))
    select!(pspikes, Not([:hit_psth]))
    select!(pspikes, Not([:miss_psth]))

    unitdepths = subset(unitdepths, :stimulus => ByRow(==("spontaneous")))
    unitdepths = subset(unitdepths, :spc => ByRow(!isnan))
    unitdepths = subset(unitdepths, :spc => ByRow(>(0)))
    unitdepths = subset(unitdepths, :sac => ByRow(!isnan))

    return (@strdict pspikes unitdepths hit_psths miss_psths)
end

begin # * Set up figure
    @unpack pspikes, unitdepths, hit_psths, miss_psths = plot_data
    begin # * Set up figure
        f = SixPanel()
        gs = subdivide(f, 3, 2)

        sf = TwoPanel()
        sgs = subdivide(sf, 1, 3)

        layerints = load(calcdir("plots", "grand_unified_layers.jld2"), "layerints")
        bins = range(0, 1, length = 11)
    end
end

begin
    begin # * Plot, for each structure, spontaneous SPC
        @info "Plotting spontaneous spike--phase coupling"
        ax = Axis(sgs[1], ylabel = "Cortical depth (%)", xlabel = "PPC",
                  title = "θ spike-phase coupling",
                  limits = ((0, nothing), (0, 1)),
                  ytickformat = depthticks,
                  xtickformat = terseticks,
                  yreversed = true)
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
            ms = upsample(_ms, 5)
            σl = upsample(_σl, 5)
            σh = upsample(_σh, 5)
            ls = lookup(ms, 1)
            band!(ax, Point2f.(collect(σl), lookup(ms, 1)),
                  Point2f.(collect(σh), lookup(ms, 1));
                  color = (structurecolormap[structure], 0.3), label = structure)
            lines!(ax, collect(ms), lookup(ms, 1),
                   color = (structurecolormap[structure], 0.7),
                   label = structure)
            scatter!(ax, collect(_ms), _ls; color = structurecolormap[structure],
                     label = structure)
        end
        l = axislegend(ax, merge = true, nbanks = 2, position = :rb, fontsize = 10)
        plotlayerints!(ax, layerints; axis = :y, flipside = true, newticks = false)
        reverselegend!(l)
    end

    begin # * Plot, for each structure, spontaneous SAC
        @info "Plotting spontaneous spike--amplitude coupling"
        ax = Axis(sgs[2], ylabel = "Cortical depth (%)", xlabel = "SAC",
                  title = "γ spike-amplitude coupling",
                  limits = ((1.1, 1.5), (0, 1)),
                  ytickformat = depthticks, xtickformat = terseticks, yreversed = true)
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
            ms = upsample(_ms, 5)
            σl = upsample(_σl, 5)
            σh = upsample(_σh, 5)
            ls = lookup(ms, 1)
            band!(ax, Point2f.(collect(σl), lookup(ms, 1)),
                  Point2f.(collect(σh), lookup(ms, 1));
                  color = (structurecolormap[structure], 0.3), label = structure)
            lines!(ax, collect(ms), lookup(ms, 1),
                   color = (structurecolormap[structure], 0.7),
                   label = structure)
            scatter!(ax, collect(_ms), _ls; color = structurecolormap[structure],
                     label = structure)
        end
        l = axislegend(ax, merge = true, nbanks = 2, position = :rb, fontsize = 10)
        plotlayerints!(ax, layerints; axis = :y, flipside = true, newticks = false)
        reverselegend!(l)
    end

    begin # * Plot, for each structure, SPC
        @info "Plotting structure-wise spike--phase coupling"
        ax = Axis(gs[1], ylabel = "Cortical depth (%)", xlabel = "PPC",
                  title = "Spike-phase coupling (θ)",
                  limits = ((0, nothing), (0, 1)),
                  ytickformat = depthticks,
                  xtickformat = terseticks,
                  yreversed = true)
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
            ms = upsample(_ms, 5)
            σl = upsample(_σl, 5)
            σh = upsample(_σh, 5)
            ls = lookup(ms, 1)
            band!(ax, Point2f.(collect(σl), lookup(ms, 1)),
                  Point2f.(collect(σh), lookup(ms, 1));
                  color = (structurecolormap[structure], 0.3), label = structure)
            lines!(ax, collect(ms), ls, color = (structurecolormap[structure], 0.7),
                   label = structure)
            scatter!(ax, collect(_ms), _ls; color = structurecolormap[structure],
                     label = structure)
        end
        l = axislegend(ax, merge = true, nbanks = 2, position = :rb)
        plotlayerints!(ax, layerints; axis = :y, flipside = true, newticks = false)
        reverselegend!(l)
    end
    begin # * Plot, for each structure, SAC
        @info "Plotting structure-wise spike--amplitude coupling"
        ax = Axis(gs[2], ylabel = "Cortical depth (%)", xlabel = "SAC",
                  title = "Spike-amplitude coupling (γ)",
                  limits = ((1.1, 1.5), (0, 1)),
                  ytickformat = depthticks,
                  xtickformat = terseticks,
                  yreversed = true)
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
            band!(ax, Point2f.(collect(σl), lookup(ms, 1)),
                  Point2f.(collect(σh), lookup(ms, 1));
                  color = (structurecolormap[structure], 0.3), label = structure)
            lines!(ax, collect(ms), ls, color = (structurecolormap[structure], 0.7),
                   label = structure)
            scatter!(ax, collect(_ms), _ls; color = structurecolormap[structure],
                     label = structure)
        end
        l = axislegend(ax, merge = true, nbanks = 2, position = :rb)
        reverselegend!(l)
        plotlayerints!(ax, layerints; axis = :y, flipside = true, newticks = false)
    end
    f
end

begin # * Preferred phases for spontaneous
    @info "Plotting preferred spike phases for spontaneous"
    ax = PolarAxis(sgs[3]; theta_as_x = false, thetalimits = (0, 1.2pi),
                   rticks = 0:0.5:1, rtickformat = depthticks,
                   title = "Layerwise PPC angle", alignmode = Outside())

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
begin
    begin # * Filter only for sensitive and significant neurons
        begin # * Only use neurons that have significant PPC's in both hit and miss cases
            ps = pspikes.hit_onset_pairwise_phase_consistency_pvalue
            ps = convert(Vector{Float64}, ps)
            ps[isnan.(ps)] .= 1.0
            ps[ps .< 0] .= 0.0
            ps = adjust(ps, BenjaminiHochberg())
            significant_idxs = ps .< SpatiotemporalMotifs.PTHR

            ps = pspikes.miss_onset_pairwise_phase_consistency_pvalue
            ps = convert(Vector{Float64}, ps)
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

        sensitive_units = unique(pspikes.ecephys_unit_id)
        sensitive_units = intersect(sensitive_units, lookup(hit_psths, Unit),
                                    lookup(miss_psths, Unit))
    end

    begin # * Plot on the same polar axis the hit and miss angles for VISp only
        mainstructure = "VISp"

        sspikes = pspikes

        bins = range(0, 1, length = 11)

        cols = [:hit_onset_pairwise_phase_consistency_angle,
            :miss_onset_pairwise_phase_consistency_angle]

        spsf = SixPanel()
        spgs = subdivide(spsf, 3, 2)
        for (i, structure) in enumerate(structures) # Generate supplementary figure
            subax = PolarAxis(spgs[i]; theta_as_x = false, thetalimits = (0, 2pi),
                              rticks = 0:0.25:1, rtickformat = depthticks,
                              title = "Layerwise PPC angle", rlimits = (0.0, 0.8))
            for col in cols
                idxs = sspikes.structure_acronym .== structure
                allsesh_sspikes = @views sspikes[idxs, :]
                allangles = map(unique(allsesh_sspikes.ecephys_session_id)) do sesh
                    _pspikes = @views allsesh_sspikes[allsesh_sspikes.ecephys_session_id .== sesh,
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

                lines!(subax, lookup(ms, 1), collect(ms), label = structure, linewidth = 7)
            end
        end
        display(spsf)
    end

    function plot_fr!(axs, i, structure)
        idxs = pspikes.structure_acronym .== [structure]
        subspikes = pspikes[idxs, :]
        vlines!.(axs, [[0.0]]; color = (:black, 0.4), linestyle = :dash)

        ics = zip(["Supra.", "Granular", "Infra."],
                  [["L1", "L2/3"], ["L4"], ["L5", "L6"]])

        for (compartment, clayers) in ics
            lidxs = indexin(clayers, "L" .* SpatiotemporalMotifs.layers)
            ss = subspikes.layer .∈ [clayers]
            layer_spikes = subspikes[ss, :]

            unitids = intersect(layer_spikes.ecephys_unit_id, lookup(hit_psths, Unit))
            hpsth = hit_psths[Unit = At(unitids)]
            hpsth = hpsth[𝑡 = intt]
            idxs = allequal.(eachslice(hpsth, dims = 2))
            hpsth = hpsth[:, .!idxs]

            _μ, (_σl, _σh) = bootstrapmedian(hpsth, dims = 2)

            μ = upsample(_μ, 5)
            σl = upsample(_σl, 5)
            σh = upsample(_σh, 5)

            ts = collect(lookup(μ, 𝑡))
            band!(axs[1], ts, parent(σl), parent(σh), label = compartment,
                  color = (mean(layercolors[lidxs]), 0.3))
            lines!(axs[1], μ, label = compartment,
                   color = mean(layercolors[lidxs]))
            scatter!(axs[1], _μ, color = mean(layercolors[lidxs]),
                     markersize = 15, label = compartment)

            unitids = intersect(layer_spikes.ecephys_unit_id, lookup(miss_psths, Unit))
            mpsth = miss_psths[Unit = At(unitids)]
            mpsth = mpsth[𝑡 = intt]
            idxs = allequal.(eachslice(mpsth, dims = 2))
            mpsth = mpsth[:, .!idxs]

            _μ, (_σl, _σh) = bootstrapmedian(mpsth, dims = 2)

            μ = upsample(_μ, 5)
            σl = upsample(_σl, 5)
            σh = upsample(_σh, 5)

            ts = collect(lookup(μ, 𝑡))
            band!(axs[2], ts, parent(σl), parent(σh), label = compartment,
                  color = (mean(layercolors[lidxs]), 0.3))
            lines!(axs[2], μ, label = compartment,
                   color = mean(layercolors[lidxs]))
            scatter!(axs[2], _μ, color = mean(layercolors[lidxs]),
                     markersize = 15, label = compartment)

            begin # * Hit-miss contrast
                commonunits = intersect(lookup(hpsth, Unit), lookup(mpsth, Unit))
                chpsth = hpsth[Unit = At(commonunits)]
                cmpsth = mpsth[Unit = At(commonunits)]
                cpsth = chpsth - cmpsth
                𝑝 = map(eachslice(cpsth, dims = 𝑡)) do x
                    x = convert(Vector{Float64}, x)
                    pvalue(HypothesisTests.SignedRankTest(x), tail = :right)
                end
                𝑝[:] .= MultipleTesting.adjust(collect(𝑝), BenjaminiHochberg())
                begin # * Plot
                    c = cpsth .+ eps() .* randn(size(cpsth))
                    _μ, (_σl, _σh) = bootstrapmedian(c, dims = 2)

                    μ = upsample(_μ, 5)
                    σl = upsample(_σl, 5)
                    σh = upsample(_σh, 5)

                    ts = collect(lookup(μ, 𝑡))
                    band!(axs[3], ts, parent(σl), parent(σh), label = compartment,
                          color = (mean(layercolors[lidxs]), 0.3))
                    lines!(axs[3], μ, label = compartment,
                           color = mean(layercolors[lidxs]))

                    sigs = 𝑝 .< SpatiotemporalMotifs.PTHR
                    scatter!(axs[3], _μ[sigs], color = mean(layercolors[lidxs]),
                             markersize = 10, label = compartment, strokewidth = 3,
                             strokecolor = mean(layercolors[lidxs]))
                    scatter!(axs[3], _μ[.!sigs], strokecolor = mean(layercolors[lidxs]),
                             markersize = 10, strokewidth = 3,
                             color = :white)
                end
            end
        end

        axislegend(axs[1]; position = :rt, merge = true, fontsize = 10)
        linkyaxes!([axs[1], axs[2]])
        return axs
    end

    begin # * Plot hit & miss firing rate in different layers (onset)
        mainstructure = "VISp"

        intt = -0.05 .. 0.2

        ax = Axis(f[3, :][1, 1], title = "Hit trials ($mainstructure)", xlabel = "Time (s)",
                  ylabel = "Normalized firing-rate",
                  limits = ((-0.01, 0.185), nothing),
                  xtickformat = terseticks,
                  ytickformat = terseticks)
        ax2 = Axis(f[3, :][1, 2], title = "Miss trials ($mainstructure)",
                   xlabel = "Time (s)",
                   ylabel = "Normalized firing-rate",
                   limits = ((-0.01, 0.185), nothing),
                   xtickformat = terseticks,
                   ytickformat = terseticks)
        ax3 = Axis(f[3, :][1, 3], title = "Hit - miss ($mainstructure)",
                   xlabel = "Time (s)",
                   ylabel = "Normalized firing rate change",
                   limits = ((-0.01, 0.185), nothing),
                   xtickformat = terseticks,
                   ytickformat = terseticks)

        plot_fr!([ax, ax2, ax3], 1, mainstructure)

        display(f)
    end
end

begin # * Hit/miss firing rates for each structure
    frsf = Figure(size = (1440, 720) .* 1.25)
    frgs = subdivide(frsf, 3, 2)

    map(enumerate(structures)) do (i, structure)
        subgs = frgs[i]

        begin # * Plot hit & miss firing rate in different layers (onset)
            idxs = pspikes.structure_acronym .== [structure]
            subspikes = pspikes[idxs, :]

            intt = -0.03 .. 0.2

            ax = Axis(subgs[1, 1], title = "Hit trials ($structure)", xlabel = "Time (s)",
                      ylabel = "Normalized firing rate",
                      limits = ((-0.01, 0.185), nothing),
                      xtickformat = terseticks,
                      ytickformat = terseticks)
            ax2 = Axis(subgs[1, 2], title = "Miss trials ($structure)", xlabel = "Time (s)",
                       ylabel = "Normalized firing rate",
                       limits = ((-0.01, 0.185), nothing),
                       xtickformat = terseticks,
                       ytickformat = terseticks)
            ax3 = Axis(subgs[1, 3], title = "Hit - miss ($structure)",
                       xlabel = "Time (s)",
                       ylabel = "Normalized firing rate change",
                       xtickformat = terseticks)

            plot_fr!([ax, ax2, ax3], i, structure)
        end
    end
    addlabels!(frgs, frsf, labelformat, recurse = [])
    display(frsf)
end

begin # * Plot on the same polar axis the hit and miss angles for VISp only
    mainstructure = "VISp"

    sspikes = pspikes
    # sunits = intersect(sensitive_units, sspikes.ecephys_unit_id)

    bins = range(0.1, 1, length = 6)

    cols = [:hit_onset_pairwise_phase_consistency_angle,
        :miss_onset_pairwise_phase_consistency_angle]

    ms = map(structures) do structure
        idxs = sspikes.structure_acronym .== structure
        allsesh_sspikes = @views sspikes[idxs, :]
        allangles = map(unique(allsesh_sspikes.ecephys_session_id)) do sesh
            _pspikes = @views allsesh_sspikes[allsesh_sspikes.ecephys_session_id .== sesh,
                                              :]
            hys = _pspikes[:, :hit_onset_pairwise_phase_consistency_angle]
            mys = _pspikes[:, :miss_onset_pairwise_phase_consistency_angle]
            ys = phasegrad(hys, mys)
            # xs = _pspikes.streamlinedepth .|> Float32
            xs = _pspikes.layer
            ys = ys .|> Float32
            # B = HistBins(xs; bins)
            # _ms = rectify(B(ys), dims = :bin)
            # ys = ToolsArray(ys, (Depth(xs),))
            _layers = 'L' .* SpatiotemporalMotifs.layers[2:end]
            _ms = map(_layers) do l
                is = xs .== [l]
                ys[is]
            end
            _ms = ToolsArray(_ms, (Depth(_layers),))
            return _ms
        end
        catppc = cat(allangles...; dims = 2)
        catppc = map(eachrow(catppc)) do x
            x = vcat(x...)
        end

        _ms = map(nansafe(circularmean), catppc) .|> only
        ms = _ms

        # * Do a test about significant differences from 0
        𝑝 = map(catppc) do x
            x = filter(!isnan, x)
            HypothesisTests.pvalue(HypothesisTests.SignedRankTest(Float64.(x)))
        end

        return ms, 𝑝
    end
    𝑝 = ToolsArray(last.(ms), (Structure(structures),)) |> stack
    𝑝[:] .= MultipleTesting.adjust(𝑝[:], BenjaminiHochberg())
    ms = ToolsArray(first.(ms), (Structure(structures),)) |> stack
end
begin # * Heatmap figure
    ax = Axis(gs[4][1, 1], ylabel = "Cortical layer",
              yticks = (1:4, SpatiotemporalMotifs.layers[2:end]),
              xticks = (1:6, structures), xticklabelrotation = pi / 8, yreversed = true,
              limits = ((0.35, 6.5), (0.5, 4.5)),
              title = "PPC angle contrast")
    colorrange = (-2 * pi / 4, 2 * pi / 4)
    function degs(x::Number)
        return string(round(Int, x * 180 / π), "°")
    end
    degs(x) = map(degs, x)

    # * Add vertical lines corresponding to the layers
    for l in 1:size(ms, Depth)
        c = SpatiotemporalMotifs.layercolors[l + 1]
        lines!(ax, [0.41, 0.41], [l - 0.5, l + 0.5];
               color = c, linewidth = 30, linecap = :butt)
    end

    p = heatmap!(ax, 1:size(ms, 2), 1:size(ms, 1), collect(ms)'; colormap = darksunset,
                 colorrange, lowclip = first(darksunset), highclip = last(darksunset))
    pp = set(𝑝, Structure => 1:size(𝑝, Structure), Depth => 1:size(𝑝, Depth))'
    sigs = findall(pp .< SpatiotemporalMotifs.PTHR)
    sigs = collect(Iterators.product(lookup(pp)...))[sigs] .|> Point2f

    ywidth = 1 #mean(diff(lookup(ms, 1)))
    xwidth = 1
    function box_corners(center::Point2f, xwidth::Real, ywidth::Real)
        x, y = center
        half_x = xwidth / 2
        half_y = ywidth / 2

        return [
                Point2f(x - half_x, y - half_y),  # Bottom-left
                Point2f(x + half_x, y - half_y),  # Bottom-right
                Point2f(x + half_x, y + half_y),  # Top-right
                Point2f(x - half_x, y + half_y)   # Top-left
                ]
    end
    boxes = map(sigs) do s
        box_corners(s, xwidth, ywidth)
    end
    poly!.(ax, boxes, color = :transparent, strokecolor = :white, strokewidth = 2)
    # plotlayerints!(ax, layerints; axis = :y, flipside = true, newticks = false)
    Colorbar(gs[4][1, 2], p; label = "Phase difference (hit - miss)",
             tickformat = degs, ticks = (-pi / 2):(pi / 4):(pi / 2), alignmode = Outside())
end

begin
    addlabels!(f, labelformat)
    addlabels!(sf, labelformat)
    wsave(plotdir("fig6", "spike_lfp.pdf"), f)
    wsave(plotdir("fig6", "spike_lfp_spontaneous.pdf"), sf)
    wsave(plotdir("fig6", "firing_rates.pdf"), frsf)
    f
end
