#! /bin/bash
# -*- mode: julia -*-
#=
exec julia +1.10.10 -t auto --color=yes "${BASH_SOURCE[0]}" "$@"
=#
using DrWatson
@quickactivate "SpatiotemporalMotifs"

using SpatiotemporalMotifs
import TimeseriesTools: freqs
using Peaks
using FileIO
using SpatiotemporalMotifs
import SpatiotemporalMotifs.layers
using Random
@preamble
set_theme!(foresight(:physics))

# begin # * So we want to load in the i) theta phases for each ii) spike time at iii) each depth, then do some comparison of these phases across layers for hit vs miss trials...
#     stimulus = r"Natural_Images"
#     vars = [:ϕ, :r]

#     session_table = load(datadir("posthoc_session_table.jld2"), "session_table")
#     oursessions = session_table.ecephys_session_id

#     path = datadir("calculations")
#     Q = calcquality(path)[Structure = At(structures)][SessionID = At(oursessions)]
#     out = load_calculations(Q; stimulus, vars)
#     Qs = calcquality(datadir("power_spectra"))[Structure = At(structures)][SessionID = At(oursessions)]
#     unitdepths = load_unitdepths(Qs)

#     begin # * Format spikes to depth dataframe. These have already been rectified
#         spikes = map(out) do o
#             regionalspikes = map(o) do _o
#                 DataFrame(collect.([keys(_o[:spiketimes]), values(_o[:spiketimes])]),
#                           [:ecephys_unit_id, :spiketimes])
#             end
#             vcat(regionalspikes...)
#         end
#         D = vcat(spikes...)
#         spikes = innerjoin(D, unitdepths, on = :ecephys_unit_id)
#         filter!(:streamlinedepth => ∈(0 .. 1), spikes) # Cortex only
#     end
# end

begin # * Load trial-by-trial spc
    file = datadir("spike_lfp.jld2")
    pspikes = load(file, "pspikes")
end

begin # * Compare distribution of PPC for hit/miss trials. Thi should be moved to calculations. But be careful; we need to use the rectified trial times
    sessionids = unique(pspikes.ecephys_session_id)
    outfile = datadir("out&stimulus=Natural_Images.jld2")
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
end

begin # * Add unit layers and rectified change times
    pspikes.rectified_change_times = Vector{Vector{<:Quantity}}(undef, nrow(pspikes))
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
end

begin # * Add trial times to pspikes. Move to calculations.
    pspikes.trial_change_times = map(eachrow(pspikes)) do row
        isempty(row[:trial_pairwise_phase_consistency]) && return
        _trials = trials[SessionID = At(row[:ecephys_session_id])]
        return _trials.change_time_with_display_delay # but be careful...the spiek times are RECTIFIED
        #!!!! NEED TO USE NON_RECTIFIED SPIKE TIMES MFFFFF
    end
end

begin # * Plot mean layer-wise PPC for hit/miss trials
    sessionid = sessionids[1]
    x = pspikes[pspikes.ecephys_session_id .== sessionid, :]
    structure = "VISl"
    x = x[x.structure_acronym .== structure, :]
    x = x[0.0 .< x.streamlinedepth .< 0.4, :]

    ppcs = x.trial_pairwise_phase_consistency
    hitmiss = x.hitmiss
    hits = map(ppcs, hitmiss) do p, h
        filter(!isnan, p[h])
    end
    misses = map(ppcs, hitmiss) do p, h
        filter(!isnan, p[.!h])
    end
end
begin
    density(filter(!isnan, vcat(hits...)))
    vlines!(mean(filter(!isnan, vcat(hits...))), color = cornflowerblue, label = "mean hit")
    density!(filter(!isnan, vcat(misses...)))
    vlines!(mean(filter(!isnan, vcat(misses...))), color = crimson, label = "mean miss")
    axislegend(; position = :lt)
    current_axis().xlabel = "PPC"
    current_axis().ylabel = "Density"
    current_figure() |> display
end

# begin # * Look at onset PPC
#     onset = pspikes.hit_onset_pairwise_phase_consistency_pvalue
#     offset = pspikes.miss_onset_pairwise_phase_consistency_pvalue
#     onset = clamp.(onset, [1e-10 .. Inf])
#     offset = clamp.(offset, [1e-10 .. Inf])
#     onset = log10.(onset)
#     offset = log10.(offset)
#     hist(filter(!isnan, onset))
#     hist!(filter(!isnan, offset))
#     current_figure() |> display
# end
# begin # * Onset preferred phase
#     hit_onset = pspikes.hit_onset_pairwise_phase_consistency
#     miss_onset = pspikes.miss_onset_pairwise_phase_consistency
#     onset = pspikes.onset_pairwise_phase_consistency
#     hit_onset_phase = pspikes.hit_onset_pairwise_phase_consistency_angle
#     miss_onset_phase = pspikes.miss_onset_pairwise_phase_consistency_angle

#     sensitive_idxs = pspikes.pairwise_phase_consistency .> 0.25

#     begin # * Only use neurons that have significant PPC's in both hit and miss cases
#         ps = pspikes.hit_onset_pairwise_phase_consistency_pvalue
#         ps[isnan.(ps)] .= 1.0
#         ps[ps .< 0] .= 0.0
#         ps = adjust(ps, BenjaminiHochberg())
#         significant_idxs = ps .< 0.01

#         ps = pspikes.miss_onset_pairwise_phase_consistency_pvalue
#         ps[isnan.(ps)] .= 1.0
#         ps[ps .< 0] .= 0.0
#         ps = adjust(ps, BenjaminiHochberg())
#         significant_idxs = significant_idxs .& (ps .< 0.01)
#     end

#     hit_onset_phase = hit_onset_phase[sensitive_idxs .& significant_idxs]
#     miss_onset_phase = miss_onset_phase[sensitive_idxs .& significant_idxs]
#     hist(filter(!isnan, hit_onset_phase), label = "hit onset phase")
#     hist!(filter(!isnan, miss_onset_phase), label = "miss onset phase")
#     axislegend(; position = :lt)
#     current_figure() |> display
# end

# begin # * Polar histogram
#     f = Figure()
#     ax = PolarAxis(f[1, 1])
#     polarhist!(ax, hit_onset_phase, bins = 30, label = "hit onset phase",
#                color = (cornflowerblue, 0.5))
#     polarhist!(ax, miss_onset_phase, bins = 30, label = "miss onset phase",
#                color = (crimson, 0.5))

#     f
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

if false # * Preferred phase, hit vs miss, across layers
    @info "Plotting preferred spike phases"

    begin # * Set up figure
        f = Figure()
        layerints = load(datadir("plots", "grand_unified_layers.jld2"), "layerints")
        bins = range(0, 1, length = 11)
    end

    ax = PolarAxis(f[1, 1]; theta_as_x = false, thetalimits = (0, 2pi),
                   rticks = 0:0.25:1, rtickformat = depthticks,
                   title = "Layerwise PPC angle", rlimits = (0.0, 0.8))

    for structure in reverse(structures)
        idxs = pspikes.structure_acronym .== structure
        allsesh_pspikes = @views pspikes[idxs, :]
        allangles = map(unique(allsesh_pspikes.ecephys_session_id)) do sesh
            _pspikes = @views allsesh_pspikes[allsesh_pspikes.ecephys_session_id .== sesh,
                                              :]
            _ys = _pspikes.hit_onset_pairwise_phase_consistency
            ys = _pspikes.hit_onset_pairwise_phase_consistency_angle
            xs = _pspikes.streamlinedepth .|> Float32
            ys = ys .|> Float32
            B = HistBins(xs; bins)
            _ms = rectify(B(ys), dims = :bin)
            c = rectify(B(_ys), dims = :bin)
            return _ms, c
        end
        css = last.(allangles)
        allangles = first.(allangles)

        catppc = cat(allangles...; dims = 2)
        catppc = map(eachrow(catppc)) do x
            x = vcat(x...)
        end
        _ms = pmap(Base.Fix1(bootstrapaverage, circularmean), catppc)
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
    display(f)
end

begin # * Plot on the same polar axis the hit and miss angles for VISp only
    structure = "VISl"

    begin # * Set up figure
        f = Figure()
        layerints = load(datadir("plots", "grand_unified_layers.jld2"), "layerints")
        bins = range(0, 1, length = 11)
    end

    ax = PolarAxis(f[1, 1]; theta_as_x = false, thetalimits = (0, 2pi),
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

using ProgressMeter
begin # * How does firing rate vary with depth? Make a psth for each theta-sensitive neuron
    binwidth = 15u"ms"
    timebins = range(SpatiotemporalMotifs.INTERVAL, step = binwidth) |> ustripall |>
               intervals
    dt0 = ustrip.(extrema(SpatiotemporalMotifs.INTERVAL))

    trial_psth_rates = @showprogress map(eachrow(pspikes)) do row # Takes about 10 minutes. Breaks if distributed...
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

begin # * Plot the mean firing rate of theta-sensitive neurons. For hit vs miss across layers
    structure = "VISp"

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

# begin # * Group over depth bins
#     depthbins = 0.05:0.1:0.95
#     depths = subspikes.streamlinedepth
#     hit_psths = ToolsArray(subspikes.hit_psth, (Depth(depths),))
#     hit_psths = groupby(hit_psths, Depth => Bins(intervals(depthbins)))
#     hit_psths = set(hit_psths, Depth => mean.(lookup(hit_psths, Depth)))
#     hit_psths = map(hit_psths) do psth
#         mean(filter(!isnothing, psth))
#     end |> stack
#     hit_psths = rectify(hit_psths, dims = 𝑡)
#     hit_psths = rectify(hit_psths, dims = Depth)

#     miss_psths = ToolsArray(subspikes.miss_psth, (Depth(depths),))
#     miss_psths = groupby(miss_psths, Depth => Bins(intervals(depthbins)))
#     miss_psths = set(miss_psths, Depth => mean.(lookup(miss_psths, Depth)))
#     miss_psths = map(miss_psths) do psth
#         mean(filter(!isnothing, psth))
#     end |> stack
#     miss_psths = rectify(miss_psths, dims = 𝑡)
#     miss_psths = rectify(miss_psths, dims = Depth)
# end

# begin
#     f = Figure()
#     ax = Axis(f[1, 1], yreversed = true)
#     # N = ZScore((hit_psths[𝑡 = -0.05 .. 0.25] .+ miss_psths[𝑡 = -0.05 .. 0.25]) / 2,
#     #    dims = 1)
#     # contrast =
#     contrast = hit_psths .- miss_psths
#     # contrast = contrast[:, 1:(end - 1)]
#     heatmap!(ax, contrast, colorrange = (-10, 10), colormap = binarysunset)
#     f |> display
# end

# begin
#     D = 0.1
#     intt = -0.0 .. 0.2
#     f = Figure(size = (400, 800))
#     ax = Axis(f[1, 1])
#     lines!(ax, hit_psths[𝑡 = intt, Depth = Near(D)])
#     lines!(ax, miss_psths[𝑡 = intt, Depth = Near(D)], color = crimson)
#     f |> display

#     D = 1.0
#     ax = Axis(f[2, 1])
#     lines!(ax, hit_psths[𝑡 = intt, Depth = Near(D)])
#     lines!(ax, miss_psths[𝑡 = intt, Depth = Near(D)], color = crimson)
#     f |> display
# end

# begin # * Calculate firing rate grouped by layers...
#     # * select neurons in layer 1
#     ss = (subspikes.layer .== "L4") .| (subspikes.layer .== "L4")
#     layer2_spikes = subspikes[ss, :]
#     layer2_hit_psth = layer2_spikes.hit_psth |> mean
#     layer2_hit_psth = rectify(layer2_psth, dims = 𝑡)
#     layer2_miss_psth = filter(!isnothing, layer2_spikes.miss_psth) |> mean
#     layer2_miss_psth = rectify(layer2_miss_psth, dims = 𝑡)

#     ss = (subspikes.layer .== "L5") .| (subspikes.layer .== "L6")
#     layer6_spikes = subspikes[ss, :]
#     layer6_hit_psth = filter(!isnothing, layer6_spikes.hit_psth) |> mean
#     layer6_hit_psth = rectify(layer6_hit_psth, dims = 𝑡)
#     layer6_miss_psth = filter(!isnothing, layer6_spikes.miss_psth) |> mean
#     layer6_miss_psth = rectify(layer6_miss_psth, dims = 𝑡)

#     f = Figure(size = (400, 800))
#     ax = Axis(f[1, 1])
#     lines!(ax, layer2_hit_psth[𝑡 = intt], label = "L2/3 Hit")
#     lines!(ax, layer2_miss_psth[𝑡 = intt], label = "L2/3 Miss")
#     axislegend(; position = :lt)

#     ax = Axis(f[2, 1])
#     lines!(ax, layer6_hit_psth[𝑡 = intt], label = "L6 Hit")
#     lines!(ax, layer6_miss_psth[𝑡 = intt], label = "L6 Miss")
#     f |> display
# end

begin
    begin # * Average psths  +  hit–miss difference (z-scored)
        structure = "VISp"

        idxs = pspikes.structure_acronym .== [structure]
        subspikes = pspikes[idxs, :]

        # * Average psths ----------------------------------------------------------
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

        # * NEW: hit – miss, then per-neuron z-score -------------------------------
        subspikes.diff_psth_z = map(eachrow(subspikes)) do row
            h, m = row.hit_psth, row.miss_psth
            (h === nothing || m === nothing) && return nothing
            d = h .- m                                    # hit – miss
            μ, σ = mean(d), std(d)
            σ == 0.0f0 && return nothing                     # guard flat traces
            return (d .- μ) ./ σ                           # z-scored difference
        end
    end

    begin # * Plot hit – miss firing rate (one axis)
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
