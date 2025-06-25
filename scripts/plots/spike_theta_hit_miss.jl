#! /bin/bash
# -*- mode: julia -*-
#=
exec julia +1.10.9 -t auto --color=yes "${BASH_SOURCE[0]}" "$@"
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

begin # * Compare distribution of PPC for hit/miss trials
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

begin # * Look at onset PPC
    onset = pspikes.hit_onset_pairwise_phase_consistency_pvalue
    offset = pspikes.miss_onset_pairwise_phase_consistency_pvalue
    onset = clamp.(onset, [1e-10 .. Inf])
    offset = clamp.(offset, [1e-10 .. Inf])
    onset = log10.(onset)
    offset = log10.(offset)
    hist(filter(!isnan, onset))
    hist!(filter(!isnan, offset))
    current_figure() |> display
end
begin # * Onset preferred phase
    hit_onset = pspikes.hit_onset_pairwise_phase_consistency
    miss_onset = pspikes.miss_onset_pairwise_phase_consistency
    onset = pspikes.onset_pairwise_phase_consistency
    hit_onset_phase = pspikes.hit_onset_pairwise_phase_consistency_angle
    miss_onset_phase = pspikes.miss_onset_pairwise_phase_consistency_angle

    sensitive_idxs = pspikes.pairwise_phase_consistency .> 0.25

    begin # * Only use neurons that have significant PPC's in both hit and miss cases
        ps = pspikes.hit_onset_pairwise_phase_consistency_pvalue
        ps[isnan.(ps)] .= 1.0
        ps[ps .< 0] .= 0.0
        ps = adjust(ps, BenjaminiHochberg())
        significant_idxs = ps .< 0.01

        ps = pspikes.miss_onset_pairwise_phase_consistency_pvalue
        ps[isnan.(ps)] .= 1.0
        ps[ps .< 0] .= 0.0
        ps = adjust(ps, BenjaminiHochberg())
        significant_idxs = significant_idxs .& (ps .< 0.01)
    end

    hit_onset_phase = hit_onset_phase[sensitive_idxs .& significant_idxs]
    miss_onset_phase = miss_onset_phase[sensitive_idxs .& significant_idxs]
    hist(filter(!isnan, hit_onset_phase), label = "hit onset phase")
    hist!(filter(!isnan, miss_onset_phase), label = "miss onset phase")
    axislegend(; position = :lt)
    current_figure() |> display
end

begin # * Polar histogram
    f = Figure()
    ax = PolarAxis(f[1, 1])
    polarhist!(ax, hit_onset_phase, bins = 30, label = "hit onset phase",
               color = (cornflowerblue, 0.5))
    polarhist!(ax, miss_onset_phase, bins = 30, label = "miss onset phase",
               color = (crimson, 0.5))

    f
end

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

begin # * Preferred phase, hit vs miss, across layers
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

begin # * How does firing rate vary with depth? Make a psth for each theta-sensitive neuron
    # ? H: Theta-sensitive neurons in superficial layers should have either delayed or
    # ? suppressed firing rates in hit trials

end
