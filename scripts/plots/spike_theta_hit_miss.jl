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

#     session_table = load(calcdir("plots", "posthoc_session_table.jld2"), "session_table")
#     oursessions = session_table.ecephys_session_id

#     path = calcdir("calculations")
#     Q = calcquality(path)[Structure = At(structures)][SessionID = At(oursessions)]
#     out = load_calculations(Q; stimulus, vars)
#     Qs = calcquality(calcdir("power_spectra"))[Structure = At(structures)][SessionID = At(oursessions)]
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
    file = calcdir("spike_lfp.jld2")
    pspikes = load(file, "pspikes")
end

# begin # * Compare distribution of PPC for hit/miss trials. Thi should be moved to calculations. But be careful; we need to use the rectified trial times
#     sessionids = unique(pspikes.ecephys_session_id)
#     outfile = calcdir("out&stimulus=Natural_Images.jld2")
#     trials = jldopen(outfile, "r") do f
#         map(sessionids) do sessionid
#             trials = f["VISp/$(sessionid)/trials"]
#         end
#     end
#     trials = ToolsArray(trials, (SessionID(sessionids),))
#     pspikes.hitmiss = map(eachrow(pspikes)) do row
#         isempty(row[:trial_pairwise_phase_consistency]) && return
#         _trials = trials[SessionID = At(row[:ecephys_session_id])]
#         return _trials.hit
#     end
# end

# begin # * Add unit layers and rectified change times
#     pspikes.rectified_change_times = Vector{Vector{<:Quantity}}(undef, nrow(pspikes))
#     pspikes.layer = Vector{String}(undef, nrow(pspikes))

#     sessionids = unique(pspikes.ecephys_session_id)
#     structs = unique(pspikes.structure_acronym)
#     jldopen(outfile, "r") do f
#         map(sessionids) do sessionid
#             map(structs) do structure
#                 lfp = f["$(structure)/$(sessionid)/V"]
#                 rectified_change_times = lookup(lfp, :changetime) |> collect
#                 idxs = pspikes.ecephys_session_id .== sessionid
#                 idxs = idxs .& (pspikes.structure_acronym .== structure)
#                 pspikes[idxs, :rectified_change_times] .= [rectified_change_times]

#                 layermap = f["$(structure)/$(sessionid)/layernames"]
#                 map(findall(idxs)) do idx
#                     depth = pspikes[idx, :probedepth]
#                     layer = layermap[Depth = Near(depth)]
#                     if layer ∈ ["or", "scwm", "cing"]
#                         i = findlast(lookup(layermap, 1) .< depth) |> last
#                         layer = layermap[i] # As in Unify Calculations.
#                     end
#                     pspikes[idx, :layer] = "L" *
#                                            SpatiotemporalMotifs.layers[parselayernum(layer)]
#                 end
#             end
#         end
#     end
# end

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
