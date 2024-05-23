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

stimulus = r"Natural_Images"
vars = [:r, :ϕ]

session_table = load(datadir("session_table.jld2"), "session_table")
oursessions = session_table.ecephys_session_id

path = datadir("calculations")
Q = calcquality(path)[structure = At(structures)]

config = @strdict stimulus vars
data, file = produce_or_load(produce_out(Q), config, datadir(); filename = savepath,
                             prefix = "out")
out = data["out"]

unitdepths = load_unitdepths(oursessions; rewrite = false)

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

begin # * Calculate spike--phase and spike--amplitude coupling across layers. # * Add preferred phase annotations
    if isfile(datadir("spike_lfp.jld2"))
        pspikes = load(datadir("spike_lfp.jld2"), "pspikes") # Delete this file to recalculate
    else
        pspikes = deepcopy(spikes)
    end
    requiredcols = [:pairwise_phase_consistency, :trial_pairwise_phase_consistency,
        :pairwise_phase_consistency_pvalue, :trial_pairwise_phase_consistency_pvalue,
        :spike_amplitude_coupling, :trial_spike_amplitude_coupling]
    if !all(hasproperty.([pspikes], requiredcols))
        pbar = ProgressBar(; refresh_rate = 5)
        job = addjob!(pbar; N = length(structures), description = "Structures")
        pjobs = [addjob!(pbar; N = length(o), description = "Session (phase)")
                 for o in out]
        rjobs = [addjob!(pbar; N = length(o), description = "Session (amplitude)")
                 for o in out]
        with(pbar) do
            for s in eachindex(structures)
                ϕ = getindex.(out[s], :ϕ)
                r = getindex.(out[s], :r)
                spc!(pspikes, ustripall.(ϕ); job = pjobs[s]) # * PPC spike--phase coupling
                sac!(pspikes, ustripall.(r); job = rjobs[s]) # * Mean normalized amplitude spike--amplitude coupling
                update!(job)
            end
        end
        save(datadir("spike_lfp.jld2"), "pspikes", pspikes)
    end
end

begin
    begin # * Set up figure
        f = TwoPanel()
        gs = subdivide(f, 1, 2)
    end

    if false # * Plot, for every structure
        ax = Axis(gs[1])
        ys = pspikes.pairwise_phase_consistency
        @assert sum(isnan.(ys)) / length(ys) < 0.05
        idxs = .!isnan.(ys) .& .!ismissing.(ys) .& (ys .> 0) # Negative values come from low sample sizes
        xs = pspikes.streamlinedepth[idxs] .|> Float32
        ys = ys[idxs] .|> Float32
        # hexbin(xs, ys)
        B = Bins(xs; bins = 20)
        ms = rectify(B(ys), dims = :bin) .|> mean
        lines!(ax, upsample(ms, 10))
    end

    begin # * Plot, for each structure
        bins = range(0, 1, length = 11)

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
                B = Bins(xs; bins)
                _ms = rectify(B(ys), dims = :bin) .|> mean
            end
            _σs = map((args...) -> nansafe(confidence)(collect(args)) |> first, mss...)
            _ms = map((args...) -> nansafe(mean)(collect(args)) |> first, mss...)

            _ls = lookup(_ms, 1)
            ms = upsample(_ms, 10)
            σs = upsample(_σs, 10)
            ls = lookup(ms, 1)
            band!(ax, ls, collect(ms .- σs), collect(ms .+ σs);
                  color = (structurecolormap[structure], 0.3), label = structure)
            lines!(ax, ls, ms, color = (structurecolormap[structure], 0.7),
                   label = structure)
            scatter!(ax, _ls, _ms; color = structurecolormap[structure], label = structure)
        end
        l = axislegend(ax, merge = true, nbanks = 2)
        reverselegend!(l)
    end
    begin # * Plot, for each structure
        bins = range(0, 1, length = 11)

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
                @assert sum(isnan.(ys)) / length(ys) < 0.3
                idxs = .!isnan.(ys) .& .!ismissing.(ys) .& (ys .> 0) # Negative values come from low sample sizes
                xs = _pspikes.streamlinedepth[idxs] .|> Float32
                ys = ys[idxs] .|> Float32
                # hexbin(xs, ys)
                B = Bins(xs; bins)
                _ms = rectify(B(ys), dims = :bin) .|> mean
            end

            _σs = map((args...) -> nansafe(confidence)(collect(args)) |> first, mss...)
            _ms = map((args...) -> nansafe(mean)(collect(args)) |> first, mss...)

            _ls = lookup(_ms, 1)
            ms = upsample(_ms, 10)
            σs = upsample(_σs, 10)
            ls = lookup(ms, 1)
            band!(ax, ls, collect(ms .- σs), collect(ms .+ σs);
                  color = (structurecolormap[structure], 0.3), label = structure)
            lines!(ax, ls, ms, color = (structurecolormap[structure], 0.7),
                   label = structure)
            scatter!(ax, _ls, _ms; color = structurecolormap[structure], label = structure)
        end
        l = axislegend(ax, merge = true, nbanks = 2)
        reverselegend!(l)
    end
    addlabels!(f)
    wsave(plotdir("spike_lfp", "spike_lfp.pdf"), f)
    f
end
