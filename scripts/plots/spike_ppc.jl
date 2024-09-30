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
vars = [:ϕ, :r]

session_table = load(datadir("posthoc_session_table.jld2"), "session_table")
oursessions = session_table.ecephys_session_id

path = datadir("calculations")
Q = calcquality(path)[Structure = At(structures)][SessionID = At(oursessions)]
out = load_calculations(Q; stimulus, vars)
Qs = calcquality(datadir("power_spectra"))[Structure = At(structures)][SessionID = At(oursessions)]
unitdepths = load_unitdepths(Qs)

begin # * Plot spontaneous PPC
    D = subset(unitdepths, :stimulus => ByRow(==(r"Natural_Images")))
    D = subset(unitdepths, :spc => ByRow(!isnan))
    D = subset(unitdepths, :spc_pvalue => ByRow(<(1e-100)))
    x = ToolsArray(D.spc, Depth(D.streamlinedepth))
    bins = range(0, 1, length = 11)
    bins = [b[1] .. b[2] for b in zip(bins[1:(end - 1)], bins[2:end])]
    xs = DimensionalData.groupby(x, Depth => Bins(bins))
    xs = mean.(xs)
    xs = rectify(set(xs, Depth => mean.(lookup(xs, Depth))); dims = Depth)
    lines(decompose(xs)...)
    # ps = D.spc_pvalue
    # ps[ps .≤ 0] .= eps()
    # MultipleTesting.adjust(collect(ps), BenjaminiHochberg())
end

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

begin
    begin # * Set up figure
        f = FourPanel()
        gs = subdivide(f, 2, 2)
        layerints = load(datadir("grand_unified_layers.jld2"), "layerints")
    end

    if false # * Plot, for every structure
        ax = Axis(gs[1])
        ys = pspikes.pairwise_phase_consistency
        @assert sum(isnan.(ys)) / length(ys) < 0.05
        idxs = .!isnan.(ys) .& .!ismissing.(ys) .& (ys .> 0) # Negative values come from low sample sizes
        xs = pspikes.streamlinedepth[idxs] .|> Float32
        ys = ys[idxs] .|> Float32
        # hexbin(xs, ys)
        B = HistBins(xs; bins = 20)
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
                B = HistBins(xs; bins)
                _ms = rectify(B(ys), dims = :bin) .|> mean
            end

            X = cat(mss...; dims = 2)
            _ms, (_σl, _σh) = bootstrapmedian(X; dims = 2)
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
            lines!(ax, ls, ms, color = (structurecolormap[structure], 0.7),
                   label = structure)
            scatter!(ax, _ls, _ms; color = structurecolormap[structure], label = structure)
        end
        l = axislegend(ax, merge = true, nbanks = 2)
        plotlayerints!(ax, layerints; axis = :x, flipside = false, newticks = false,
                       bgcolor = Makie.RGBA(0, 0, 0, 0))
        reverselegend!(l)
    end
    begin # * Plot, for each structure
        bins = range(0, 1, length = 11)

        ax = Axis(gs[2], xlabel = "Cortical depth (%)", ylabel = "SAC",
                  title = "Spike-amplitude coupling (γ)",
                  limits = ((0.05, 0.95), (1.05, 1.45)),
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
                B = HistBins(xs; bins)
                _ms = rectify(B(ys), dims = :bin) .|> mean
            end

            X = cat(mss...; dims = 2)
            _ms, (_σl, _σh) = bootstrapmedian(X; dims = 2)
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
            lines!(ax, ls, ms, color = (structurecolormap[structure], 0.7),
                   label = structure)
            scatter!(ax, _ls, _ms; color = structurecolormap[structure], label = structure)
        end
        l = axislegend(ax, merge = true, nbanks = 2)
        reverselegend!(l)
        plotlayerints!(ax, layerints; axis = :x, flipside = false, newticks = false,
                       bgcolor = Makie.RGBA(0, 0, 0, 0))
    end
end

begin # * Preferred phases
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
            @assert sum(isnan.(ys)) / length(ys) < 0.3
            idxs = .!isnan.(_ys) .& .!ismissing.(_ys) .& (_ys .> 0) # Negative values come from low sample sizes
            xs = _pspikes.streamlinedepth[idxs] .|> Float32
            ys = ys[idxs] .|> Float32
            # hexbin(xs, ys)
            B = HistBins(xs; bins)
            _ms = rectify(B(ys), dims = :bin) .|> circularmean
            c = rectify(B(_ys), dims = :bin) .|> mean
            return _ms, c
        end
        css = last.(mss)
        mss = first.(mss)
        X = cat(mss...; dims = 2)
        C = cat(css...; dims = 2)

        _ms, (_σl, _σh) = bootstrapaverage(circularmean, X; dims = 2)
        idxs = .!isnan.(_ms) .& .!isnan.(_σl) .& .!isnan.(_σh)
        _ms = _ms[idxs]

        _cs, _ = bootstrapmedian(C; dims = 2)
        _cs = _cs[idxs]

        _ls = lookup(_ms, 1)

        x = unwrap(_ms)
        x = upsample(x, 5)
        ms = SpatiotemporalMotifs.wrap.(x; domain = (-π, π))

        cs = upsample(_cs, 10)
        cs = MinMax(cs)(cs)

        c = seethrough(structurecolormap[structure])
        # band!(ax, lookup(mu, 1), l, h; color = muc |> collect, colormap = c, label = s)
        lines!(ax, lookup(ms, 1), ms, color = cs |> collect, label = s, colormap = c,
               linewidth = 7)
    end
end

begin
    bins = range(0, 1, length = 11)
    ax = Axis(gs[4], xlabel = "Cortical depth (%)", ylabel = "Mean firing rate",
              title = "Firing rate",
              limits = ((0.05, 0.95), (0, nothing)),
              xtickformat = depthticks)
    for structure in reverse(structures)
        idxs = pspikes.structure_acronym .== structure
        allsesh_pspikes = @views pspikes[idxs, :]
        mss = map(unique(allsesh_pspikes.ecephys_session_id)) do sesh
            _pspikes = @views allsesh_pspikes[allsesh_pspikes.ecephys_session_id .== sesh,
                                              :]
            # ys = _pspikes.firing_rate .|> Float32
            ys = map(_pspikes.spiketimes) do s
                N = length(s) / (maximum(s) - minimum(s)) # Mean firing rate
            end .|> Float32
            xs = _pspikes.streamlinedepth .|> Float32
            # hexbin(xs, ys)
            B = HistBins(xs; bins)
            _ms = rectify(B(ys), dims = :bin) .|> mean
        end
        X = cat(mss...; dims = 2)
        _ms, (_σl, _σh) = bootstrapmedian(X; dims = 2)
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
        lines!(ax, ls, ms, color = (structurecolormap[structure], 0.7),
               label = structure)
        scatter!(ax, _ls, _ms; color = structurecolormap[structure], label = structure)
    end
    l = axislegend(ax, merge = true, nbanks = 2, position = :lt)
    reverselegend!(l)

    plotlayerints!(ax, layerints; axis = :x, flipside = false, newticks = false,
                   bgcolor = Makie.RGBA(0, 0, 0, 0))

    addlabels!(f)
    wsave(plotdir("spike_lfp", "spike_lfp.pdf"), f)
    f
end

# begin # * Spike MUA spectrum. Currently for stimulustask
#     dfs = DataFrames.groupby(spikes, :ecephys_session_id)

#     df = dfs[1]
#     df = df[df.structure_acronym .== "VISp", :]
#     S = AN.Session(df.ecephys_session_id |> unique |> only)
#     spontimes = AN.getspiketimes(S, df.id)
#     df.spontimes = getindex.([spontimes], df.id)
#     sort!(df, :streamlinedepth)
#     # df = df[df.streamlinedepth .> 0.8, :]
#     spiketrains = spiketrain.(df.spontimes)
#     epochs = AN.getepochs(S)
#     spontdur = epochs.interval[findlast(epochs.stimulus_name .== "spontaneous")]
#     df.spontimes = map(df.spontimes) do t
#         t = t[t .∈ [spontdur]]
#     end
#     if false
#         # * Calculate correlations
#         O = progressmap(x -> stoic(x...),
#                         Iterators.product(spiketrains, spiketrains);
#                         parallel = true)
#         O = reshape(O, length(spiketrains), length(spiketrains))
#         F = eigen(O)
#         w = F.vectors[:, end]
#         is = partialsortperm(w, 1:20; rev = true) # * Top 25 neurons
#     end
# end
# begin
#     dt = 0.001
#     ts = range(minimum(minimum.(spikes.spiketimes)), maximum(maximum.(spikes.spiketimes)),
#                step = dt)
#     ts = [a .. b for (a, b) in zip(ts[1:(end - 1)], ts[2:end])]

#     ss = progressmap(spiketrains; parallel = true) do s
#         length.(DimensionalData.groupby(s, 𝑡 =>Bins(ts)))
#     end

#     is = partialsortperm(sum.(ss), 1:10; rev = true)
#     a = ss[is[1]]
#     b = ss[is[2]]
#     y = crosscor(a, b, 1:Int(1 ÷ dt))
#     lines(range(dt, step = dt, length = length(y)), y)

#     # s = sum(ss)
#     # s = set(s, 𝑡 =>mean.(times(s)))
#     # s = rectify(s, dims = 𝑡)
#     # S = spectrum(s, 0.05; padding = 100)
#     # lines((TimeseriesTools.freqs(S)[2:end]), (parent(S)[2:end]);
#     #       axis = (; xscale = log10, yscale = log10, limits = ((1, nothing), nothing)))
#     # * Now convert to mua
#     # spiketrains = map(spiketrains) do s
#     #     a = DimensionalData.groupby(s, 𝑡 =>Bins(ts))
#     # end
#     # bs = progressmap(collect(values(spikes))[1:100]; parallel = true) do x
#     #     x = x[minimum(ts) .< x .< maximum(ts)]
#     #     b = B(x)(x)
#     #     b = length.(b)
#     # end
#     # b = set(bs[1], sum(bs))
#     # b = set(b, Dim{:bin} => 𝑡(ts[1:(end - 1)]))

#     # s = log10.(spectrum(b, 0.5)[𝑓 = 1 .. 50])
#     # lines(freqs(s), s)
#     # ax = current_axis()
#     # ax.xlabel = "Hz"
#     # ax.ylabel = "power"
#     # ax.title = "MUA spectrum"
#     # current_figure()
# end

# begin # * MUA spectrum of spontaneous. Is there theta in spiking, suggesting this is not a volume conduction from hippocampus?
#     sessionid = 1108528422
#     S = AN.Session(sessionid)
#     df = AN.getepochs(S)
# end
# begin
#     structure = "VISpm"
#     spontimes = AN.getspiketimes(S, structure)
#     epoch = df[4, :].interval
#     map(collect(keys(spontimes))) do id
#         t = spontimes[id]
#         t = t[t .∈ [epoch]]
#         spontimes[id] = t
#     end
#     spiketrains = spiketrain.(values(spontimes))
#     spiketrains = spiketrains[.!isempty.(spiketrains)]

#     dt = 0.005
#     ts = range(epoch, step = dt)
#     tints = [a .. b for (a, b) in zip(ts[1:(end - 1)], ts[2:end])]
#     ss = pmap(spiketrains) do s
#         length.(DimensionalData.groupby(s, 𝑡 => Bins(tints)))
#     end
#     x = sum(ss)
#     x = rectify(set(x, 𝑡 => mean.(times(x))), dims = 𝑡)
# end

# begin
#     f = Figure()
#     ax = Axis(f[1, 1]; yscale = log10)
#     s = spectrum(x)
#     lines!(ax, decompose(s[1:200])...)
#     vlines!.(ax, 1.33 .* (1:7))
#     f
# end
# begin
#     plot(spectrum(x .- mean(x)))
# end
