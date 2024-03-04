#! /bin/bash
#=
exec julia -t auto "${BASH_SOURCE[0]}" "$@"
=#
using DrWatson
@quickactivate "SpatiotemporalMotifs"

import SpatiotemporalMotifs as SM
import TimeseriesTools: freqs
using Peaks
using FileIO
SM.@preamble
set_theme!(foresight(:physics))

stimulus = "spontaneous"
xtickformat = x -> string.(round.(Int, x))
theta = 4 .. 10
gamma = 40 .. 100

session_table = load(datadir("session_table.jld2"), "session_table")
oursessions = session_table.ecephys_session_id

path = datadir("PowerSpectra")
Q = SM.calcquality(path)[stimulus = At(stimulus)]

begin # * Load data
    S = map(lookup(Q, :structure)) do structure
        out = map(lookup(Q, :sessionid)) do sessionid
            if Q[sessionid = At(sessionid), structure = At(structure)] == 0
                return nothing
            end
            filename = joinpath(path,
                                savename((@strdict sessionid structure stimulus), "jld2",
                                         connector = SM.connector))
            S = load(filename, "S")
        end
        out = filter(!isnothing, out)
        out = filter(x -> maximum(DimensionalData.metadata(x)[:streamlinedepths]) > 0.90,
                     out) # Remove sessions that don't have data to a reasonable depth

        m = DimensionalData.metadata.(out)
        sessions = getindex.(m, :sessionid)

        streamlinedepths = getindex.(m, :streamlinedepths)
        layerinfo = getindex.(m, :layerinfo)

        unidepths = SM.commondepths(streamlinedepths)
        out = map(out, streamlinedepths, layerinfo) do o, s, l
            o = set(o, Dim{:depth} => s)
            layernames = DimArray(l[1], Dim{:depth}(lookup(o, Dim{:depth})))
            layernums = DimArray(l[3], Dim{:depth}(lookup(o, Dim{:depth})))
            o = o[Dim{:depth}(Near(unidepths))]
            layernames = layernames[Dim{:depth}(Near(unidepths))]
            layernums = layernums[Dim{:depth}(Near(unidepths))]
            # @assert length(unique(lookup(o, Dim{:depth}))) == length(unidepths)
            @assert issorted(lookup(o, Dim{:depth}))
            push!(o.metadata, :layernames => layernames)
            push!(o.metadata, :layernums => layernums)
            o = set(o, Dim{:depth} => unidepths)
        end
        layernames = DimArray(stack(getindex.(DimensionalData.metadata.(out), :layernames)),
                              (Dim{:depths}(unidepths), Dim{:sessionid}(sessions)))
        layernums = DimArray(stack(getindex.(DimensionalData.metadata.(out), :layernums)),
                             (Dim{:depths}(unidepths), Dim{:sessionid}(sessions)))
        S = stack(Dim{:sessionid}(sessions), out, dims = 3)
        layernums = SM.parselayernum.(layernames)
        return S, layernames, layernums
    end
    S, layernames, layernums = zip(S...)
end

begin # * Format layers
    meanlayers = map(layernums) do l
        round.(Int, mean(l, dims = 2))
    end
    S̄ = map(S, meanlayers) do s, l
        s = set(s, Dim{:depth} => Dim{:layer}(SM.layernum2name.(parent(l)[:])))
        s = set(s, :layer => DimensionalData.Irregular)
    end
    S̄ = map(S̄) do s
        ss = map(unique(lookup(s, :layer))) do l
            ls = s[Dim{:layer}(At(l))]
            if hasdim(ls, :layer)
                ls = mean(ls, dims = :layer)
            end
            ls
        end
        cat(ss..., dims = :layer)
    end
    S̄ = DimArray(S̄ |> collect, (Dim{:structure}(lookup(Q, :structure)),))
    S = DimArray(S |> collect, (Dim{:structure}(lookup(Q, :structure)),))
    S̄ = map(S̄) do s
        N = UnitEnergy(s, dims = 1)
        N(s)
    end
end;

begin # * Mean power spectrum in V1
    f = Figure()
    structure = "VISp"
    ax = Axis(f[1, 1]; xscale = log10, yscale = log10,
              limits = ((3, 300), (nothing, nothing)),
              xlabel = "Frequency (Hz)",
              xgridvisible = true,
              ygridvisible = true,
              xgridstyle = :dash,
              ygridstyle = :dash,
              xtickformat)

    # * Band annotations
    vspan!(ax, extrema(theta)..., color = (crimson, 0.22),
           label = "𝛉 ($(theta.left) – $(theta.right) Hz)")
    vspan!(ax, extrema(gamma)..., color = (cornflowerblue, 0.22),
           label = "𝛄 ($(gamma.left) – $(gamma.right) Hz)")

    s = S̄[structure = At(structure)]
    μ = mean(s, dims = (:sessionid, :layer))
    μ = dropdims(μ, dims = (:sessionid, :layer)) |> ustripall
    σ = std(s, dims = (:sessionid, :layer)) ./ 2
    σ = dropdims(σ, dims = (:sessionid, :layer)) |> ustripall
    lines!(ax, freqs(μ), μ; color = cucumber)
    band!(ax, freqs(μ), collect.([μ - σ, μ + σ])...; color = (cucumber, 0.32))

    # * Find peaks
    pks, proms = findpeaks(μ, 2; N = 2)
    scatter!(ax, freqs(pks), pks .* 1.25, color = :black,
             markersize = 10, marker = :dtriangle)
    text!(ax, freqs(pks), pks;
          text = string.(round.(freqs(pks), digits = 1)) .* [" Hz"],
          align = (:center, :bottom), color = :black, rotation = 0,
          fontsize = 16,
          offset = (0, 15))

    # * Fooof fit
    fooof = x -> AN.aperiodicfit(x, [3, 300]; aperiodic_mode = "fixed", max_n_peaks = 5,
                                 peak_threshold = 1, peak_width_limits = [1, 50])
    ff, ps = fooof(μ)
    lines!(ax, freqs(μ), ff.(freqs(μ)), color = (:gray, 0.8), linestyle = :dash)
    text!(ax, 14, exp10(-2.5); text = "𝛂 = $(round(ps[:χ], sigdigits=3))", fontsize = 16)

    axislegend(ax, position = :rt, labelsize = 13)
    f
end
begin # * Calculate the channel-wise fits
    Sl = map(S) do s
        s = dropdims(mean(s, dims = :sessionid), dims = :sessionid)
    end
    L = map(Sl) do s
        fooof.(eachslice(ustripall(s), dims = (2,)))
    end
    χ = [getindex.(last.(l), :χ) for l in L]
    L = [first.(l) for l in L]
end
begin # * Plot the exponent
    f = Figure()
    ax = Axis(f[1, 1]; xlabel = "Cortical depth (%)", ylabel = "1/f exponent",
              limits = ((0, 1), (nothing, nothing)))
    for (i, l) in enumerate(upsample.(χ, 10))
        lines!(ax, lookup(l, 1), l; color = (SM.structurecolors[i], 1),
               label = SM.structures[i])
    end
    axislegend(ax, position = :lt, nbanks = 2)
    display(f)
end

begin # * Plot fooof residuals
    f = Figure()
    ax2 = Axis(f[1, 1]; xscale = log10,
               limits = ((3, 300), (-0.2, 0.9)), xtickformat) # xticksvisible = false, yaxisposition = :right,
    #    xticklabelsvisible = false,
    Sr_log = deepcopy(ustripall.(Sl))
    Sr_log = map(Sr_log, meanlayers) do s, m
        s = set(s, Dim{:depth} => Dim{:layer}(SM.layernum2name.(parent(m)[:])))
        s = set(s, :layer => DimensionalData.Irregular)
    end
    Sr_log = DimArray(Sr_log |> collect,
                      (Dim{:structure}(lookup(Q, :structure)),))
    map(Sr_log, L) do s, l
        for i in eachindex(l)
            s[:, i] .= log10.(s[:, i]) .- log10.(l[i].(freqs(s[:, i])))
        end
    end
    for (c, l) in reverse(collect(zip(SM.layercolors, SM.layers)))
        s = Sr_log[structure = At(structure)][layer = At(l)][Freq(3 .. 300)]
        display(s)
        if ndims(s) > 1
            μ = mean(s, dims = :layer)
            μ = dropdims(μ, dims = :layer) |> ustripall
        else
            μ = s
        end
        # σ = std(s, dims = :layer) ./ 2
        # σ = dropdims(σ, dims = :layer) |> ustripall
        lines!(ax2, TimeseriesTools.freqs(μ), μ; color = (c, 0.7))
        # band!(ax, TimeseriesTools.freqs(μ), collect(μ .- σ), collect(μ .+ σ);
        #       color = (c, 0.32))
    end
    C = Colorbar(f[1, 2], colormap = cgrad(SM.layercolors, categorical = true),
                 ticks = (range(0, 1, length = (2 * length(SM.layercolors) + 1))[2:2:end],
                          ["VISp"] .* SM.layers))
    # linkxaxes!(ax, ax2)
    f
end

begin # * Calculate the residual power in each band
    Sr = deepcopy(ustripall.(Sl))
    map(Sr, L) do s, l
        for i in eachindex(l)
            s[:, i] .= (s[:, i]) .- (l[i].(freqs(s[:, i])))
        end
    end
end

begin # * Plot the total residual theta power across channels
    f = Figure()
    ax = Axis(f[1, 1]; xlabel = "Cortical depth (%)", yticks = WilkinsonTicks(4),
              ylabel = "Residual 𝜽 ($(theta.left) – $(theta.right) Hz) [$(unit(eltype(S[1][1])))]",
              yticklabelrotation = π / 2)

    for (i, s) in enumerate(SM.structures)
        ss = Sr[structure = At(s)][Freq(theta)]
        ss = upsample(ss, 5)
        lines!(ax, lookup(ss, 2), sum(parent(ss), dims = 1)[:] ./ step(lookup(ss, 1));
               color = SM.structurecolors[i], label = SM.structures[i])
    end

    axislegend(ax, position = :lt, nbanks = 2)
    display(f)
end
begin # * Residual gamma power across channels
    f = Figure()
    ax = Axis(f[1, 1]; xlabel = "Cortical depth (%)", yticks = WilkinsonTicks(4),
              ylabel = "Residual 𝜸 ($(gamma.left) – $(gamma.right) Hz) [$(unit(eltype(S[1][1])))]",
              yticklabelrotation = π / 2)

    for (i, s) in enumerate(SM.structures)
        ss = Sr[structure = At(s)][Freq(gamma)]
        ss = upsample(ss, 5)
        lines!(ax, lookup(ss, 2), sum(parent(ss), dims = 1)[:] ./ step(lookup(ss, 1));
               color = SM.structurecolors[i], label = SM.structures[i])
    end

    axislegend(ax, position = :rt, nbanks = 2)
    display(f)
end

begin # * Plot the spectral width of the gamma band
    f = Figure()
    ax = Axis(f[1, 1]; xlabel = "Cortical depth (%)", yticks = WilkinsonTicks(4),
              ylabel = "𝜸 spectral width ($(gamma.left) – $(gamma.right) Hz) [Hz]",
              yticklabelrotation = π / 2)

    for (i, s) in enumerate(SM.structures)
        ss = Sr[structure = At(s)][Freq(gamma)]
        df = ustrip(step(lookup(Sr[1], Freq)))
        N = parent(ss) ./ (sum(parent(ss), dims = 1) ./ df) # A density
        fs = collect(lookup(ss, 1))
        μ = sum(fs .* N, dims = 1)[:] ./ df
        σ = sqrt.(sum((fs .- μ') .^ 2 .* N, dims = 1)[:] ./ df)
        σ = DimArray(σ, (Dim{:depth}(lookup(ss, 2)),))
        σ = upsample(σ, 5)
        lines!(ax, lookup(σ, 1), σ[:];
               color = SM.structurecolors[i], label = SM.structures[i])
    end

    axislegend(ax, position = :lb, nbanks = 2)
    display(f)
end
