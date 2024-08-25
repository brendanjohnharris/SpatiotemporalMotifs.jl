#! /bin/bash
#=
exec julia -t auto "${BASH_SOURCE[0]}" "$@"
=#
using DrWatson
@quickactivate "SpatiotemporalMotifs"

import TimeseriesTools: freqs
using FileIO
using Unitful
using Normalization
using SpatiotemporalMotifs
@preamble
set_theme!(foresight(:physics))

stimulus = r"Natural_Images"
vars = [:ϕ, :r]

session_table = load(datadir("posthoc_session_table.jld2"), "session_table")
oursessions = session_table.ecephys_session_id

path = datadir("calculations")
Q = calcquality(path)[structure = At(structures)]
Q = Q[sessionid = At(oursessions)]
quality = mean(Q[stimulus = At(stimulus)])

begin # * Set up main figure
    mf = TwoPanel()
    mgs = subdivide(mf, 1, 2)
end

stimuli = ["r\"Natural_Images\"", "spontaneous", "flash_250ms"]
pQ = calcquality(datadir("power_spectra"))
for stimulus in stimuli
    _Q = pQ[stimulus = At(stimulus), structure = At(structures)]
    subsessions = intersect(oursessions, lookup(_Q, :sessionid))
    if length(subsessions) < length(oursessions)
        @warn "Power spectra calculations are incomplete, proceeding regardless"
    end
    _Q = _Q[Dim{:sessionid}(At(subsessions))]
    filebase = stimulus == "spontaneous" ? "" : "_$stimulus"
    f = Figure(size = (900, 1080))

    begin # * Load data
        S = map(lookup(_Q, :structure)) do structure
            out = map(lookup(_Q, :sessionid)) do sessionid
                if _Q[sessionid = At(sessionid), structure = At(structure)] == 0
                    return nothing
                end
                filename = savepath((@strdict sessionid structure stimulus), "jld2",
                                    datadir("power_spectra"))
                C = load(filename, "C")
                # S = load(filename, "sC")
                # return (C .- S) ./ median(S)
            end
            idxs = .!isnothing.(out)
            out = out[idxs]

            m = DimensionalData.metadata.(out)
            out = map(out) do o
                dropdims(mean(o, dims = Chan), dims = Chan)
            end
            out = stack(Dim{:sessionid}(lookup(_Q, :sessionid)[idxs]), out, dims = 3)
            return out
        end
    end

    begin # * Supplemental average comodulograms
        f = SixPanel()
        gs = subdivide(f, 3, 2)
        map(gs, structures, S) do g, structure, s
            ax = Axis(g[1, 1]; title = structure, xlabel = "Phase frequency (Hz)",
                      ylabel = "Amplitude frequency (Hz)")
            s = dropdims(mean(s, dims = :sessionid); dims = :sessionid)
            s = upsample(s, 5, 1)
            s = upsample(s, 5, 2)
            p = heatmap!(ax, s; colormap = seethrough(reverse(sunrise)), rasterize = 5)
            Colorbar(g[1, 2], p; label = "Modulation index")
        end
        addlabels!(f)
        display(f)
        wsave(plotdir("comodulograms", "comodulograms_$stimulus.pdf"), f)
    end

    if stimulus == "spontaneous" # * Plot into main figure
        structure = "VISl"
        ax = Axis(mgs[1][1, 1]; title = structure, xlabel = "Phase frequency (Hz)",
                  ylabel = "Amplitude frequency (Hz)")
        s = S[lookup(_Q, :structure) .== structure] |> only
        display(s)
        s = dropdims(mean(s, dims = :sessionid); dims = :sessionid)
        s = upsample(s, 5, 1)
        s = upsample(s, 5, 2)
        p = heatmap!(ax, s; colormap = seethrough(reverse(sunrise)), rasterize = 5)
        Colorbar(mgs[1][1, 2], p; label = "Mean PAC")
    end
end

config = @strdict stimulus vars
data, file = produce_or_load(produce_uni, config, datadir(); filename = savepath)
uni = data["uni"]

unidepths = getindex.(uni, :unidepths)
layerints = getindex.(uni, :layerints)
layernames = getindex.(uni, :layernames)
layernums = getindex.(uni, :layernums)

begin # * Normalize amplitudes and generate a burst mask
    r = [abs.(uni[i][:r]) for i in eachindex(uni)]
    r_h = [r[i][:, :, lookup(r[i], :trial) .== true] for i in eachindex(r)]
    r_m = [r[i][:, :, lookup(r[i], :trial) .== false] for i in eachindex(r)]

    ϕ = [uni[i][:ϕ] for i in eachindex(uni)]
    ϕ = [mod2pi.(x .+ pi) .- pi for x in ϕ]
    ϕ_h = [ϕ[i][:, :, lookup(ϕ[i], :trial) .== true] for i in eachindex(ϕ)]
    ϕ_m = [ϕ[i][:, :, lookup(ϕ[i], :trial) .== false] for i in eachindex(ϕ)]

    uni = []
    GC.gc()
end

begin # * Pac
    PAC = progressmap(ϕ, r) do ϕ, r
        pac(ϕ, r; dims = :trial)
    end
    PAC_h = progressmap(ϕ_h, r_h) do ϕ, r
        pac(ϕ, r; dims = :trial)
    end
    PAC_m = progressmap(ϕ_m, r_m) do ϕ, r
        pac(ϕ, r; dims = :trial)
    end
end

begin # * Supplemental figure: spatiotemporal PAC over all regions
    cmax = maximum(PAC[1])
    f = SixPanel()
    gs = subdivide(f, 3, 2)
    for (g, l, P) in zip(gs, layerints, PAC)
        s = metadata(P)[:structure]
        ax = Axis(g[1, 1]; title = s, yreversed = true,
                  limits = (nothing, (extrema(lookup(P, :depth)))))
        p = plotlayermap!(ax, ustripall(P[Ti = SpatiotemporalMotifs.INTERVAL]), l) |> first
        Colorbar(g[1, 2], p, label = "PAC")

        if s == "VISl"
            ax = Axis(mgs[2][1, 1]; title = s, yreversed = true,
                      limits = (nothing, (extrema(lookup(P, :depth)))))
            p = plotlayermap!(ax, ustripall(P[Ti = SpatiotemporalMotifs.INTERVAL]), l) |>
                first
            Colorbar(mgs[2][1, 2], p, label = "PAC")
        end
    end
    addlabels!(f)
    f
end

begin # * Save figure
    addlabels!(mf)
    wsave(plotdir("pac_task", "average_pac.pdf"), mf)
end

# begin
#     mask_h = [r̂_h[i] .> 2.0 for i in eachindex(r̂_h)]
#     mask_m = [r̂_m[i] .> 2.0 for i in eachindex(r̂_m)]

#     ϕ = [abs.(uni[i][:ϕ]) for i in eachindex(uni)]

#     ϕ_h = [ϕ[i][:, :, lookup(ϕ[i], :trial) .== true] for i in eachindex(ϕ)]
#     ϕ̂_h = deepcopy(ϕ_h)
#     [ϕ̂_h[i][.!mask_h[i]] .= NaN for i in eachindex(mask_h)]

#     ϕ_m = [ϕ[i][:, :, lookup(ϕ[i], :trial) .== false] for i in eachindex(ϕ)]
#     ϕ̂_m = deepcopy(ϕ_m)
#     [ϕ̂_m[i][.!mask_m[i]] .= NaN for i in eachindex(mask_m)]

#     ϕ = []
#     uni = []
#     GC.gc()
# end

# begin # * Heatmap of the number of samples
#     pn_h = dropdims.(sum.(mask_h, dims = 3), dims = 3)
#     pn_m = dropdims.(sum.(mask_m, dims = 3), dims = 3)

#     f = Figure(size = (720, 1440))

#     for i in eachindex(pn_h)
#         # * Hit
#         ax = Axis(f[i, 1], yreversed = true)
#         structure = DimensionalData.metadata(pn_h[i])[:structure]
#         ax.title = structure * ": hit"
#         m = ustripall(pn_h[i])
#         colorrange = (nmin, max(maximum(m), maximum(ustripall(pn_m[i]))))
#         ints = layerints[i]
#         p = plotlayermap!(ax, m, ints; arrows = false, colorrange,
#                           colormap = :bone, lowclip = :crimson) |>
#             first
#         if i > 5
#             ax.xlabel = "Time (s)"
#         end

#         # * Miss
#         ax = Axis(f[i, 2], yreversed = true)
#         ax.title = structure * ": miss"
#         m = ustripall(pn_m[i])
#         p = plotlayermap!(ax, m, ints; arrows = false, colorrange,
#                           colormap = :bone, lowclip = :crimson) |>
#             first
#         c = Colorbar(f[i, 3], p)
#         c.label = "Num. of bursts"
#     end

#     wsave(plotdir("pac_task", "supplemental_numsamples.pdf"), f)
# end

# begin # * Heatmap of the PAC entropy
#     function pacentropy(ϕ) # ?! Normalised or unnormalized? Domain is always the same, so...?
#         _p = dropdims(nansafe(x -> StateSpaceSet(x); dims = 3)(ϕ), dims = 3)
#         function ent(x)
#             if isempty(x)
#                 return 3 # Some entropy higher than the maximum possible, 1
#             else
#                 ComplexityMeasures.information_normalized(Shannon(),
#                                                           ValueBinning(RectangularBinning(20)),
#                                                           StateSpaceSet(x))
#             end
#         end
#         return ent.(_p)
#     end

#     PAC_h = pacentropy.(ϕ̂_h)
#     PAC_m = pacentropy.(ϕ̂_m)

#     f = Figure(size = (720, 1440))

#     for i in eachindex(PAC_h)
#         # * Hit
#         ax = Axis(f[i, 1], yreversed = true)
#         structure = DimensionalData.metadata(pn_h[i])[:structure]
#         ax.title = structure * ": hit"
#         m = ustripall(PAC_h[i])
#         colorrange = [extrema(m)..., extrema(ustripall(PAC_m[i]))...] |> extrema
#         m[pn_h[i] .< nmin] .= 2
#         ints = layerints[i]
#         p = plotlayermap!(ax, m, ints; arrows = false, colormap = :bone, colorrange,
#                           highclip = :crimson) |>
#             first
#         if i > 5
#             ax.xlabel = "Time (s)"
#         end

#         # * Miss
#         ax = Axis(f[i, 2], yreversed = true)
#         ax.title = structure * ": miss"
#         m = ustripall(PAC_m[i])
#         m[pn_m[i] .< nmin] .= 2
#         p = plotlayermap!(ax, m, ints; arrows = false, colorrange, highclip = :crimson,
#                           colormap = :bone) |>
#             first
#         c = Colorbar(f[i, 3], p)
#         c.label = "PAC entropy"
#     end
#     wsave(plotdir("pac_task", "supplemental_PAC.pdf"), f)
# end
