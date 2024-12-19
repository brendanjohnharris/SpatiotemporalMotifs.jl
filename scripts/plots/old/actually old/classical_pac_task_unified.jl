#! /bin/bash
#=
exec julia -t auto "${BASH_SOURCE[0]}" "$@"
=#
using DrWatson
@quickactivate "SpatiotemporalMotifs"

import TimeseriesTools: freqs
using FileIO
using Unitful
using ModulationIndices
using Normalization
using SpatiotemporalMotifs
@preamble
set_theme!(foresight(:physics))

stimulus = r"Natural_Images"
vars = [:ϕ, :r]

session_table = load(datadir("session_table.jld2"), "session_table")
oursessions = session_table.ecephys_session_id

path = datadir("calculations")
Q = calcquality(path)[Structure = At(structures)]
quality = mean(Q[stimulus = At(stimulus)])

begin
    config = @strdict stimulus vars
    data, file = produce_or_load(produce_uni, config, datadir(); filename = savepath)
    uni = data["uni"]
    unidepths = getindex.(uni, :unidepths)
    layerints = getindex.(uni, :layerints)
    layernames = getindex.(uni, :layernames)
    layernums = getindex.(uni, :layernums)
    r = [uni[i][:r] for i in eachindex(uni)]
    ϕ = [uni[i][:ϕ] for i in eachindex(uni)]
    ϕ = [mod2pi.(x .+ pi) .- pi for x in ϕ]
    uni = []
    GC.gc()

    ϕ_h = [ϕ[i][:, :, lookup(ϕ[i], Trial) .== true] for i in eachindex(ϕ)]
    r_h = [r[i][:, :, lookup(r[i], Trial) .== true] for i in eachindex(r)]
    ϕ_m = [ϕ[i][:, :, lookup(ϕ[i], Trial) .== false] for i in eachindex(ϕ)]
    r_m = [r[i][:, :, lookup(r[i], Trial) .== false] for i in eachindex(r)]
end

begin # * PAC over depths
    pacc = [pac(ϕ[i], r[i]; dims = 𝑡) for i in eachindex(ϕ)]
    pac_h = [pac(ϕ_h[i], r_h[i]; dims = 𝑡) for i in eachindex(ϕ_h)]
    pac_m = [pac(ϕ_m[i], r_m[i]; dims = 𝑡) for i in eachindex(ϕ_m)]
    outs = map(eachindex(pacc)) do i
        s = structures[i]
        μ, (σl, σh) = bootstrapmedian(pacc[i] |> ustripall; dims = 2)
    end
    μs = first.(outs)
    σs = last.(outs)
end

begin # * Set up figure
    f = FourPanel()
    gs = subdivide(f, 2, 2)
end
begin # * Plot over layers
    ax = Axis(gs[1][1, 1], limits = ((0, 1), (nothing, nothing)),
              xlabel = "Cortical depth (%)", ylabel = "PAC", backgroundcolor = :transparent,
              xtickformat = depthticks)
    for i in eachindex(pacc)
        s = structures[i]
        μ = μs[i]
        (σl, σh) = σs[i]
        mu = upsample(μ, 5)
        l = upsample(σl, 5) |> parent
        h = upsample(σh, 5) |> parent
        band!(ax, lookup(mu, 1), l, h; color = (structurecolormap[s], 0.2),
              label = s)
        lines!(ax, lookup(mu, 1), mu, color = structurecolormap[s], label = s,
               alpha = 0.8,
               linewidth = 4)
        scatter!(ax, lookup(μ, 1), μ, color = structurecolormap[s], label = s,
                 markersize = 10, alpha = 0.8)
    end
    axislegend(ax; position = :lt, nbanks = 2, merge = true)

    # ax2 = Axis(f[2, 1], limits = ((0, 100), (nothing, nothing)),
    #            xlabel = "Cortical depth (%)",
    #            ylabel = "PAC difference (hit - miss)")
    # for i in eachindex(pac_m)
    #     lines!(ax2, 100 * lookup(pac_m[i], Depth), ustripall(pac_h[i] .- pac_m[i]),
    #            color = (structurecolors[i], 0.8), label = structures[i])
    # end
    # axislegend(ax2)

    if true # * Add roseflowers
        idx = Depth(Near(0.25))
        struc = 1
        tortinset!(f[1, 1], ϕ[struc][idx], r[struc][idx],
                   colormap = seethrough(structurecolors[struc], 0.1, 1),
                   halign = 0.3, valign = 0.6, color = structurecolors[struc])
        scatter!(ax, [25], [μs[struc][idx]], color = structurecolors[struc],
                 markersize = 15,
                 strokecolor = :white, strokewidth = 2)

        idx = Depth(Near(0.3))
        struc = 2
        tortinset!(f[1, 1], ϕ[struc][idx], r[struc][idx],
                   colormap = seethrough(structurecolors[struc], 0.1, 1),
                   halign = 0.01, valign = 0.4, color = structurecolors[struc])
        scatter!(ax, [30], [μs[struc][idx]], color = structurecolors[struc],
                 markersize = 15,
                 strokecolor = :white, strokewidth = 2)

        idx = Depth(Near(0.9))
        struc = 6
        tortinset!(f[1, 1], ϕ[struc][idx], r[struc][idx],
                   colormap = seethrough(structurecolors[struc], 0.1, 1),
                   halign = 0.87, valign = 0.95, color = structurecolors[struc])
        scatter!(ax, [88], [μs[struc][idx]], color = structurecolors[struc],
                 markersize = 15,
                 strokecolor = :white, strokewidth = 2)

        idx = Depth(Near(0.9))
        struc = 3
        tortinset!(f[1, 1], ϕ[struc][idx], r[struc][idx],
                   colormap = seethrough(structurecolors[struc], 0.1, 1),
                   halign = 0.9, valign = 0.48, color = structurecolors[struc])
        scatter!(ax, [92], [μs[struc][idx]], color = structurecolors[struc],
                 markersize = 15,
                 strokecolor = :white, strokewidth = 2)
    end
end

begin # * Polar plot of coupling angle peaks (find peaks? what if there is more than 1? )
    function phipeak(r, ϕ; n = 50)
        ϕ = ModulationIndices.tortbin(ϕ; n)
        h = [mean(r[ϕ .== i]) for i in 1:n]
        angles = range(start = -π + π / n, stop = π - π / n, length = n) |> collect
        _, i = findmax(h)
        phimax = angles[i]
    end
    function phipeaks(r, ϕ; kwargs...)
        peaks = progressmap(eachslice(r, dims = Depth),
                            eachslice(ϕ, dims = Depth); parallel = true) do a, b
            phipeak(a, b; kwargs...)
        end
        out = ϕ[1, :, 1] # Just a depth slice
        out .= peaks
        return out
    end
    peaks = phipeaks.(r, ϕ; n = 50)
end
begin
    f = Figure()
    ax = PolarAxis(f[1, 1])
    # for (i, p) in enumerate(peaks)
    i = 1
    p = peaks[1]
    lines!(ax, collect(p) .+ pi, collect(lookup(p, 1)); label = structures[i],
           color = structurecolors[i])
    # end
    f
end

# ---------------------------------- Hit/miss prediction --------------------------------- #
# begin # * LDA input data and downsampling
#     H = [getindex.(Og, s) for s in eachindex(Og[1])] # Grouped by subject

#     # * Make sure temporal dims same size
#     ts = intersect([intersect(lookup.(h, 𝑡)...) for h in H]...)
#     H = [[_h[Ti = At(ts)] for _h in h] for h in H]
#     H = stack.([Structure(structures)], H, dims = 3)
#     H = permutedims.(H, [(1, 3, 2)])

#     H = [h[1:10:end, :, :] for h in H]
# end

if true
    begin # * Classify using classical PAC over depths
        regcoef = 0.01
        folds = 5
        repeats = 10
        # bac = progressmap(pacc) do p
        #     classify_kfold(p; regcoef, k = folds, repeats)
        # end
        # p = vcat(pacc...)
        p = pacc[6]
        # p = ToolsArray(p, (Dim{:d}(1:size(p, 1)), dims(pacc[1], 2)))
        bac = classify_kfold(p[:, 1:200]; regcoef = 1, k = folds, repeats)
        x = dropdims(mean(p[:, lookup(p, Trial)]; dims = 2); dims = 2)
        x = (x .- dropdims(mean(p[:, .!lookup(p, Trial)]; dims = 2); dims = 2)) ./ x
        plot(x)
        classifier(p; regcoef = 0.5)
    end

    begin # * Single-subject classifications, returning 5-fold balanced accuracy
        regcoef = 0.5
        folds = 5
        repeats = 10
        bac_pred = pmap(H) do h
            h = h[Ti = -0.25u"s" .. 0.0u"s"]
            bac = classify_kfold(h; regcoef, k = folds, repeats)
        end

        bac_pre = pmap(H) do h
            h = h[Ti = -0.25u"s" .. 0.25u"s"]
            bac = classify_kfold(h; regcoef, k = folds, repeats)
        end

        bac_post = pmap(H) do h
            h = h[Ti = 0.25u"s" .. 0.75u"s"]
            bac = classify_kfold(h; regcoef, k = folds, repeats)
        end

        bac_sur = pmap(H) do h
            h = h[Ti = -0.25u"s" .. 0.25u"s"]
            idxs = randperm(size(h, Trial))
            h = set(h, Trial => lookup(h, Trial)[idxs])
            bac = classify_kfold(h; regcoef, k = folds, repeats)
        end

        D = @strdict bac_pred bac_pre bac_post bac_sur regcoef folds repeats
    end

    begin # * Map of region-wise weightings
        W = map(H) do h
            h = h[Ti = -0.25u"s" .. 0.75u"s"] # !!!
            N, M = classifier(h; regcoef = 0.5) # !!!
            W = projection(M)
            W = reshape(W, size(h)[1:2])
            return ToolsArray(W, dims(h)[1:2])
        end
        W = stack(SessionID(sessionids), W, dims = 3)
        W = W ./ maximum(abs.(W))
    end
    push!(D, "W" => w)
    # tagsave(datafile, D)
end

begin # * Save
    display(f)
    wsave(plotdir("classical_pac_task", "pac_mean.pdf"), f)
end

# begin # * Single-trial PAC
#     function classicalpac(ϕ, r; kwargs...)
#         pac = dropdims(sum(ϕ, dims = (Ti,)); dims = (Ti))
#         pac .= 0.0
#         mp = zip(eachslice(parent(ϕ); dims = (2, 3)), eachslice(parent(r); dims = (2, 3)))
#         Threads.@threads for (i, (_p, _a)) in collect(enumerate(mp))
#             pac[i] = ModulationIndices.tort2010(_p[:], _a[:]; kwargs...)
#         end
#         return pac
#     end

#     pacc = classicalpac.(ϕ, r; n = 20)
#     pac_h = [pacc[i][:, lookup(pac[i], Trial) .== true] for i in eachindex(pacc)]
#     pac_m = [pacc[i][:, lookup(pac[i], Trial) .== false] for i in eachindex(pacc)]
# end
# begin
#     f = Figure(size = (560, 720))
#     ax = Axis(f[1, 1], limits = ((0, 100), (nothing, nothing)),
#               xlabel = "Cortical depth (%)", ylabel = "Mean PAC")
#     for i in eachindex(pac)
#         μ = dropdims(mean(pac[i], dims = Trial), dims = Trial) |> ustripall
#         σ = dropdims(std(pacc[i], dims = Trial), dims = Trial) |> ustripall
#         σ = σ ./ sqrt(size(pacc[i], 2)) # SEM

#         band!(ax, 100 * lookup(pacc[i], Depth), collect(μ .- σ),
#               collect(μ .+ σ);
#               color = (structurecolors[i], 0.3))
#         lines!(ax, 100 * lookup(pacc[i], Depth), μ,
#                color = (structurecolors[i], 0.8), label = structures[i])
#     end
#     axislegend(ax; position = :ct, nbanks = 2)

#     ax = Axis(f[2, 1], limits = ((0, 100), (nothing, nothing)),
#               xlabel = "Cortical depth (%)",
#               ylabel = "PAC difference (Cohen's d)")
#     for i in eachindex(pacc)
#         μ_h = dropdims(mean(pac_h[i], dims = Trial), dims = Trial) |> ustripall
#         μ_m = dropdims(mean(pac_m[i], dims = Trial), dims = Trial) |> ustripall
#         σ = dropdims(std(pacc[i], dims = Trial), dims = Trial) |> ustripall
#         d = (μ_h .- μ_m) ./ σ # * Cohen's d effect size

#         lines!(ax, 100 * lookup(d, Depth), d,
#                color = (structurecolors[i], 0.8), label = structures[i])
#     end
#     axislegend(ax; position = :rt, nbanks = 2)

#     display(f)
#     wsave(plotdir("classical_pac_task", "pac_singletrial.pdf"), f)
# end
