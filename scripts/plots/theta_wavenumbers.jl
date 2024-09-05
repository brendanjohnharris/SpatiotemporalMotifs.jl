#! /bin/bash
#=
exec julia -t auto "${BASH_SOURCE[0]}" "$@"
=#
# ? Expected execution time: 30 mins
using DrWatson
@quickactivate "SpatiotemporalMotifs"

import TimeseriesTools: freqs
using FileIO
using Unitful
import DimensionalData: metadata
using MultivariateStats
using SpatiotemporalMotifs
@preamble
set_theme!(foresight(:physics))
Random.seed!(32)

vars = [:k]
INTERVAL = SpatiotemporalMotifs.INTERVAL
mainstructure = "VISl"
maincolorrange = [-2.5, 2.5]
session_table = load(datadir("posthoc_session_table.jld2"), "session_table")
oursessions = session_table.ecephys_session_id
path = datadir("calculations")

begin # * Set up main figure
    mf = TwoPanel()
    mfs = subdivide(mf, 1, 4)
end

stimulus = r"Natural_Images"
uni = load_uni(; stimulus, vars)

begin # * Supplemental material: average wavenumbers in each region
    f = Figure(size = (720, 1440))

    for i in eachindex(uni)
        k = uni[i][:k][:, 2:end, :]
        k = uconvert.(u"mm^-1", k)

        # * Hit
        ax = Axis(f[i, 1], yreversed = true)
        ax.limits = (nothing, (0.1, 0.95))
        structure = metadata(k)[:structure]
        ax.title = structure * ": hit"
        m = median(k[:, :, lookup(k, Trial) .== true], dims = Trial)
        m = dropdims(m, dims = Trial)[𝑡(SpatiotemporalMotifs.INTERVAL)]
        colorrange = maximum(abs.(ustripall(m))) * [-1, 1]
        ints = uni[i][:layerints]
        p = plotlayermap!(ax, m, ints; arrows = true, colorrange) |> first
        if i > 5
            ax.xlabel = "Time (s)"
        end

        if structure == mainstructure # * Plot into main figure
            ax = Axis(mfs[1], yreversed = true)
            ax.xlabel = "Time (s)"
            ax.limits = (nothing, (0.1, 0.95))
            ax.title = structure * ": hit"
            m = median(k[:, :, lookup(k, Trial) .== true], dims = Trial)
            m = dropdims(m, dims = Trial)[𝑡(SpatiotemporalMotifs.INTERVAL)]
            ints = uni[i][:layerints]
            p = plotlayermap!(ax, m, ints; arrows = true, colorrange) |> first
            c = Colorbar(mfs[end]; colorrange = maincolorrange, colormap = defaultcolormap,
                         highclip = defaultcolormap[end], lowclip = defaultcolormap[1])
            c.label = "θ wavenumber ($(unit(eltype(k))))"
            # push!(maincolorrange, colorrange)
        end

        # * Miss
        ax = Axis(f[i, 2], yreversed = true)
        ax.limits = (nothing, (0.1, 0.95))
        ax.title = structure * ": miss"
        m = median(k[:, :, lookup(k, Trial) .== false], dims = Trial)
        m = dropdims(m, dims = Trial)[𝑡(SpatiotemporalMotifs.INTERVAL)]
        ints = uni[i][:layerints]
        p = plotlayermap!(ax, m, ints; arrows = true, colorrange) |> first
        c = Colorbar(f[i, 3], p)
        c.label = "θ wavenumber ($(unit(eltype(k))))"

        if structure == mainstructure # * Plot into main figure
            ax = Axis(mfs[2], yreversed = true)
            ax.xlabel = "Time (s)"
            ax.limits = (nothing, (0.1, 0.95))
            ax.title = structure * ": miss"
            m = median(k[:, :, lookup(k, Trial) .== false], dims = Trial)
            m = dropdims(m, dims = Trial)[𝑡(SpatiotemporalMotifs.INTERVAL)]
            ints = uni[i][:layerints]
            p = plotlayermap!(ax, m, ints; arrows = true, colorrange = maincolorrange) |>
                first
        end
    end
    addlabels!(f)
    display(f)
    wsave(plotdir("theta_wavenumbers", "supplemental_wavenumber.pdf"), f)
end

# ? Flashes
stimulus = "flash_250ms"
uni = load_uni(; stimulus, vars)

begin # * Supplemental material: average wavenumbers in each region
    f = SixPanel()
    gs = subdivide(f, 3, 2)

    for i in eachindex(uni)
        k = uni[i][:k][:, 2:end, :, :]
        k = uconvert.(u"mm^-1", k)
        structure = metadata(k)[:structure]

        ax = Axis(gs[i][1, 1], yreversed = true, xlabel = "Time (s)")
        ax.limits = (nothing, (0.1, 0.95))
        ax.title = structure
        m = median(k, dims = (Trial, SessionID))
        m = dropdims(m, dims = (Trial, SessionID))[𝑡(SpatiotemporalMotifs.INTERVAL)]
        ints = uni[i][:layerints]
        p = plotlayermap!(ax, m, ints; arrows = true) |> first
        c = Colorbar(gs[i][1, 2], p)
        c.label = "θ wavenumber ($(unit(eltype(k))))"
        if structure == mainstructure
            colorrange = maincolorrange
            ax = Axis(mfs[3], yreversed = true, xlabel = "Time (s)")
            ax.limits = (nothing, (0.1, 0.95))
            ax.title = structure * ": flashes"
            m = median(k, dims = (Trial, SessionID))
            m = dropdims(m, dims = (Trial, SessionID))[𝑡(SpatiotemporalMotifs.INTERVAL)]
            ints = uni[i][:layerints]
            p = plotlayermap!(ax, m, ints; arrows = true, colorrange) |> first
        end
    end
    addlabels!(f)
    display(f)
    wsave(plotdir("theta_wavenumbers", "supplemental_wavenumber_flash.pdf"), f)
end

begin # * Save main figure
    addlabels!(mf)
    wsave(plotdir("theta_wavenumbers", "theta_wavenumber.pdf"), mf)
    display(mf)
end
# begin # * Set up main figure
#     fig = TwoPanel()
#     gs = subdivide(fig, 1, 2)
# end

# begin # * Calculate a global order parameter at each time point
#     Og = map(out) do o
#         O = map(o) do p
#             k = p[:k]
#             ret = dropdims(mean(sign.(k), dims = Depth), dims = Depth)
#             trials = p[:trials][1:size(k, :changetime), :]
#             trialtimes = trials.change_time_with_display_delay
#             @assert all(isapprox.(ustripall(lookup(k, :changetime)), trialtimes,
#                                   atol = 0.05)) # These won't match exactly, because the data change times have been adjusted for rectification. This is OK.
#             set(ret, Dim{:changetime} => Trial(trials.hit))
#         end
#     end

#     Og_h = [[o[:, lookup(o, Trial) .== true] for o in O] for O in Og]
#     Og_m = [[o[:, lookup(o, Trial) .== false] for o in O] for O in Og]
# end

# begin # * Plot the mean order parameter across time, and the hit/miss contrast
#     function mf(Og)
#         map(Og) do O
#             o = dropdims.(mean.(O; dims = Trial), dims = Trial)
#             sessionids = [metadata(_o)[:sessionid] for _o in O]
#             o = [_o[(end - minimum(length.(o)) + 1):end] for _o in o]
#             stack(SessionID(sessionids), o)
#         end
#     end
#     Ō = mf(Og)
#     Ō_h = mf(Og_h)
#     Ō_m = mf(Og_m)

#     ax = Axis(gs[1]; xlabel = "Time (s)", ylabel = "Mean order parameter",
#               xautolimitmargin = (0, 0), xminorticksvisible = true,
#               xminorticks = IntervalsBetween(5), yminorticksvisible = true,
#               yminorticks = IntervalsBetween(5))
#     hlines!(ax, [0]; color = (:black, 0.5), linestyle = :dash, linewidth = 2)
#     vlines!(ax, [0, 0.25]; color = (:black, 0.5), linestyle = :dash, linewidth = 2)
#     for (i, O) in reverse(collect(enumerate(Ō)))
#         structure = metadata(O)[:structure]
#         O = O[𝑡(SpatiotemporalMotifs.INTERVAL)]
#         μ = dropdims(mean(O, dims = SessionID), dims = SessionID)
#         σ = dropdims(std(O, dims = SessionID), dims = SessionID)
#         σ = σ ./ 2
#         bargs = [times(μ), μ .- σ, μ .+ σ] .|> ustripall .|> collect
#         band!(ax, bargs..., color = (structurecolors[i], 0.3), label = structure)
#         lines!(ax, times(μ) |> ustripall, μ |> ustripall |> collect,
#                color = (structurecolors[i], 0.7),
#                label = structure)
#     end
#     l = axislegend(ax, position = :lt, nbanks = 2, framevisible = true, labelsize = 12,
#                    merge = true)
#     reverselegend!(l)

#     # f = Figure()
#     # ax = Axis(f[1, 1]; xlabel = "Time (s)", ylabel = "Mean order parameter (hit - miss)",
#     #           xautolimitmargin = (0, 0), xminorticksvisible = true,
#     #           xminorticks = IntervalsBetween(5), yminorticksvisible = true,
#     #           yminorticks = IntervalsBetween(5))
#     # hlines!(ax, [0]; color = (:black, 0.5), linestyle = :dash, linewidth = 2)
#     # vlines!(ax, [0, 0.25]; color = (:black, 0.5), linestyle = :dash, linewidth = 2)
#     # for (i, (O_h, O_m)) in enumerate(zip(Ō_h, Ō_m))
#     #     structure = metadata(O_h)[:structure]
#     #     x = O_h .- O_m
#     #     μ = dropdims(mean(x, dims = SessionID), dims = SessionID)
#     #     σ = dropdims(std(x, dims = SessionID), dims = SessionID)
#     #     σ = σ ./ sqrt(size(O_h, 2)) # SEM
#     #     bargs = [times(μ), μ .- σ, μ .+ σ] .|> ustripall .|> collect
#     #     band!(ax, bargs..., color = (structurecolors[i], 0.3))
#     #     lines!(ax, μ |> ustripall, color = (structurecolors[i], 0.7), label = structure)
#     # end
#     # axislegend(ax, position = :rb, framevisible = true, labelsize = 12)
#     # display(f)
#     # wsave(plotdir("theta_waves_task", "regional_orderparameter.pdf"), f)
# end
