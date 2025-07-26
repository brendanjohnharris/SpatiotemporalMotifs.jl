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
import SpatiotemporalMotifs.layers
using Random
@preamble
set_theme!(foresight(:physics))

begin # * Load the CSD flashes
end
# begin # * Load the CSD passive natural images
#     session_table = load(calcdir("posthoc_session_table.jld2"), "session_table")
#     oursessions = session_table.ecephys_session_id
#     path = calcdir("calculations")
#     Q = calcquality(path)[Structure = At(structures)]
#     Q = Q[SessionID(At(oursessions))]
#     @assert mean(Q[stimulus = At("Natural_Images_passive_nochange")]) == 1
#     out = load_calculations(Q; stimulus = "Natural_Images_passive_nochange",
#                                      vars = [:csd])
# end

# begin # * Save average CSD for each session
#     csd = map(out) do out_structure
#         csd = map(out_structure) do out_session
#             csd_session = mapslices(mean, out_session[:csd], dims = :changetime)
#             csd_session = dropdims(csd_session, dims = :changetime)
#             return csd_session
#         end
#         return ToolsArray(csd, (SessionID(oursessions),))
#     end
#     csd = ToolsArray(csd, (Structure(structures),))
#     D = @strdict csd
#     tagsave(calcdir("average_csd.jld2"), D)
# end

begin # * Save average CSD for each session
    conf = Dict("stimulus" => "Natural_Images_passive_nochange")

    data, data_file = produce_or_load(copy(conf), calcdir("plots");
                                      filename = savepath("csd")) do conf
        stimulus = conf["stimulus"]

        session_table = load(calcdir("posthoc_session_table.jld2"), "session_table")
        oursessions = session_table.ecephys_session_id
        path = calcdir("calculations")
        Q = calcquality(path)[Structure = At(structures)]
        Q = Q[SessionID(At(oursessions))]
        @assert mean(Q[stimulus = At(stimulus)]) == 1

        out = load_calculations(Q; stimulus, vars = [:csd])

        csd = map(out) do out_structure
            csd = map(out_structure) do out_session
                csd_session = mapslices(mean, out_session[:csd], dims = :changetime)
                csd_session = dropdims(csd_session, dims = :changetime)
                return csd_session
            end
        end
        csd = map(csd) do x
            filter(x) do _x
                metadata(_x)[:sessionid] in oursessions
            end
        end
        csd = map(csd) do x
            sessionids = map(x -> metadata(x)[:sessionid], x)
            x = ToolsArray(x, (SessionID(sessionids),))
        end

        # session_table = load(calcdir("posthoc_session_table.jld2"), "session_table")
        # oursessions = session_table.ecephys_session_id
        # path = calcdir("calculations")
        # Q = calcquality(path)[Structure = At(structures)]
        # Q = Q[SessionID(At(oursessions))]
        # @assert mean(Q[stimulus = At(r"Natural_Images")]) == 1
        # out = load_calculations(Q; stimulus = r"Natural_Images", vars = [])
        streamlinedepths = map(out) do O
            depths = map(O) do o
                o[:streamlinedepths]
            end
            depths = ToolsArray(depths, (SessionID(oursessions),))
        end
        layernames = map(out) do O
            names = map(O) do o
                o[:layernames]
            end
            names = ToolsArray(names, (SessionID(oursessions),))
        end
        return Dict("csd" => csd,
                    "streamlinedepths" => streamlinedepths,
                    "layernames" => layernames)
    end
    csd = data["csd"]
    streamlinedepths = data["streamlinedepths"]
    layernames = data["layernames"]

    anatomical_l4s = map(layernames) do O
        l4s = map(O) do o
            l4 = contains.(o, ["4"])
            l4 = lookup(o[l4], 1) |> mean # 'center' of l4
        end
        # l4s = ToolsArray(l4s, (SessionID(oursessions),))
    end
end

begin
    function find_L4_center(csd::MultivariateTimeSeries, doplot = false)

        # * Go back to probe depths
        depths = sort(collect(values(metadata(csd)[:depths])))
        csd = set(csd, Depth = depths)

        # * Get 0-100 ms window
        csd = csd[𝑡 = 0u"s" .. 0.1u"s"]
        # * Normalize CSD so that proms are meaningful
        zcsd = ZScore(csd)(csd)

        dt = samplingperiod(csd)
        w = 5u"ms"
        w = ceil(Int, uconvert(NoUnits, w / dt))
        minima = map(eachslice(.-zcsd, dims = Depth)) do x
            findpeaks(x, w, minprom = 0.5)
        end
        proms = getindex.(minima, 2)
        minima = first.(minima) # 3rd is widths

        _minima = map(minima, lookup(minima, Depth)) do m, d
            ts = times(m)
            idxs = map(ts) do t
                ustrip(csd[𝑡 = At(t), Depth = At(d)]) < 0 # Only select true 'sinks'
            end
            isempty(idxs) ? [] : m[idxs]
        end
        _minima = filter(!isempty, minima)

        allmins = map(_minima, lookup(_minima, Depth)) do m, d
            zip(lookup(m, 𝑡), d) |> collect
        end
        allmins = vcat(allmins...)

        # * Find the depth of the first sink
        t, idx = findmin(first.(allmins))
        if t > 0.05u"s"
            @warn "No L4 center found in the first 50 ms, using first minimum at $(t)"
        end
        L4 = last(allmins[idx])

        # * Cross check against layer names
        layernames = metadata(csd)[:layernames]
        if !contains(layernames[idx], "4")
            @warn "The first minimum at $(t) is not in L4, but in $(layernames[idx])"
        end

        if doplot
            fax = heatmap(csd |> ustripall, colormap = :turbo,
                          colorrange = symextrema(csd |> ustripall))
            map(ustripall(minima), lookup(minima)[1] .|> ustrip,
                ustripall(proms)) do m, depth, prom
                scatter!(lookup(m)[1], fill(depth, length(m)), color = juliapurple,
                         markersize = prom * 10)
            end
            Colorbar(fax.figure[1, 2], fax.plot)
            fax.axis.yreversed = true

            layernames = metadata(csd)[:layernames]
            if last(layernames) ∈ ["or", "scwm", "cing"]
                layernames[end] = layernames[end - 1]
            end
            plotlayerints!(fax.axis, layernames; width = 0.04)
            # fax.axis.yticks = (lookup(csd, Depth), layernames)
            # fax.axis.ylabel = "Depth (μm)"

            sessionid = metadata(csd)[:sessionid]
            structure = metadata(csd)[:structure]
            fax.axis.title = "$(structure) CSD (session $(sessionid))"
            mkpath(calcdir("plots", "csd", "$(sessionid)"))
            save(calcdir("plots", "csd", "$(sessionid)",
                         "$(structure).pdf"), fax)
        end
        return L4
    end
end

begin # * Layer identification for all sessions and structures
    l4s = progressmap(csd, structures) do c, structure
        @info "Finding L4 center for $(structure)"
        l4s = map(c) do csd_session
            find_L4_center(csd_session, true)
        end
        # l4s = ToolsArray(l4s, (SessionID(lookup(csd_structure, :sessionid)),))
        # D = @strdict l4s
        # tagsave(calcdir("csd", "l4_depth.jld2"), D; structure = structure)
    end
end

begin # * Compare inferred L4 location to anatomical L4
    structure = 1

    x = anatomical_l4s[structure] |> parent
    y = l4s[structure] |> parent

    cor(x, y)

    scatter(x, y) |> display
    f = Figure()
    ax = Axis(f[1, 1], title = "CSD vs Anatomical L4 depth")
    density!(ax, y,
             color = :red, label = "CSD L4")
    density!(ax, x,
             color = :blue, label = "Anatomical L4")
    axislegend(ax)
    display(f)
end

begin
    unidepths = 0:40:900
    Ls = map(layernames) do L
        map(L) do l
            if last(l) ∈ ["or", "scwm", "cing"]
                l[end] = l[end - 1]
            end
            l[Depth = Near(unidepths)]
        end
    end
    layernums = [[ToolsArray(SpatiotemporalMotifs.parselayernum.(parent(p)), dims(p))
                  for p in L] for L in Ls]
    layerints = map(layernums) do layernums
        map(unique(vcat(parent.(layernums)...))) do l
            depths = map(enumerate(layernums)) do (sid, ls)
                if sum(ls .== l) == 0
                    ma = sum(ls .== l - 1) == 0 ? 1.0 :
                         maximum(lookup(ls[ls .== l - 1], Depth))
                    mi = sum(ls .== l + 1) == 0 ? 0.0 :
                         minimum(lookup(ls[ls .== l + 1], Depth))
                    m = mean([ma, mi])
                    @info "Element $sid has no layer $(l)"
                    [m, m]
                else
                    collect(extrema(lookup(ls, Depth)[ls .== l]))
                end
            end
            Interval(extrema(vcat(depths...))...)
        end
    end
end

begin # * Paired change plots
    f = SixPanel()
    gs = subdivide(f, 3, 2)
    axs = map(enumerate(gs)) do (structure, g)
        y = anatomical_l4s[structure] |> parent
        x = l4s[structure] |> parent

        r = corspearman(x, y)
        ax = Axis(g, xticks = ([1, 2], ["CSD L4", "Anatomical L4"]),
                  ylabel = "Depth(μm)",
                  title = structures[structure], yreversed = true)

        _x = 1.0 .+ 0.15 .* abs.(randn(length(x)))
        _y = 2.0 .- 0.15 .* abs.(randn(length(y)))

        axx = Axis(g; limits = ((0.5, 2.6), (0, 800)), yaxisposition = :right,
                   yreversed = true)

        hidespines!(axx)
        hidedecorations!(axx; ticklabels = false)
        axx.xticklabelsvisible = false
        linkxaxes!(ax, axx)
        linkyaxes!(ax, axx)
        plotlayerints!(axx, layerints[structure], flipside = true)

        scatter!(ax, _y, y; color = structurecolors[structure], markersize = 10)
        scatter!(ax, _x, x; color = :grey, markersize = 10)

        segs = zip(zip(_x, _y) .|> Point2f, zip(x, y) .|> Point2f)
        map(segs) do seg
            lines!(ax, seg...; color = :black, linewidth = 0.25)
        end
        violin!(ax, fill(2 - 0.025, length(y)), y; side = :right, width = 1.0,
                show_median = true, color = (structurecolors[structure], 0.8),
                bandwidth = 15)
        violin!(ax, fill(1 + 0.025, length(x)), x; side = :left, width = 1.0,
                show_median = true, color = (:grey, 0.5))
        text!(ax, (0.6, 0.84), text = "ρ = $(round(r, digits = 2))", space = :relative,
              fontsize = 22)
        return ax
    end
    addlabels!(f, labelformat)
    display(f)
    wsave(plotdir("figS6", "figS6.pdf"), f)
end
