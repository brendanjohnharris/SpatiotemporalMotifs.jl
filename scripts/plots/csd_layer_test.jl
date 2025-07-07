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

begin # * Load the CSD flashes
end
# begin # * Load the CSD passive natural images
#     session_table = load(datadir("posthoc_session_table.jld2"), "session_table")
#     oursessions = session_table.ecephys_session_id
#     path = datadir("calculations")
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
#     tagsave(datadir("average_csd.jld2"), D)
# end

begin # * Save average CSD for each session
    conf = Dict("stimulus" => "Natural_Images_passive_nochange")

    csd, csd_file = produce_or_load(copy(conf), datadir("plots");
                                    filename = savepath,
                                    prefix = "csd") do conf
        stimulus = conf["stimulus"]

        session_table = load(datadir("posthoc_session_table.jld2"), "session_table")
        oursessions = session_table.ecephys_session_id
        path = datadir("calculations")
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
        return Dict(structures .=> csd)
    end
end

begin # * Try just averaging nochange results
    heatmap(csd["VISp"][2] |> ustripall,
            axis = (; yreversed = true, limits = ((0, 0.1), nothing)),
            colorrange = (-1e-8, 1e-8), colormap = :turbo)
end

function find_L4_center(csd::MultivariateTimeSeries, doplot = false)
    # * Get 0-100 ms window
    csd = csd[𝑡 = 0u"s" .. 0.1u"s"]
    # * Normalize CSD so that proms are meaningful
    zcsd = ZScore(csd)(csd)

    dt = samplingperiod(csd)
    w = 15u"ms"
    w = ceil(Int, uconvert(NoUnits, w / dt))
    minima = map(eachslice(.-zcsd, dims = Depth)) do x
        findpeaks(x, w, minprom = 0.3)
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
        fax.axis.yticks = (lookup(csd, Depth), layernames)

        sessionid = metadata(csd)[:sessionid]
        structure = metadata(csd)[:structure]
        fax.axis.title = "$(structure) CSD (session $(sessionid))"
        mkpath(datadir("plots", "csd", "$(sessionid)"))
        save(datadir("plots", "csd", "$(sessionid)",
                     "$(structure).pdf"), fax)
    end
    return L4
end

begin # * Identify Layer 4
    x = csd["VISp"][2]
    l4 = find_L4_center(x, true)
end

begin # * Layer identification for all sessions and structures
    l4s = map(structures) do structure
        @info "Finding L4 center for $(structure)"
        csd_structure = csd[structure]
        l4s = map(csd_structure) do csd_session
            find_L4_center(csd_session, true)
        end
        # l4s = ToolsArray(l4s, (SessionID(lookup(csd_structure, :sessionid)),))
        # D = @strdict l4s
        # tagsave(datadir("csd", "l4_depth.jld2"), D; structure = structure)
    end
end
begin # * Pull out anatomical l4s
    session_table = load(datadir("posthoc_session_table.jld2"), "session_table")
    oursessions = session_table.ecephys_session_id
    path = datadir("calculations")
    Q = calcquality(path)[Structure = At(structures)]
    Q = Q[SessionID(At(oursessions))]
    @assert mean(Q[stimulus = At(r"Natural_Images")]) == 1
    out = load_calculations(Q; stimulus = r"Natural_Images", vars = [])
    true_l4s = map(out) do O
        l4s = map(O) do o
            l4 = contains.(o[:layernames], ["4"])
            l4 = o[:streamlinedepths][l4] |> mean # 'center' of l4
        end
        l4s = ToolsArray(l4s, (SessionID(oursessions),))
    end
end
begin # * Compare inferred L4 location to anatomical L4
    structure = 1

    x = l4s[structure][SessionID = At(oursessions)] |> parent
    y = true_l4s[structure][SessionID = At(oursessions)] |> parent

    cor(x, y)

    scatter(x, y) |> display
    f = Figure()
    ax = Axis(f[1, 1], title = "Inferred vs Anatomical L4 depth")
    hist!(ax, x,
          color = :blue, label = "Inferred L4")
    hist!(ax, y,
          color = :red, label = "Anatomical L4")
    axislegend(ax)
    display(f)
end

begin
    s = 2
    i = 19
    begin # * Our csd
        csdf = outf[s][i][:csd][𝑡 = 1:1562]
        csdf = mapslices(median, csdf, dims = :changetime)
        csdf = dropdims(csdf, dims = :changetime)
        heatmap(csdf |> ustripall,
                axis = (; yreversed = true, limits = ((0, 0.1), nothing)))
    end
    begin # * Our csd
        csd = out[s][i][:csd][𝑡 = 1:1562]
        csd = mapslices(median, csd, dims = :changetime)
        csd = dropdims(csd, dims = :changetime)
        heatmap(csd |> ustripall,
                axis = (; yreversed = true, limits = ((0, 0.1), nothing)))
    end

    begin
        f = Figure()
        ax = Axis(f[1, 1], yreversed = true)
        x = mean([csd, csdf])
        heatmap!(ax, ustripall(x))
        ax.limits = ((0, 0.1), nothing)
        display(f)
    end
end

begin
    l4 = find_L4_center(parent(csd), lookup(csd, :Depth); window = (0, 50), fs = 1250)
    f = Figure()
    ax = Axis(f[1, 1], title = "L4 center depth", yreversed = true)
    heatmap!(ax, ustripall(csd))
    hlines!(ax, [l4], color = :red, label = "L4 center")
    ax.limits = ((0, 0.25), nothing)
    display(f)
end
