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

begin # * Load the CSD flashes
    session_table = load(datadir("posthoc_session_table.jld2"), "session_table")
    oursessions = session_table.ecephys_session_id
    path = datadir("calculations")
    Q = calcquality(path)[Structure = At(structures)]
    Q = Q[SessionID(At(oursessions))]
    @assert mean(Q[stimulus = At("flash_250ms")]) == 1
    out = load_calculations(Q; stimulus = "flash_250ms", vars = [:csd])
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
begin # * Save average CSD for each session for flashes
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
    csd = ToolsArray(csd, (Structure(structures),))
    D = @strdict csd
    tagsave(datadir("average_csd.jld2"), D)
end

begin # * Try just averaging nochange results
    csd = load(datadir("average_csd.jld2"), "csd")
    heatmap(csd[2][15] |> ustripall,
            axis = (; yreversed = true, limits = ((0, 0.1), nothing)),
            colorrange = (-2e-8, 2e-8))
end

function find_L4_center(csd::MultivariateTimeSeries;
                        window::Tuple{Real, Real} = (0, 50),
                        fs::Real = 1_000)
    # # * Sigmoid normalize csd in each depth
    # N = Sigmoid(csd; dims = 2)
    # csd = N(csd)
    # * Get 0-100 ms window
    csd = csd[𝑡 = 0u"s" .. 0.1u"s"]
    dt = samplingperiod(csd)
    w = 10u"ms"
    w = ceil(Int, uconvert(NoUnits, w / dt))
    minima = map(eachslice(.-csd, dims = Depth)) do x
        findpeaks(x, w)
    end
    proms = getindex.(minima, 2)
    minima = first.(minima) # 3rd is widths
    fax = heatmap(csd |> ustripall, colormap = :turbo,
                  colorrange = symextrema(csd |> ustripall))
    map(ustripall(minima), lookup(minima)[1] .|> ustrip, ustripall(proms)) do m, depth, prom
        scatter!(lookup(m)[1], fill(depth, length(m)), color = cornflowerblue,
                 markersize = prom * 2e9)
    end
    Colorbar(fax.figure[1, 2], fax.plot)
    fax.axis.yreversed = true
    fax.axis.yticks = (lookup(csd, Depth), metadata(csd)[:layernames])
    display(fax)
end

begin # * Identify Layer 4
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
