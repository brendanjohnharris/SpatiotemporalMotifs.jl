#! /bin/bash
#=
exec julia -t auto "${BASH_SOURCE[0]}" "$@"
=#
using DrWatson
@quickactivate "SpatiotemporalMotifs"

import TimeseriesTools: freqs
using FileIO
using Unitful
using JSON
import CairoMakie.Axis
using SpatiotemporalMotifs
@preamble
set_theme!(foresight(:physics))

stimulus = r"Natural_Images"
vars = [:r]

session_table = load(datadir("session_table.jld2"), "session_table")
oursessions = deepcopy(session_table)

begin
    path = datadir("calculations")
    Q = calcquality(path)[structure = At(structures)]
    quality = mean(Q[stimulus = At(stimulus)])
    config = @strdict stimulus vars
    data, file = produce_or_load(produce_out(Q), config, datadir(); filename = savepath,
                                 prefix = "out")
    out = data["out"]
end

begin # * Layer consistency
    sessions = [getindex.(o, :sessionid) for o in out]
    layers = [getindex.(o, :layernames) for o in out]
    depths = [getindex.(o, :streamlinedepths) for o in out]

    layerints = map(layers, depths) do ls, ds
        map(ls, ds) do l, d
            map(SpatiotemporalMotifs.layers) do _l
                finds = d[occursin.([_l], l)]
                isempty(finds) ? (NaN, NaN) : extrema(finds)
            end |> Iterators.flatten |> collect
        end
    end

    function deviation(x) # * Distance in IQRs from the median
        μ = median(x)
        σ = nansafe(iqr)(x)
        return abs.(x .- μ) ./ σ
    end
    layerdeviations = map(layerints) do ints
        devs = map(eachindex(ints[1])) do i
            deviation(getindex.(ints, i))
        end
        devs = cat(devs..., dims = 2)
        devs = nansafe(maximum; dims = 2)(devs)[:]
    end
    decision = maximum(hcat(layerdeviations...), dims = 2)[:] # Maximum value of IQR distance from median for any layer boundary in any region
    @assert all([sessions[1]] .== sessions)
    decision = decision .≤ 3.0
    goodsessions = sessions[1][decision .& nomissinglayers]
end

begin # * Task performance
    path = datadir("calculations")
    performance = load_performance(; path)
    performantsessions = performance.sessionid[performance.mean_dprime .> 1]
    goodsessions = [g for g in goodsessions if g in performantsessions]
end

begin # * Manual blacklist
    blacklist = [1128520325] # * This session has no channels in layer 3
    goodsessions = [g for g in goodsessions if g ∉ blacklist]
end

begin
    newsessions = subset(oursessions, :ecephys_session_id => ByRow(∈(goodsessions)))
    tagsave(datadir("posthoc_session_table.jld2"), Dict("session_table" => newsessions))
    write(datadir("posthoc_session_table.json"), JSON.json(newsessions))
end
# Read the dataframe as read("$(@__DIR__)/../session_table.json", String) |> JSON.parse |>
# DataFrame