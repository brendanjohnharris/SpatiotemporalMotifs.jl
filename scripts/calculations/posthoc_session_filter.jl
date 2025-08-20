#! /bin/bash
# -*- mode: julia -*-
#=
exec julia +1.10.10 -t auto --color=yes "${BASH_SOURCE[0]}" "$@"
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

session_table = load(calcdir("plots", "session_table.jld2"), "session_table")
oursessions = deepcopy(session_table)

begin
    Q = calcquality()[Structure = At(structures)]
    quality = mean(Q[stimulus = At(stimulus)])
    out = load_calculations(Q; stimulus, vars = [:r])
end

begin # * Layer consistency
    sessions = [getindex.(o, :sessionid) for o in out]
    layers = [getindex.(o, :layernames) for o in out]
    depths = [getindex.(o, :streamlinedepths) for o in out]

    # ? Has most depths
    has_good_depths = reduce(.&,
                             [(maximum.(d) .> 0.90) .& (minimum.(d) .< 0.10)
                              for d in depths])

    # ? No missing layers
    nomissinglayers = [[unique(parselayernum.(_l)) for _l in unique.(l)] for l in layers]
    nomissinglayers = [[all(1:5 .∈ [_l]) for _l in _ls] for _ls in nomissinglayers]
    nomissinglayers = reduce(.&, nomissinglayers)

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
    goodsessions = sessions[1][decision .& nomissinglayers .& has_good_depths]
end

begin # * Task performance
    path = calcdir("calculations")
    performance = load_performance(; path)
    performantsessions = performance.sessionid[performance.mean_dprime .> 1]
    goodsessions = [g for g in goodsessions if g in performantsessions]
end

begin # * Manual blacklist
    blacklist = [1128520325] # * This session has no channels in layer 3 (after downsampling for unified depths)
    # push!(blacklist, 1093642839) # * This session has no theta
    # push!(blacklist, 1087992708) # * This session has no theta
    goodsessions = [g for g in goodsessions if g ∉ blacklist]
end

begin
    newsessions = subset(oursessions, :ecephys_session_id => ByRow(∈(goodsessions)))
    mkpath(calcdir("plots"))
    tagsave(calcdir("plots", "posthoc_session_table.jld2"),
            Dict("session_table" => newsessions))
    write(calcdir("plots", "posthoc_session_table.json"), JSON.json(newsessions))
end
# Read the dataframe as read("$(@__DIR__)/../plots/session_table.json", String) |> JSON.parse |> DataFrame

begin # * Formatted subject table
    using Dates
    file = calcdir("plots", "posthoc_session_table.jld2")
    newsessions = load(file, "session_table")
    experimental_model_table = newsessions[:,
                                           [:mouse_id,
                                               :ecephys_session_id,
                                               :date_of_acquisition,
                                               :equipment_name,
                                               :genotype,
                                               :sex,
                                               :age_in_days,
                                               :session_number
                                           ]]

    # * Add number of LFP channels
    D = map(out, structures) do O, structure
        d = map(O) do o
            r = o[:r]
            num_units = length(o[:spiketimes])
            m = DimensionalData.metadata(r)
            DataFrame(; structure, ecephys_session_id = m[:sessionid],
                      num_channels = size(r, Depth), num_units)
        end
        vcat(d...)
    end
    D = vcat(D...)

    D = leftjoin(experimental_model_table, D, on = :ecephys_session_id)
    D.date_of_acquisition = map(D.date_of_acquisition) do x
        x = split(x, " ") |> first
        x = DateTime(x, "yyyy-mm-dd")
    end
    D = DataFrames.groupby(D, :genotype)

    meanstd(x) = "$(round(mean(x), sigdigits=2)) ± $(round(std(x), sigdigits=2))"
    D = DataFrames.combine(D, :genotype => :genotype,
                           :num_channels => meanstd => "Num. channels",
                           :num_units => meanstd => "Num. units",
                           :date_of_acquisition => extrema => "Dates of acquisition",
                           :age_in_days => meanstd => "Age (days)",
                           [:mouse_id, :sex] => ((x, y) -> length(unique(x[y .== "M"]))) => "Num. male",
                           [:mouse_id, :sex] => ((x, y) -> length(unique(x[y .== "F"]))) => "Num. female") |>
        unique
    D.var"Dates of acquisition" .= [Dates.format.(d, ["yyyy-mm-dd"])
                                    for d in D.var"Dates of acquisition"]
    D.var"Dates of acquisition" .= join.(D.var"Dates of acquisition", [" to "])
    open(calcdir("experimental_model_table.csv"), "w") do io
        write(io, join(names(D), ","))
        write(io, "\n")
        for d in eachrow(D)
            write(io, join(string.(collect(d)), ","))
            write(io, "\n")
        end
    end
end
