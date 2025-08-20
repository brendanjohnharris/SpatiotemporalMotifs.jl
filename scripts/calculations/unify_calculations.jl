#! /bin/bash
#=
exec julia +1.10.10 -t auto "${BASH_SOURCE[0]}" "$@"
=#
# ? Expected execution time: 15 mins
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

begin # * Load quality
    path = calcdir("calculations")
    Q = calcquality(path)[Structure = At(structures)]
end

begin
    stimulus = r"Natural_Images"

    session_table = load(calcdir("plots", "posthoc_session_table.jld2"), "session_table")
    oursessions = session_table.ecephys_session_id

    # * Theta wavenumber
    vars = [:k]
    uni = load_uni(; stimulus, vars)
    begin # * Grand  unified layers
        layerints = getindex.(uni, :layerints)
        layerints = map(zip(layerints...)) do ints
            mi = minimum.(ints) |> minimum
            ma = maximum.(ints) |> maximum
            return mi .. ma
        end
    end

    tagsave(calcdir("plots", "grand_unified_layers.jld2"), (@strdict layerints))
end

# * Now do unified data
begin
    map(lookup(Q, :stimulus), eachslice(Q, dims = :stimulus)) do stimulus, q
        if !contains(stimulus |> string, "omission") &&
           !contains(stimulus |> string, "nochange") &&
           !contains(stimulus |> string, "passive") &&
           stimulus !== r"Natural_Images"
            @info "Calculating unified data for stimulus: $stimulus"
            uni = load_uni(; stimulus, vars = [:csd])
        end
        uni = []
        GC.gc()
    end
end
