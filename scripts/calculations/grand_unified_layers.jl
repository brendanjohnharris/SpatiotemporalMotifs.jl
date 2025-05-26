#! /bin/bash
#=
exec julia +1.10.9 -t auto "${BASH_SOURCE[0]}" "$@"
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

stimulus = r"Natural_Images"
datafile = datadir("theta_waves_task.jld2")

session_table = load(datadir("posthoc_session_table.jld2"), "session_table")
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

tagsave(datadir("grand_unified_layers.jld2"), (@strdict layerints))
tagsave(datadir("plots", "grand_unified_layers.jld2"), (@strdict layerints))
