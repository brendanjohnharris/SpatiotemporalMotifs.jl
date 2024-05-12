#! /bin/bash
#=
exec julia -t auto "${BASH_SOURCE[0]}" "$@"
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
structurecolors, structures, commondepths, parselayernum,
load_calculations, unify_calculations,
plotlayerints!, plotlayermap!,
@preamble
set_theme!(foresight(:physics))
Random.seed!(32)

stimulus = r"Natural_Images"
datafile = datadir("theta_waves_task.jld2")

session_table = load(datadir("session_table.jld2"), "session_table")
oursessions = session_table.ecephys_session_id

path = datadir("calculations")
Q = calcquality(path)[structure = At(structures)]
quality = mean(Q[stimulus = At(stimulus)])

# * Theta wavenumber
vars = [:k]
config = @strdict stimulus vars
data, file = produce_or_load(produce_uni, config, datadir(); filename = savepath)
uni = data["uni"]
begin # * Grand  unified layers
    layerints = getindex.(uni, :layerints)
    layerints = map(zip(layerints...)) do ints
        mi = minimum.(ints) |> minimum
        ma = maximum.(ints) |> maximum
        return mi .. ma
    end
end

tagsave(datadir("grand_unified_layers.jld2"), (@strdict layerints))
