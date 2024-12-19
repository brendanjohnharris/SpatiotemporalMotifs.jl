#! /bin/bash
#=
exec julia -t auto "${BASH_SOURCE[0]}" "$@"
=#
using DrWatson
@quickactivate "SpatiotemporalMotifs"

import DimensionalData: metadata
using SpatiotemporalMotifs
@preamble
set_theme!(foresight(:physics))
Random.seed!(42)
using CairoMakie
using GraphMakie
using GraphMakie.Graphs
using SimpleWeightedGraphs
using Downloads
using LinearAlgebra
using FileIO

## This file follows the code from Siegle et al. 2021: 'Survey of spiking in the mouse visual cortex reveals functional hierarchy'
# ? See https://github.com/AllenInstitute/neuropixels_platform_paper/blob/master/Figure2/comparison_anatomical_functional_connectivity_final.ipynb

## * The anatomical hierarchy
thr = 8
areas = ["VISp", "VISl", "VISrl", "VISal", "VISpm", "VISam"]
ascore = [-0.357, -0.093, -0.059, 0.152, 0.327, 0.441]
amatrix = [abs(a - b) for (a, b) in Iterators.product(values(ascore), values(ascore))] |>
          Symmetric
amatrix = amatrix ./ maximum(amatrix)
amatrix = 1.0 .- amatrix
amatrix = amatrix .- Diagonal(amatrix)
while sum(amatrix .> 0) > thr * 2
    amatrix[findmin(x -> x == 0 ? Inf : x, amatrix)[2]] = 0
end
amatrix = amatrix .- minimum(filter(>(0), amatrix)) / 1.25
amatrix = amatrix ./ maximum(amatrix)
ag = SimpleWeightedGraph(amatrix)

## * Load functional data
dir = tempdir()
baseurl = "https://github.com/AllenInstitute/neuropixels_platform_paper/raw/master/data/processed_data"
f = "FFscore_grating_10.npy"
Downloads.download(joinpath(baseurl, f), joinpath(dir, f))

## * Load functional data
FF_score, FF_score_b, FF_score_std = eachslice(load(joinpath(dir, f)), dims = 1)
fmatrix = Symmetric((FF_score .+ FF_score') ./ 2)
fmatrix = abs.(fmatrix)
fmatrix = fmatrix ./ maximum(fmatrix)
fmatrix = 1.0 .- fmatrix
fmatrix = fmatrix .- Diagonal(fmatrix)
while sum(fmatrix .> 0) > thr * 2
    fmatrix[findmin(x -> x == 0 ? Inf : x, fmatrix)[2]] = 0
end
fmatrix = fmatrix .- minimum(filter(>(0), fmatrix)) / 1.25
fmatrix = fmatrix ./ maximum(fmatrix)
fg = SimpleWeightedGraph(fmatrix)

## * Plot the hierarchy
begin
    colormap = getindex.([cgrad(:inferno)], (1:6) ./ (7))
    f = Figure()
    ax = Axis(f[1, 1])
    plotstructurecenters!(ax, ag; colormap, structures = areas, arrow_size = 50)
    f
end

graphplot(SimpleWeightedDiGraph([0 -1; 1 0]), arrow_size = 50, arrow_shift = :end)
