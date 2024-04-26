module SpatiotemporalMotifs
import SpatiotemporalMotifs as SM
using DrWatson.Dates
using DimensionalData
using PythonCall
using Conda
using StatsBase
using JLD2
using CairoMakie
using Foresight
set_theme!(foresight(:physics))
using Statistics
using IntervalSets
using Normalization
using Distributed
using DSP
using Peaks
using TimeseriesTools
using TimeseriesTools.Unitful
using ModulationIndices
using DataFrames
using ComplexityMeasures
using JSON
using DrWatson
import AllenNeuropixels as AN
import AllenNeuropixelsBase as ANB
using Term.Progress

import DimensionalData: metadata

include("Utils.jl")
include("Calculations.jl")
include("Plots.jl")

export classifier,
       classify_kfold,
       commondepths,
       crossvalidate,
       load_calculations,
       load_performance,
       parselayernum,
       plotlayerints!,
       plotlayermap!,
       structurecolormap,
       structurecolors,
       structures,
       lfpcolormap,
       amplitudecolormap,
       phasecolormap,
       unify_calculations,
       calcquality,
       layernum2name,
       plotdir,
       savepath,
       symextrema,
       produce_out,
       produce_uni,
       plotstructurecenters!,
       load_unitdepths,
       regressor,
       regress_kfold,
       hierarchy_scores,
       visual_cortex_layout,
       ppc, spc!, sac, sac!, Bins,
       @preamble

end # module
