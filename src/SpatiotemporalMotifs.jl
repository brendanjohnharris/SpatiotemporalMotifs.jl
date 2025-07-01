module SpatiotemporalMotifs
import SpatiotemporalMotifs as SM
using DrWatson.Dates
using DimensionalData
using PythonCall
using Conda
# using StatsBase
using JLD2
using CairoMakie
using Foresight
set_theme!(foresight(:physics))
using Statistics
# using IntervalSets
using Normalization
using Distributed
using DSP
using Peaks
using Term
using TimeseriesSurrogates
# using TimeseriesTools
# using TimeseriesTools.Unitful
using ModulationIndices
using DataFrames
using ComplexityMeasures
using JSON
# using DrWatson
import AllenNeuropixels as AN
import AllenNeuropixelsBase as ANB
using Term.Progress

# * Easier to reexport these for workers
using Reexport
@reexport using DimensionalData
@reexport using StatsBase
@reexport using DrWatson
@reexport using Unitful
@reexport using IntervalSets

import DimensionalData: metadata

include("Utils.jl")
include("Calculations.jl")
include("Plots.jl")

export classifier,
       tuneclassifier,
       classify_kfold,
       commondepths,
       crossvalidate,
       load_calculations,
       load_performance,
       send_calculations,
       send_powerspectra,
       parselayernum,
       plotlayerints!,
       plotlayermap!,
       depthticks,
       structurecolormap,
       structurecgrad,
       structurecolors,
       layers,
       layercolors,
       structures,
       lfpcolormap,
       defaultcolormap,
       amplitudecolormap,
       phasecolormap,
       unify_calculations,
       calcquality,
       layernum2name,
       plotdir,
       savepath,
       symextrema,
       load_calculations,
       load_uni,
       plotstructurecenters!,
       load_unitdepths,
       regressor,
       regress_kfold,
       hierarchy_scores,
       visual_cortex_layout,
       plot_visual_cortex!,
       ppc, spc!, sac, sac!, HistBins, pac,
       confidence, quartiles, bootstrapmean,
       bootstrapmedian, bootstrapaverage, hierarchicalkendall,
       progressmap,
       plotspectrum!, fooof, tortinset!,
       @preamble

end # module
