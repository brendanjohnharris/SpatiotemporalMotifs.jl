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
using Term
using TimeseriesSurrogates
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
       send_calculations,
       send_thalamus_calculations,
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
       produce_out,
       produce_uni,
       plotstructurecenters!,
       load_unitdepths,
       regressor,
       regress_kfold,
       hierarchy_scores,
       visual_cortex_layout,
       plot_visual_cortex!,
       ppc, spc!, sac, sac!, Bins, pac,
       confidence, quartiles,
       bootstrapmedian, bootstrapaverage,
       progressmap,
       plotspectrum!, fooof, tortinset!,
       @preamble

end # module
