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
@preamble
set_theme!(foresight(:physics))
Random.seed!(32)

sessionid = 1092466205
session = AN.Session(sessionid)
stimulus = "spontaneous"

Rs = map(SpatiotemporalMotifs.structures) do structure
    LFP = AN.formatlfp(session; structure, stimulus)
    channels = lookup(LFP, Chan)
    depths = AN.getchanneldepths(session, LFP; method = :probe)
    LFP = set(LFP, Chan => Depth(depths))
    LFP = rectify(LFP; dims = Depth)

    θ = bandpass(LFP, SpatiotemporalMotifs.THETA.*u"Hz")
    ϕ = analyticphase(θ)
    k = -centralderiv(ϕ, dims = Depth, grad = phasegrad)
    R = dropdims(mean(sign.(k), dims = Depth); dims = Depth)
end
begin
    f = Figure()
    ax = Axis(f[1, 1])
    boxplot!.([ax], fill.(1:length(Rs), length.(Rs)), collect.(Rs))
    # plot!((mean.(Rs)); color=SpatiotemporalMotifs.structurecolors)
    f
end
