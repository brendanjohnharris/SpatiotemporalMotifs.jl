#! /bin/bash
#=
exec julia -t auto "${BASH_SOURCE[0]}" "$@"
=#
using DrWatson
@quickactivate "SpatiotemporalMotifs"
using SpatiotemporalMotifs
using IntervalSets
using Normalization
using Cairomakie # ! Needs to be installed in home project
import CairoMakie.record
using DimensionalData
using Foresight
using TimeseriesTools
using AllenNeuropixels
import AllenNeuropixels as AN
import AllenNeuropixelsBase as ANB

# * Get data
begin
    sessionid = SpatiotemporalMotifs.DEFAULT_SESSION_ID
    trial = SpatiotemporalMotifs.DEFAULT_TRIAL_NUM
    stimulus = r"NaturalImages"
end
begin
    files = map(structures) do structure
        file = savepath((; sessionid, stimulus, structure), "jld2")
        file = datadir("calculations", file)
    end
    θs = map(files) do file
        file = jldopen(file, "r")
        θ = file["θ"][:, :, trial]
    end
end
begin
    start = 100
    step = 4
    N = 4000
    LFP = ANB.rectifytime.(_LFP)
    LFP = lowpass.(LFP, pass)
    LFP = getindex.(LFP, [range(; start, step, length = N)], [:])
    channels = dims.(LFP, 2) .|> collect
    flashes = AN.getstimuli(session)[AN.getstimuli(session).stimulus_name .== "flash_250ms",
                                     :]
    flashes = [Interval(x...) for x in zip(flashes.start_time, flashes.end_time)]
end

begin # Set up plot
    f = Figure(; size = (1080, 900), light_direction = Vec3(100000.0, 100000.0, 0.0))
    ax = Axis3(f[1, 1]; aspect = :data, title = "0 s", titlegap = -100,
               xzpanelvisible = false, yzpanelvisible = false, xypanelvisible = false)
    hidedecorations!(ax)
    ax.xspinesvisible = ax.yspinesvisible = ax.zspinesvisible = false
    probeids, c, p = AN.Plots.plotbrain!(ax, session; dark = false,
                                         probestyle = :meshscatter,
                                         channels, markersize = 100.0, fontsize = 15.0)
    f
end

begin # ? Animate activity from an LFP trace. An example is saved in PlotBrain.jld2
    framerate = 60
    ts = size(LFP[1], 1)
    truets = times(LFP[1])
    idxs = indexin(probeids, getindex.(DimensionalData.metadata.(LFP), :probeid))
    X = collect.(LFP[idxs]) # Sort to probe order in plot

    Nc = length.([_c[] for _c in c])

    colormap = cgrad(:binary)
    X = [MinMax(x, dims = 1)(x) for x in X]
    Nt = unique(size.(X, 1))
    length(Nt) > 1 && error("All provided arrays must have the same number of rows")
    Nt = first(Nt)
    @withprogress record(f, "PlotBrain.mp4", 1:Nt; framerate) do t
        for i in eachindex(c)
            c[i][] = colormap[X[i][t, :]]
        end
        if t < Nt ÷ 2
            ax.azimuth[] += 0.0015
            ax.elevation[] += 0.0005
        else
            ax.azimuth[] -= 0.0015
            ax.elevation[] -= 0.0005
        end
        if any(truets[t] .∈ flashes)
            f.scene.backgroundcolor[] = darkbg * 2
        else
            f.scene.backgroundcolor[] = darkbg
        end
        ax.title[] = "$(round(truets[t], digits=1)) s"
        @logprogress t / Nt
    end
end

# ? Do it in a loop for GLMakie interactivity
if false
    map(ts) do t
        for i in eachindex(c)
            c[i][] = colormap[X[i][t, :]]
        end
    end
    sleep(1 / fs)
end
