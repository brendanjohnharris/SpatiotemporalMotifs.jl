#! /bin/bash
#=
exec julia -t auto "${BASH_SOURCE[0]}" "$@"
=#
using DrWatson
@quickactivate "SpatiotemporalMotifs"
using SpatiotemporalMotifs
@preamble
set_theme!(foresight(:physics))

using CairoMakie
using CairoMakie.FileIO
using Downloads
using Foresight
using AllenNeuropixels
import AllenNeuropixels as AN
import AllenNeuropixelsBase as ANB
set_theme!(foresight())

begin
    sessionid = 1128520325
    session = AN.Session(sessionid)
    probeids = AN.getprobes(session).id
    chrome = FileIO.load(Downloads.download("https://raw.githubusercontent.com/nidorx/matcaps/master/1024/C8D1DC_575B62_818892_6E747B.png"))
end

begin # Set up plot
    f = Figure(; size = (2560, 1440))
    ax = Axis3(f[1, 1];
               aspect = :data,
               xzpanelvisible = false,
               yzpanelvisible = false,
               xypanelvisible = false)
    ax.scene.lights = [DirectionalLight(RGBf(6, 6, 6), Vec3f(2, -1, 1))]
    # ax.scene.lights[2] = PointLight(Makie.RGB(1.0, 1.0, 1),
    #                                       Vec3(14530.263, 11445.357, 9516.364))
    hidedecorations!(ax)
    ax.xspinesvisible = ax.yspinesvisible = ax.zspinesvisible = false
    # push!(structurecolors, Foresight.keppel)
    probeids, c, p = AN.Plots.plotbrain!(ax, session;
                                         dark = false,
                                         probestyle = :meshscatter,
                                         markersize = 100.0,
                                         fontsize = 0.0,
                                         matcap = false,
                                         shading = Makie.MultiLightShading)
    ax.azimuth = 2.6
    ax.elevation = 0.24
    ss = AN.getprobestructures(session, structures)
    ss = getindex.([ss], probeids)
    is = indexin(ss, s)

    [c[i][] = fill(structurecolors[is[i]], length(c[i][])) for i in eachindex(is)]
    f
    wsave(plotdir("schematic", "plotbrain_static.png"), f)
end
