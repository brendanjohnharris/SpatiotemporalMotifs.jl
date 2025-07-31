
#! /bin/bash
#=
exec julia +1.10.10 -t auto "${BASH_SOURCE[0]}" "$@"
=#
using DrWatson
@quickactivate "SpatiotemporalMotifs"
import AllenNeuropixels as AN
import AllenNeuropixelsBase as ANB
using Statistics
using JLD2
using TimeseriesTools
using SpatiotemporalMotifs
using GeometryBasics, RPRMakie, RadeonProRender# Must be installed in home project
RPRMakie.activate!(; plugin = RPRMakie.Northstar, iterations = 500)
using Colors, FileIO
using Colors: N0f8

begin
    using Distributed
    using USydClusters
    project = Base.active_project()
    procs = USydClusters.Physics.addprocs(6; ncpus = 12, mem = 50,
                                          walltime = 4, project, qsub_flags = `-q l40s`)
end

@everywhere begin
    using DrWatson
    import AllenNeuropixels as AN
    import AllenNeuropixelsBase as ANB
    using Statistics
    using JLD2
    using TimeseriesTools
    using SpatiotemporalMotifs
    using GeometryBasics, RPRMakie, RadeonProRender# Must be installed in home project
    RPRMakie.activate!(; plugin = RPRMakie.Northstar, iterations = 500)
    using Colors, FileIO
    using Colors: N0f8

    begin
        sessionid = SpatiotemporalMotifs.DEFAULT_SESSION_ID
        trial = SpatiotemporalMotifs.DEFAULT_TRIAL_NUM
        stimulus = r"Natural_Images"
        if !isdefined(Main, :session) || !(AN.getid(session) == sessionid)
            session = AN.Session(sessionid)
        end
        # * Find the trial time by loading VISp calculation
        structure = "VISp"
        file = Dict("sessionid" => sessionid, "stimulus" => stimulus,
                    "structure" => structure)
        file = savepath(file, "jld2")
        file = jldopen(calcdir("calculations", file), "r")
        trials = file["trials"]
        times = trials[trial, :change_time_with_display_delay]
        times = times - 0.25 .. times + 0.75
        probes = AN.getprobes(session)
        LFP = map(probes.id) do probeid
            AN.getlfp(session, probeid; times)
        end
        LFP = matchdim(LFP; dims = Ti)
        LFP = cat(LFP..., dims = Dim{:channel})
        LFP = bandpass(LFP, 4 .. 10)
        LFP = normalize(LFP, MinMax; dims = 1) .* 1.5
        LFP = LFP[1:5:end, :]
    end

    function rot(angle, azimuth)
        return [cos(angle)*cos(azimuth) -sin(angle) cos(angle)*sin(azimuth);
                sin(angle)*cos(azimuth) cos(angle) sin(angle)*sin(azimuth);
                -sin(azimuth) 0 cos(azimuth)]
    end
    function render_brainframe(frame;
                               format = plotdir("videoS1", "glass_brain.jpg"))
        begin # * Parameters
            radiance = 2
            scale = 1e-5
            translation = [0.02, 0.0, 0.00]
            centre = Vec3f(0, 0, 0.05)
            eyeposition = centre .+ rot(-π / 3, π / 6) * [-0.2, 0, 0] |> Vec3f
            ambient = 0.5
            background = load(RPR.assetpath("studio032.exr"))
            lights = [EnvironmentLight(ambient, background),
                PointLight(3.0 .* rot(-0.9π, π / 6) * [-0.3, 0, 0] |> Vec3f,
                           RGBf(radiance, radiance, radiance * 1.1))]
            # markersize = 1000
        end

        fig = Figure(; size = (1920, 1080)) #(3840, 2160));
        ax = LScene(fig[1, 1]; show_axis = false, scenekw = (; lights = lights))
        ax = Scene(; size = (1920, 1080), lights = lights)

        screen = RPRMakie.Screen(ax)
        matsys = screen.matsys

        function stainedglass(color = :white, transparency = 0.85; kwargs...)
            RPR.Glass(matsys;
                      reflection_mode = 1,
                      refraction_thin_surface = false,
                      #   refraction_caustics = true,
                      transparency = Vec4f(transparency),
                      refraction_absorption_distance = 0.0,
                      refraction_ior = 1.6,
                      reflection_ior = 1.6,
                      reflection_color = Makie.to_color(color),
                      refraction_color = Makie.to_color(color),
                      refraction_absorption_color = Makie.to_color(color),
                      refraction_weight = Vec4f(1.0),
                      reflection_weight = Vec4f(1.0),
                      kwargs...)
        end
        function whitehot(T, color, emission_color = color; kwargs...)
            stainedglass(color, 0.0; emission_weight = Vec4f(T),
                         emission_color = Makie.to_color(emission_color), kwargs...)
        end

        # chrome = RPR.Chrome(matsys)
        # diffuse = RPR.DiffuseMaterial(matsys)
        # plastic = RPR.Plastic(matsys)

        backdrop = load(Makie.assetpath("matball_floor.obj"))
        backdrop.position .= map(backdrop.position) do x
            x = x .* [1, -1, 1] .+ [2.7, -3, -0.05]
        end
        mesh!(ax, backdrop; color = RGBf(0.8))

        begin # Plot structures
            D = ANB.getstructureidmap()
            # ? Plot the whole brain
            material = stainedglass(:white, 0.95; refraction_thin_surface = false)
            pw = AN.Plots.plotbrainstructure!(ax, 997; hemisphere = :both, scale, material)
            offset = .-mean(mean.(collect.(pw[1][]))) .+ translation .+ centre # Place the brain front and center
            translate!(pw, offset)

            targets = unique(["VISp",
                                 "VISl",
                                 "VISrl",
                                 "VISal",
                                 "VISpm",
                                 "VISam",
                                 "CA",
                                 "TH"])
            D = ANB.getstructureidmap()
            ids = getindex.((D,), targets)
            for id in ids # Plot structures
                color = AN.getstructurecolor(id)
                p = AN.Plots.plotbrainstructure!(ax, id; hemisphere = :right, scale,
                                                 material = stainedglass(color))
                translate!(p, offset)
            end
        end

        begin # * Plot probes
            channels = AN.getchannels(session)
            findchannels = unique(channels.probe_id)
            findchannels = Dict(findchannels .=> [channels[channels.probe_id .== p, :].id
                                 for p in findchannels])

            markercolor = :black
            glowcolor = colorant"#f9f37c" #:white
            # frame = Observable(0)
            # t = lift(i) do i
            #     lookup(LFP, Ti)[frame + 1]
            # end
            t = Observable(lookup(LFP, Ti)[frame + 1])
            T = lift(t) do t
                LFP[Ti = At(t)]
            end
            materials = lift(T) do T
                whitehot.(T, markercolor, glowcolor)
            end
            cp = map(enumerate(collect(findchannels))) do (i, fc)
                probeid = first(fc)
                xyz = AN.getchannelcoordinates(session, probeid)
                xyz = filter(xyz) do cx
                    first(cx) in lookup(LFP, :channel)
                end
                channels = keys(xyz) |> collect
                points = Point3f.(values(xyz))
                points = map(points) do x
                    AN.Plots.ccftransform(x) .* scale .+ offset
                end

                idxs = indexin(channels, lookup(LFP, :channel))
                channelmaterials = map(idxs) do i
                    lift(materials -> getindex(materials, i), materials)
                end
                function plotfunc()
                    map(points, channelmaterials) do p, m
                        meshscatter!(ax, p; markersize = 1.5e-3,
                                     material = m)
                    end
                end
                return channelmaterials, plotfunc
            end
            plotfuncs = last.(cp)
            ps = [pf() for pf in plotfuncs]
        end
        # t[] = lookup(LFP, Ti)[300]

        cam3d!(ax)
        cam = cameracontrols(ax) # ax.scene
        cam.eyeposition[] = eyeposition
        cam.lookat[] = centre
        cam.upvector[] = Vec3f(0, 0, 1)
        cam.fov[] = 32

        image = colorbuffer(screen)
        file, ext = splitext(format)
        file = file * "$(lpad(frame, 4, '0'))" * ext
        @info "Saving to $file"
        wsave(file, image)
    end
end

if false
    @spawnat 2 render_brainframe(200)
end
begin
    pmap(render_brainframe, 1:250; distributed = true)
end
