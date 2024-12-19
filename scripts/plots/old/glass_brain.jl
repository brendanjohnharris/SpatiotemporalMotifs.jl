
#! /bin/bash
#=
exec julia -t auto "${BASH_SOURCE[0]}" "$@"
=#
using DrWatson
@quickactivate "SpatiotemporalMotifs"
import AllenNeuropixels as AN
import AllenNeuropixelsBase as ANB
using Statistics
using JLD2
using LinearAlgebra
using Term
using TimeseriesTools
using SpatiotemporalMotifs
using GeometryBasics, RPRMakie, RadeonProRender# Must be installed in home project
RPRMakie.activate!(; plugin = RPRMakie.Northstar, iterations = 500)
using Colors, FileIO
using Colors: N0f8

if false
    using Distributed
    using USydClusters
    project = Base.active_project()
    procs = USydClusters.Physics.addprocs(6; ncpus = 4, mem = 8,
                                          walltime = 4, project, qsub_flags = `-q l40s`)
end

begin # @everywhere
    using DrWatson
    import AllenNeuropixels as AN
    import AllenNeuropixelsBase as ANB
    using Statistics
    using JLD2
    using LinearAlgebra
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
        file = jldopen(datadir("calculations", file), "r")
        trials = file["trials"]
        _times = trials[trial, :change_time_with_display_delay]
        times = _times - 1 .. _times + 1.5
        probes = AN.getprobes(session)
        LFP = map(probes.id) do probeid
            AN.getlfp(session, probeid; times)
        end
        LFP = matchdim(LFP; dims = 𝑡)
        LFP = cat(LFP..., dims = Chan)
        LFP = .-bandpass(LFP, 4 .. 10)
        LFP = normalize(LFP, MinMax; dims = 1) #.* 1.5
        times = _times - 0.25 .. _times + 0.75
        LFP = LFP[𝑡(times)]
        LFP = LFP[1:5:end, :]
    end

    function rot(angle, azimuth)
        return [cos(angle)*cos(azimuth) -sin(angle) cos(angle)*sin(azimuth);
                sin(angle)*cos(azimuth) cos(angle) sin(angle)*sin(azimuth);
                -sin(azimuth) 0 cos(azimuth)]
    end

    function rot_eyepos(eyeposition, lookat, azimuth_rad)
        direction = lookat - eyeposition
        direction_norm = normalize(direction)
        up = [0.9, 0, -0.4359]
        right = cross(up, direction_norm)
        rotation_matrix = [cos(azimuth_rad)+right[1]^2 * (1 - cos(azimuth_rad)) right[1] * right[2] * (1 - cos(azimuth_rad))-right[3] * sin(azimuth_rad) right[1] * right[3] * (1 - cos(azimuth_rad))+right[2] * sin(azimuth_rad);
                           right[2] * right[1] * (1 - cos(azimuth_rad))+right[3] * sin(azimuth_rad) cos(azimuth_rad)+right[2]^2 * (1 - cos(azimuth_rad)) right[2] * right[3] * (1 - cos(azimuth_rad))-right[1] * sin(azimuth_rad);
                           right[3] * right[1] * (1 - cos(azimuth_rad))-right[2] * sin(azimuth_rad) right[3] * right[2] * (1 - cos(azimuth_rad))+right[1] * sin(azimuth_rad) cos(azimuth_rad)+right[3]^2 * (1 - cos(azimuth_rad))]
        rotated_direction = rotation_matrix * direction
        new_eyeposition = lookat - rotated_direction
        return new_eyeposition
    end

    function render_brainframe(frame;
                               format = plotdir("glass_brain", "glass_brain.jpg"))
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

        # fig = Figure(; size = (1920, 1080)) #(3840, 2160));
        # ax = LScene(fig[1, 1]; show_axis = false, scenekw = (; lights = lights))
        ax = Scene(; size = (2160, 1440), lights = lights)

        screen = RPRMakie.Screen(ax)
        matsys = screen.matsys

        function stainedglass(color = :white, transparency = 0.8; kwargs...)
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
        function whitehot(T, color, emission_color = RGBf(T * 4.0, T * 0.8, T * 0.6);
                          kwargs...)
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
            material = stainedglass(:white, 0.9; refraction_thin_surface = false)
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
            # glowcolor = colorant"#f9f37c" #:white
            # frame = Observable(0)
            # t = lift(i) do i
            #     lookup(LFP, 𝑡)[frame + 1]
            # end
            t = Observable(lookup(LFP, 𝑡)[frame + 1])
            T = lift(t) do t
                LFP[Ti = At(t)]
            end
            materials = lift(T) do T
                whitehot.(T, markercolor)# , glowcolor)
            end
            cp = map(enumerate(collect(findchannels))) do (i, fc)
                probeid = first(fc)
                xyz = AN.getchannelcoordinates(session, probeid)
                xyz = filter(xyz) do cx
                    first(cx) in lookup(LFP, Chan)
                end
                channels = keys(xyz) |> collect
                points = Point3f.(values(xyz))
                points = map(points) do x
                    AN.Plots.ccftransform(x) .* scale .+ offset
                end

                idxs = indexin(channels, lookup(LFP, Chan))
                channelmaterials = map(idxs) do i
                    lift(materials -> getindex(materials, i), materials)
                end
                function plotfunc()
                    map(points, channelmaterials) do p, m
                        meshscatter!(ax, p; markersize = 1.5e-3,
                                     material = m, colormap = [:black, :brown, :white])
                    end
                end
                return channelmaterials, plotfunc
            end
            plotfuncs = last.(cp)
            ps = [pf() for pf in plotfuncs]
        end
        # t[] = lookup(LFP, 𝑡)[300]

        cam3d!(ax)
        cam = cameracontrols(ax) # ax.scene
        cam.eyeposition[] = rot_eyepos(eyeposition, centre, -((frame - 125) / 160)^2 + 0.3)
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

if true
    # @spawnat 2 render_brainframe(200)
    # render_brainframe(125)
    render_brainframe(249)
end
begin
    # pmap(render_brainframe, 1:250; distributed = true)
    # progressmap(render_brainframe, 0:249)
end
