using RPRMakie, FileIO
RPRMakie.activate!(; plugin = RPRMakie.Northstar, iterations = 500)

img = begin
    radiance = 3
    lights = [
        EnvironmentLight(1.5, load(RPR.assetpath("studio026.exr"))),
        PointLight([1, -1, 1.5] |> Vec3f,
                   RGBf(radiance, radiance, radiance * 1.5))
    ]
    fig = Figure(; size = (1080, 720))
    ax = LScene(fig[1, 1]; show_axis = false, scenekw = (lights = lights,))
    screen = RPRMakie.Screen(ax.scene)
    matsys = screen.matsys

    material = RPR.Matx(matsys, "./Glass.mtlx")

    mesh!(ax, load(Makie.assetpath("matball_floor.obj"));
          color = :white)
    matball!(ax, material; color = nothing)
    cam = cameracontrols(ax.scene)
    cam.eyeposition[] = Vec3f(0.0, -1, 0.7)
    cam.lookat[] = Vec3f(0)
    cam.upvector[] = Float32[0.0, -0.01, 1.0]
    update_cam!(ax.scene, cam)
    colorbuffer(screen)
end
save("tempfig.png", img)
