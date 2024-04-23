using GraphMakie
using GraphMakie.Graphs
using SimpleWeightedGraphs
using LinearAlgebra
using JSON

plotdir(args...) = projectdir("plots", args...)
export plotdir

const structures = ["VISp", "VISl", "VISrl", "VISal", "VISpm", "VISam"]
const connector = "-"
const layers = ["1", "2/3", "4", "5", "6"]
# const layercolors = reverse(binarysunset[range(start = 0, stop = 1,
const layercolors = (cgrad(:inferno)[range(start = 0, stop = 0.8,
                                           length = length(layers))])

colors = [Foresight.cornflowerblue,
    Foresight.crimson,
    Foresight.cucumber,
    Foresight.california,
    Foresight.juliapurple]
const structurecolors = [colors..., Makie.RGB(0.5, 0.5, 0.5)]
const structurecolormap = Dict(structures .=> structurecolors)
const lfpcolormap = darksunset
const amplitudecolormap = :bone
const phasecolormap = cyclic

function plotlayerints!(ax, ints; dx = 0.02, width = dx)
    acronyms = "L" .* layers
    ticks = mean.(ints)
    freeze!(ax)
    ax.yticks = (ticks, acronyms)
    ax.yticksvisible = false
    xmins = [!Bool(mod(i, 2)) * dx for i in 1:length(acronyms)]
    xmaxs = xmins .+ width
    hspan!(ax, extrema(vcat(collect.(extrema.(ints))...))...; xmin = 0, xmax = dx * 2,
           color = :white)
    ps = map(ints, layercolors, xmins, xmaxs) do i, color, xmin, xmax
        hspan!(ax, extrema(i)...; xmin, xmax, color)
    end
    return ps
end
function plotlayerints!(ax, ints::DimArray{<:String, 1}; width = 0.02)
    ints = parselayernum.(ints)
    ds = diff(lookup(ints, :depth)) / 2
    append!(ds, first(ds))
    ints = map(unique(ints), ds) do i, d
        x = lookup(ints[ints .== i], :depth)
        return ustrip(minimum(x) - d) .. ustrip(maximum(x) + d)
    end
    return plotlayerints!(ax, ints; width, dx = 0)
end

function wrap(x; domain)
    domain_min, domain_max = extrema(domain)
    domain_range = domain_max - domain_min
    y = mod(x - domain_min, domain_range) + domain_min
    if y < domain_min
        y += domain_range
    end
    return y
end

function plotlayermap!(ax, m, ints; arrows = false,
                       colorrange = maximum(abs.(ustripall(m))) * [-1, 1],
                       colormap = binarysunset, doupsample = true, domain = nothing,
                       rasterize = 10, kwargs...)
    if doupsample
        if isnothing(domain)
            x = upsample(ustripall(m), 5, 2)
        else
            x = mapslices(unwrap, ustripall(m), dims = 2)
            x = upsample(x, 5, 2)
            x = wrap.(x; domain)
        end
    else
        x = ustripall(m)
    end
    p = heatmap!(ax, x; colormap, colorrange, rasterize, kwargs...)
    if arrows
        q = ustripall((m))[1:100:end, 1:3:end]
        q ./= maximum(abs.(q))
        arrows!(ax, lookup(q, 1), lookup(q, 2), zeros(size(q)), parent(q);
                lengthscale = 0.07,
                normalize = false, color = (:black, 0.4))
    end

    ps = vlines!(ax, [0, 0.25]; color = (:white, 0.5), linestyle = :dash, linewidth = 3)
    pl = plotlayerints!(ax, ints)
    return p, ps, pl
end

function tortinset!(gl, ϕ, r; width = Relative(0.15),
                    height = Relative(0.2),
                    halign = 0.05,
                    valign = 0.6, n = 40, colormap = :binary, color = :gray, kwargs...)
    inset_ax = Axis(gl;
                    width,
                    height,
                    halign,
                    valign,
                    xticklabelsvisible = false,
                    yticklabelsvisible = false,
                    limits = ((nothing, nothing), (-1.75, 1.75)),
                    kwargs...)
    hidedecorations!(inset_ax)
    hidespines!(inset_ax)
    inset_ax.scene.backgroundcolor = Makie.RGBA(0, 0, 0, 0)
    # translate!(inset_ax.scene, Vec3f(0, 0, -100000))

    b = ModulationIndices.tortbin(ϕ; n)
    a = r
    h = [mean(a[b .== i]) for i in 1:n]
    angles = range(start = -π + π / n, stop = π - π / n, length = n)
    _, i = findmax(h)
    # lines!(inset_ax, [angles[i], angles[i]], [0, 1], color = :black)
    # if scatter
    #     scatter!(inset_ax, angles, h ./ maximum(h), color = :black)
    # end
    ex = 1.3
    xs = ((-π * ex):0.01:(((π * ex)))) .+ π .+ angles[i]
    cs = collect(xs)

    for i in eachindex(cs)
        cs[i] = h[findmin(abs.(phasegrad.(angles .+ π, xs[i])))[2]]
    end
    ys = -cos.(xs)
    lines!(inset_ax, xs, ys; color = cs, colormap)
    i = findmax(cs)[2]
    scatter!(inset_ax, [xs[i]], [ys[i]]; color, markersize = 15,
             strokecolor = :white, strokewidth = 2)
end

function plotstructurecenters!(ax,
                               ag = SimpleWeightedDiGraph(zeros(length(structures),
                                                                length(structures)));
                               colormap = getindex.([structurecolormap], structures),
                               structures = structures, arrow_size = 10,
                               arrow_shift = :end, curve_distance = 20,
                               curve_distance_usage = true,
                               kwargs...)
    function offset(x, y, p)
        d = norm(x - y)
        x = (x + y) ./ 2
        u = reverse(x - y) |> collect
        u[1] = -u[1]
        u = u ./ norm(u)
        x = x .+ u .* p .* d
    end

    layout = Point2f[(350, 350), # VISp
                     (170, 310), # VISl
                     (300, 130), # VISrl
                     (180, 195), # VISal
                     (475, 240), # VISpm
                     (450, 140)] # VISam
    fwaypoints = Dict(1 => [offset(layout[2], layout[3], 0.4)],
                      2 => [offset(layout[4], layout[3], 0.1)],
                      6 => [offset(layout[4], layout[6], 0.3)],
                      3 => [offset(layout[1], layout[5], -0.1)],
                      4 => [offset(layout[4], layout[5], -0.1)])

    if !(ag isa Observable)
        ag = Observable(ag)
    end
    if !isempty(edges(ag[]))
        ew = lift(ag -> weight.(edges(ag)), ag)
        aw = lift(ew -> ew .* arrow_size, ew)
    else
        ew = 10
        aw = arrow_size .* ew
    end
    p = graphplot!(ax, ag;
                   edge_width = ew,
                   selfedge_size = 0,
                   layout,
                   edge_plottype = :beziersegments,
                   curve_distance,
                   curve_distance_usage,
                   node_size = 30,
                   ilabels = structures,
                   ilabels_fontsize = 9,
                   node_color = colormap,
                   # waypoints=fwaypoints,
                   arrow_size = aw,
                   arrow_shift,
                   arrow_color = (colorant"#0072BD", 0.7),
                   edge_color = (colorant"#0072BD", 0.7), kwargs...)
end

const visual_cortex_file = datadir("visual_cortex.json")
const _visdraworder = ["VISpm", "VISam", "VISrl", "VISal", "VISl", "VISp"]
const _vishierarchyorder = ["VISp", "VISl", "VISrl", "VISal", "VISpm", "VISam"]
function load_visual_cortex(file = visual_cortex_file)
    D = JSON.parsefile(file)
    outline = D["outline"]
    fill = D["fill"]
    fill = map(collect(fill)) do (k, v)
        v = [Float32.(_v) for _v in v]
        v = Makie.Polygon(Point{2}.(v))
        k => v
    end |> Dict
    fill = Dict(_visdraworder .=> getindex.([fill], _visdraworder))
    outline = Dict(k => Vec2f.(v) for (k, v) in collect(outline))
    outline = Dict(_visdraworder .=> getindex.([outline], _visdraworder))
    return fill, outline
end

function plot_visual_cortex!(ax = Axis(Figure()[1, 1], aspect = 1, yreversed = true),
                             showoutline = true,
                             strokewidth = 0,
                             draworder = _visdraworder,
                             colororder = _vishierarchyorder,
                             legendorder = colororder,
                             colormap = cgrad(:inferno, length(legendorder),
                                              categorical = true),
                             alpha = 0.5,
                             kwargs...)
    fill, outline = load_visual_cortex()
    colormap = cgrad(colormap, length(legendorder), categorical = true)
    _colormap = cgrad(colormap, length(legendorder); categorical = true, alpha)
    k = keys(fill) |> collect
    legendorder = indexin(legendorder, k)

    colors = indexin(k, colororder)
    colors = _colormap[colors] # To original ordering of structures
    v = values(outline) |> collect
    p1 = map(legendorder) do i
        poly!(v[i]; label = k[i], strokewidth, color = colors[i], kwargs...)
    end
    if showoutline
        draworder = indexin(draworder, k)
        colors = indexin(k, colororder)
        colors = colormap[colors]
        p2 = map(draworder) do i
            lines!(v[i]; linewidth = 5, color = colors[i])
        end
    else
        p2 = []
    end
    axislegend(ax)
    return p1, p2
end
function plot_visual_cortex(; kwargs...)
    f = Figure()
    ax = Axis(f[1, 1], aspect = 1, yreversed = true)
    p1, p2 = plot_visual_cortex!(ax; kwargs...)
    display(f)
    return f, ax, [p1, p2]
end
