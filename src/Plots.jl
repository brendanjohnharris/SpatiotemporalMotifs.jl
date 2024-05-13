using GraphMakie
using GraphMakie.Graphs
using SimpleWeightedGraphs
using LinearAlgebra
using JSON
import TimeseriesTools: freqs

plotdir(args...) = projectdir("plots", args...)
export plotdir

const connector = "&"
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
const structurecgrad = cgrad(structurecolors; categorical = true)
const defaultcolormap = binarysunset
const lfpcolormap = darksunset
const amplitudecolormap = :bone
const phasecolormap = cyclic
DEFAULT_SESSION_ID = 1140102579
DEFAULT_TRIAL_NUM = 26

const visual_cortex_layout = Dict("VISp" => [350, 350],
                                  "VISl" => [170, 310],
                                  "VISrl" => [300, 130],
                                  "VISal" => [180, 195],
                                  "VISpm" => [475, 240],
                                  "VISam" => [450, 140])
const hierarchy_scores = Dict("VISp" => -0.357, "VISl" => -0.093, "VISrl" => -0.059,
                              "VISal" => 0.152, "VISpm" => 0.327, "VISam" => 0.441) # * Anatomical hierarchy scores from Siegle 2021
function depthticks(x)
    if all(x .< 50) # Assume percentages
        x = x .* 100
    end
    return string.(round.(Int, x))
end
function plotlayerints!(ax, ints; dx = 0.02, width = dx, axis = :y, newticks = true,
                        flipside = false)
    acronyms = "L" .* layers
    ticks = mean.(ints)
    # freeze!(ax)
    if axis === :y && newticks
        ax.yticks = (ticks, acronyms)
        ax.yticksvisible = false
    elseif axis === :x && newticks
        ax.xticks = (ticks, acronyms)
        ax.xticksvisible = false
    end
    xmins = [!Bool(mod(i, 2)) * dx for i in 1:length(acronyms)]
    xmaxs = xmins .+ width
    if flipside
        xmins = 1 .- xmins
        xmaxs = 1 .- xmaxs
    end
    if axis === :y
        hspan!(ax, extrema(vcat(collect.(extrema.(ints))...))...;
               xmin = minimum([xmins xmaxs]),
               xmax = maximum([xmins xmaxs]),
               color = :white)
        ps = map(ints, layercolors, xmins, xmaxs) do i, color, xmin, xmax
            hspan!(ax, extrema(i)...; xmin, xmax, color)
        end
    elseif axis === :x
        vspan!(ax, extrema(vcat(collect.(extrema.(ints))...))...;
               ymin = minimum([xmins xmaxs]),
               ymax = maximum([xmins xmaxs]),
               color = :white)
        ps = map(ints, layercolors, xmins, xmaxs) do i, color, xmin, xmax
            vspan!(ax, extrema(i)...; ymin = xmin, ymax = xmax, color)
        end
    end

    return ps
end
function plotlayerints!(ax, ints::DimArray{<:String, 1}; width = 0.02, kwargs...)
    ints = parselayernum.(ints)
    ds = diff(lookup(ints, :depth)) / 2
    append!(ds, first(ds))
    ints = map(unique(ints), ds) do i, d
        x = lookup(ints[ints .== i], :depth)
        return ustrip(minimum(x) - d) .. ustrip(maximum(x) + d)
    end
    return plotlayerints!(ax, ints; width, dx = 0, kwargs...)
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
function plotlayermap!(ax, m; arrows = false,
                       colorrange = extrema(ustripall(m)),
                       colormap = binarysunset, doupsample = true, domain = nothing,
                       rasterize = 10, stimulus = [0, 0.25], kwargs...)
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
    if arrows == true
        arrows = (100, 3)
    end
    if length(arrows) == 1
        arrows = (arrows, 3)
    end
    if first(arrows) != 0
        q = ustripall((m))[1:arrows[1]:end, 1:arrows[2]:end]
        q ./= maximum(abs.(q))
        arrows!(ax, lookup(q, 1), lookup(q, 2), zeros(size(q)), parent(q);
                lengthscale = 0.07, arrowsize = 0.2,
                normalize = false, color = (:black, 0.4))
    end
    if !isempty(stimulus) && !isnothing(stimulus)
        ps = vlines!(ax, stimulus; color = (:white, 0.5), linestyle = :dash, linewidth = 3)
    else
        ps = []
    end
    return p, ps
end

function plotlayermap!(ax, m, ints; kwargs...)
    p, ps = plotlayermap!(ax, m; kwargs...)
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
function load_visual_cortex(file = visual_cortex_file; scale = [1.8, 2], offset = [-32, 0])
    D = JSON.parsefile(file)
    outline = D["outline"]
    fill = D["fill"]
    fill = map(collect(fill)) do (k, v)
        v = [Float32.((_v .+ offset) .* scale) for _v in v]
        v = Makie.Polygon(Point{2}.(v))
        k => v
    end |> Dict
    fill = Dict(_visdraworder .=> getindex.([fill], _visdraworder))
    outline = map(collect(outline)) do (k, v)
        k => map(v -> Vec2f.((v .+ offset) .* scale), v)
    end |> Dict
    outline = Dict(_visdraworder .=> getindex.([outline], _visdraworder))
    return fill, outline
end

function plot_visual_cortex!(ax = Axis(Figure()[1, 1], aspect = 1, yreversed = true);
                             showoutline = true,
                             strokewidth = 0,
                             draworder = _visdraworder,
                             colororder = _vishierarchyorder,
                             legendorder = colororder,
                             colormap = cgrad(:inferno, length(legendorder),
                                              categorical = true),
                             fillalpha = 0.5,
                             strokealpha = 1,
                             scale = [1.8, 2],
                             offset = [32, 0],
                             kwargs...)
    fill, outline = load_visual_cortex(; scale, offset)
    _colormap = cgrad(colormap, length(legendorder); categorical = true, alpha = fillalpha)
    colormap = brighten.(colormap, 1 - strokealpha)
    colormap = cgrad(colormap, length(legendorder), categorical = true)
    k = keys(fill) |> collect
    legendorder = indexin(legendorder, k)

    colors = indexin(k, colororder)
    colors = _colormap[colors] # To original ordering of structures
    v = values(fill) |> collect
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
    l = axislegend(ax, position = :lt)
    return p1, p2, l
end
function plot_visual_cortex(; kwargs...)
    f = Figure()
    ax = Axis(f[1, 1], aspect = 1, yreversed = true)
    p1, p2 = plot_visual_cortex!(ax; kwargs...)
    display(f)
    return f, ax, [p1, p2]
end
fooof = x -> AN.aperiodicfit(x, [3, 300]; aperiodic_mode = "fixed", max_n_peaks = 5,
                             peak_threshold = 1, peak_width_limits = [1, 50])
function plotspectrum!(ax, s::AbstractDimArray;
                       textposition = (14, exp10(-2.9)), annotations = [:peaks, :fooof],
                       color = cucumber, label = nothing)
    μ = mean(s, dims = (:sessionid, :layer))
    μ = dropdims(μ, dims = (:sessionid, :layer)) |> ustripall
    σ = std(s, dims = (:sessionid, :layer)) ./ 2
    σ = dropdims(σ, dims = (:sessionid, :layer)) |> ustripall
    p = lines!(ax, freqs(μ), μ; color = (color, 0.8), label)
    band!(ax, freqs(μ), collect.([max.(μ - σ, eps()), μ + σ])...; color = (color, 0.32))

    # * Find peaks
    if :peaks in annotations
        pks, proms = findpeaks(μ, 2; N = 2)
        scatter!(ax, freqs(pks), pks .* 1.25, color = :black,
                 markersize = 10, marker = :dtriangle)
        text!(ax, freqs(pks), pks;
              text = string.(round.(freqs(pks), digits = 1)) .* [" Hz"],
              align = (:center, :bottom), color = :black, rotation = 0,
              fontsize = 16,
              offset = (0, 15))
    end

    # * Fooof fit
    ff, ps = fooof(μ)
    lines!(ax, freqs(μ), ff.(freqs(μ)), color = (color, 0.5), linestyle = :dash,
           linewidth = 5)
    if :fooof in annotations
        text!(ax, textposition...; text = L"𝛂 = $(round(ps[:χ], sigdigits=3))",
              fontsize = 16, align = (:right, :center))
    end
    return p, ps[:χ]
end
