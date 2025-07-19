using GraphMakie
using GraphMakie.Graphs
using SimpleWeightedGraphs
using LinearAlgebra
using JSON
import TimeseriesTools: freqs

const connector = "&"
const layers = ["1", "2/3", "4", "5", "6"]
# const layercolors = reverse(binarysunset[range(start = 0, stop = 1,
const layercolormap = reverse(cgrad(:roma)) #Foresight.pelagic
const layercolors = (layercolormap[range(start = 0.0, stop = 1,
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
const DEFAULT_SESSION_ID = 1140102579
const DEFAULT_TRIAL_NUM = 14

labelformat(n) = "$(Char(n + 96))"
export labelformat

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
function plotlayerints!(ax, ints; dx = 0.03, width = dx, axis = :y, newticks = true,
                        flipside = false, bgcolor = :transparent)
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
    bmin = minimum([xmins xmaxs])
    bmax = maximum([xmins xmaxs])
    if axis === :y
        hspan!(ax, extrema(vcat(collect.(extrema.(ints))...))...;
               xmin = bmin,
               xmax = bmax,
               color = bgcolor)
        ps = map(ints, layercolors, xmins, xmaxs) do i, color, xmin, xmax
            hspan!(ax, extrema(i)...; xmin, xmax, color)
        end
    elseif axis === :x
        vspan!(ax, extrema(vcat(collect.(extrema.(ints))...))...;
               ymin = bmin,
               ymax = bmax,
               color = bgcolor)
        ps = map(ints, layercolors, xmins, xmaxs) do i, color, xmin, xmax
            vspan!(ax, extrema(i)...; ymin = xmin, ymax = xmax, color)
        end
    end

    return ps
end
function plotlayerints!(ax, ints::ToolsArray{<:String, 1}; width = 0.04, kwargs...)
    ints = parselayernum.(ints)
    ds = diff(lookup(ints, Depth)) / 2
    append!(ds, first(ds))
    ints = map(unique(ints), ds) do i, d
        x = lookup(ints[ints .== i], Depth)
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
                       colormap = defaultcolormap, doupsample = true, domain = nothing,
                       rasterize = 5, stimulus = [0, 0.25], lengthscale = 0.07,
                       arrowsize = 5, arrowcolor = (:black, 0.4), kwargs...)
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
                lengthscale, arrowsize,
                normalize = false, color = arrowcolor)
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
function tortinset!(gl, ϕ; width = Relative(0.15),
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
    h = [sum(b .== i) for i in 1:n]
    angles = range(start = -π + π / n, stop = π - π / n, length = n)
    _, i = findmax(h)
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
                               scale = [-1.8, 2],
                               offset = [32, 0],
                               kwargs...)
    function woffset(x, y, p)
        d = norm(x - y)
        x = (x + y) ./ 2
        u = reverse(x - y) |> collect
        u[1] = -u[1]
        u = u ./ norm(u)
        x = x .+ u .* p .* d
    end

    function polygon_centroid(vertices)
        n = length(vertices)
        area = 0.0
        Cx = 0.0
        Cy = 0.0

        for i in 1:n
            x0, y0 = vertices[i]
            x1, y1 = vertices[mod1(i + 1, n)]
            cross = x0 * y1 - x1 * y0
            area += cross
            Cx += (x0 + x1) * cross
            Cy += (y0 + y1) * cross
        end

        area /= 2.0
        Cx /= (6.0 * area)
        Cy /= (6.0 * area)

        return (Cx, Cy)
    end

    outlines = load_visual_cortex(; scale, offset)[1]
    outlines = getindex.([outlines], structures)
    outlines = Makie.GeometryBasics.coordinates.(outlines)
    layout = outlines .|> polygon_centroid .|>
             Point2f
    fwaypoints = Dict(1 => [woffset(layout[2], layout[3], 0.4)],
                      2 => [woffset(layout[4], layout[3], 0.1)],
                      6 => [woffset(layout[4], layout[6], 0.3)],
                      3 => [woffset(layout[1], layout[5], -0.1)],
                      4 => [woffset(layout[4], layout[5], -0.1)])

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

const visual_cortex_file = projectdir("scripts/plots/visual_cortex/visual_cortex.json")
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
                             strokealpha = 1.0,
                             scale = [-1.8, 2],
                             offset = [32, 0],
                             colors = nothing,
                             kwargs...)
    fill, outline = load_visual_cortex(; scale, offset)
    v = values(fill) |> collect
    k = keys(fill) |> collect
    legendorder = indexin(legendorder, k)
    if isnothing(colors)
        _colormap = cgrad(colormap, length(legendorder); categorical = true,
                          alpha = fillalpha)
        colormap = brighten.(colormap, 1 - strokealpha)
        colormap = cgrad(colormap, length(legendorder), categorical = true)
        colors = indexin(k, colororder)
        colors = _colormap[colors] # To original ordering of structures
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
    else
        colors = lift(c -> c[indexin(k, colororder)], colors)
        p1 = map(legendorder) do i
            poly!(v[i]; label = k[i], strokewidth,
                  color = lift(colors -> colors[i], colors),
                  alpha = fillalpha,
                  colormap, kwargs...)
        end
        if showoutline
            draworder = indexin(draworder, k)
            p2 = map(draworder) do i
                lines!(v[i]; linewidth = 5, color = lift(colors -> colors[i], colors),
                       colormap, kwargs...)
            end
        else
            p2 = []
        end
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
function fooof(x; kwargs...)
    AN.aperiodicfit(x, [3, 300]; aperiodic_mode = "fixed", max_n_peaks = 8,
                    peak_threshold = 1, peak_width_limits = [1, 50], kwargs...)
end
function trial_1onf(x::UnivariateRegular)
    minfreq = 10.0u"Hz"
    maxfreq = 300.0u"Hz" # Should match upper limit of the 'pass' band in Calculations.jl
    S = spectrum(x, 5u"Hz")
    S = S[𝑓 = minfreq .. maxfreq] |> ustripall
    f, D = fooof(S; peak_width_limits = [10, 100], max_n_peaks = 2)

    # lines(freqs(S), S |> parent; axis = (; xscale = log10, yscale = log10),
    #       label = "spectrum")
    # lines!(freqs(S), f.(freqs(S)))
    # current_figure() |> display
    return D
end
function plotspectrum!(ax, s::AbstractToolsArray;
                       textposition = (14, exp10(-2.9)), annotations = [:peaks, :fooof],
                       color = cucumber, label = nothing)
    μ = mean(s, dims = (SessionID, :layer))
    μ = dropdims(μ, dims = (SessionID, :layer)) |> ustripall
    σ = std(s, dims = (SessionID, :layer)) ./ 2
    σ = dropdims(σ, dims = (SessionID, :layer)) |> ustripall
    p = lines!(ax, freqs(μ), collect(μ); color = (color, 0.8), label)
    band!(ax, freqs(μ), collect.([max.(μ - σ, eps()), μ + σ])...; color = (color, 0.32),
          label)

    # * Find peaks
    if :peaks in annotations
        pks, proms = findpeaks(μ, 2; N = 2)
        scatter!(ax, collect(freqs(pks)), collect(pks .* 1.25), color = :black,
                 markersize = 10, marker = :dtriangle)
        text!(ax, collect(freqs(pks)), collect(pks);
              text = string.(round.(collect(freqs(pks)), digits = 1)) .* [" Hz"],
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
