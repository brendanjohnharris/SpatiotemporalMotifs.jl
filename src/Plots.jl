
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
