#! /bin/bash
# -*- mode: julia -*-
#=
exec julia +1.10.10 -t auto --color=yes "${BASH_SOURCE[0]}" "$@"
=#
using DrWatson
@quickactivate "SpatiotemporalMotifs"

using SpatiotemporalMotifs
import TimeseriesTools: freqs
using Peaks
using FileIO
using SpatiotemporalMotifs
import SpatiotemporalMotifs.layers
import SpatiotemporalMotifs.PTHR
using Random
using Distributed
@preamble
set_theme!(foresight(:physics))

stimuli = [r"Natural_Images", "spontaneous", "flash_250ms"]
mkpath(plotdir("fano_factor"))

plot_data, data_file = produce_or_load(Dict(), calcdir("plots");
                                       filename = savepath("fano_factor")) do _
    session_table = load(calcdir("plots", "posthoc_session_table.jld2"), "session_table")
    oursessions = session_table.ecephys_session_id
    path = calcdir("madev")
    QQ = calcquality(path)
    plot_data = map(stimuli) do stimulus
        Q = QQ[stimulus = At(stimulus),
               Structure = At(structures),
               SessionID(At(oursessions))]
        @assert mean(Q) > 0.9

        begin # * Load data
            unitdepths = map(lookup(Q, Structure)) do structure
                out = map(lookup(Q, SessionID)) do sessionid
                    if Q[SessionID = At(sessionid), Structure = At(structure)] == 0
                        return nothing
                    end
                    filename = savepath((@strdict sessionid structure stimulus), "jld2",
                                        path)
                    jldopen(filename, "r") do f
                        unitdepths = f["unitdepths"]
                        mad = f["mad"]
                        layermap = mad.metadata[:layerinfo]
                        layermap = ToolsArray(layermap[3], (Depth(layermap[2]),))

                        unitlayers = map(unitdepths.probedepth) do depth
                            layermap[Depth = Near(depth)]
                        end
                        unitdepths.layer = unitlayers
                        return unitdepths
                    end
                end
                out = filter(!isnothing, out)
                out = filter(!isempty, out)
                # out = filter(out) do x # Remove sessions that don't have data to a reasonable depth
                #     maximum(DimensionalData.metadata(x[1])[:streamlinedepths]) > 0.90
                # end
            end
        end

        plot_data = @strdict unitdepths
        return plot_data
    end

    return Dict(string.(stimuli) .=> plot_data)
end

begin
    fanorange = 10^(1.5) .. 1e3
    stimulus = "spontaneous"
    structure = "VISp"
end

begin # * Plot fano factor for VISp, layer 2/3, spontaneous
    unitdepths = plot_data[stimulus]["unitdepths"][findfirst(structures .== structure)]
    unitdepths = map(unitdepths) do units
        filter(units) do unit
            unit.layer == 2 # Choose layer 2/3
        end
    end
    unitdepths = filter(!isempty, unitdepths)
    oursessions = map(unitdepths) do u
        only(unique(u.ecephys_session_id))
    end
    # end
    # begin # * Plot mean fano factor for 1 subject
    f = TwoPanel()
    ax = Axis(f[1, 1];
              xlabel = "Time lag (ms)",
              ylabel = "Fano Factor",
              title = "Fano Factor VISp L2/3",
              xscale = log10, yscale = log10)
    fanos = map(unitdepths) do units
        fanos = units.fano_factor
        fanos = ToolsArray(fanos, (Unit(units.ecephys_unit_id),)) |> stack
        mean_fano = mean(fanos, dims = Unit)
        mean_fano = dropdims(mean_fano, dims = Unit)
    end
    fanos = ToolsArray(fanos, (SessionID(oursessions),)) |> stack
    fanos = fanos[𝑡 = 1 .. 1000]
    mu = median(fanos, dims = SessionID)
    mu = dropdims(mu, dims = SessionID)
    # s = std(fanos, dims = SessionID)
    sl = mapslices(x -> quantile(x, 0.25), fanos; dims = SessionID)
    sl = dropdims(sl, dims = SessionID)
    su = mapslices(x -> quantile(x, 0.75), fanos; dims = SessionID)
    su = dropdims(su, dims = SessionID)

    fano_range = mu[fanorange] .|> log10
    t_range = times(fano_range) .|> log10
    xx = hcat(ones(length(t_range)), t_range)
    m = xx \ fano_range

    # for x in fanos
    #     lines!(ax, x, color = cornflowerblue, alpha = 0.3, linewidth = 1)
    # end
    band!(ax, times(mu), sl, su, color = cornflowerblue, alpha = 0.3)

    lines!(ax, mu, color = cornflowerblue, linewidth = 3, label = "Mean Fano Factor")

    lines!(ax, exp10.(t_range), exp10.(m[1] .+ m[2] .* t_range);
           color = :red, linestyle = :dash, linewidth = 3)

    text!(ax, [10^0.5], [10^0.4]; text = "c = $(round(m[2], digits=2))",
          align = (:left, :center), fontsize = 20)
    display(f)
end

begin # * Now do slope variation over layers boxplot.
    unitdepths = plot_data[stimulus]["unitdepths"][findfirst(structures .== structure)]
    layer_slopes = map(Dim{:layer}(2:5)) do l# L2/3 to L6
        layer_slopes = map(unitdepths) do units
            units = subset(units, :layer => ByRow(==(l)))
            if isempty(units)
                return nothing
            else
                slopes = map(units.fano_factor) do fano
                    fano_range = fano[fanorange] .|> log10
                    t_range = times(fano_range) .|> log10
                    X = hcat(ones(length(t_range)), t_range)
                    m = X \ fano_range

                    return last(m) # Fano slope
                end

                return ToolsArray(slopes, (Unit(units.ecephys_unit_id),)) |> stack
            end
        end
        idxs = findall(!isnothing, layer_slopes)
        layer_slopes = layer_slopes[idxs]
        oursessions = map(unitdepths) do ls
            only(unique(ls.ecephys_session_id))
        end
        oursessions = oursessions[idxs]
        layer_slopes = mean.(layer_slopes, dims = Unit)
        layer_slopes = dropdims.(layer_slopes, dims = Unit)
        layer_slopes = ToolsArray(layer_slopes, (SessionID(oursessions),)) |> stack
    end
    layer_slopes = ToolsArray(layer_slopes, (Dim{:layer}(2:5),))
end
begin # * plot
    begin
        ax = Axis(f[1, 2]; ylabel = "Cortical layer",
                  xlabel = "Median fano factor exponent",
                  title = "Fano factor exponent", yreversed = true,
                  yticks = (2:5, ["L2/3", "L4", "L5", "L6"]))

        _b = layer_slopes
    end

    begin # * Do statistical test against value of 0.5
        𝑝 = map(_b) do b
            𝑝 = HypothesisTests.SignedRankTest(collect(b))
            𝑝 = HypothesisTests.pvalue(𝑝, tail = :right)
        end
        𝑝[:] .= MultipleTesting.adjust(collect(𝑝)[:], BenjaminiHochberg())
    end
    for l in lookup(_b, :layer)
        if l > 1
            x = _b[layer = At(l)]
            p = 𝑝[layer = At(l)]
            _x = filter(!isnan, collect(x)[:])
            if p < SpatiotemporalMotifs.PTHR
                boxplot!(ax, fill(l, length(_x)), _x;
                         orientation = :horizontal, show_outliers = false,
                         color = Foresight.colororder[l - 1],
                         strokecolor = Foresight.colororder[l - 1], strokewidth = 3)
            else
                boxplot!(ax, fill(l, length(_x)), _x;
                         orientation = :horizontal, show_outliers = false,
                         color = :transparent,
                         strokecolor = Foresight.colororder[l - 1], strokewidth = 3)
            end
        end
    end
    display(f)
    wsave(plotdir("fano_factors.pdf"), f)
end
