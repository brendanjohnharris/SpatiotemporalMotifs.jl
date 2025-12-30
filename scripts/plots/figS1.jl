#! /bin/bash
# -*- mode: julia -*-
#=
exec julia +1.10.10 -t auto --color=yes "${BASH_SOURCE[0]}" "$@"
=#
using DrWatson
@quickactivate "SpatiotemporalMotifs"

using SpatiotemporalMotifs
import TimeseriesTools: freqs
using SpatiotemporalMotifs
using Distributed

using Random
@preamble
set_theme!(foresight(:physics))

config = (; shuffles = 1000,
          regcoef = 0.1,
          repeats = 5,
          folds = 20,
          nsur = 10)
layerints = load(calcdir("plots", "grand_unified_layers.jld2"), "layerints")

if nprocs() == 1 && !isfile(calcdir("plots", savepath("figS1", config, "jld2"))) # * Initialize workers
    addprocs(10)
    @everywhere using DrWatson
    @everywhere DrWatson.@quickactivate "SpatiotemporalMotifs"
    @everywhere using SpatiotemporalMotifs
    @everywhere @preamble
end

plot_data, data_file = produce_or_load(config, calcdir("plots");
                                       filename = savepath("figS1")) do config
    @unpack shuffles, regcoef, repeats, folds, nsur = config
    begin # * Load the trial LFP for natural images
        session_table = load(calcdir("plots", "posthoc_session_table.jld2"),
                             "session_table")
        oursessions = session_table.ecephys_session_id
        path = calcdir("calculations")
        Q = calcquality(path)[Structure = At(structures)]
        Q = Q[SessionID(At(oursessions))]
        @assert mean(Q[stimulus = At(r"Natural_Images")]) == 1
        out = load_calculations(Q; stimulus = r"Natural_Images", vars = [:V])
    end

    begin # * 1/f across all sessions. Should take about 15 minutes over 50 workers
        oneoneff = map(out) do out_structure
            @info "Calculating 1/f exponent for $(metadata(out_structure[1][:V])[:structure])"
            pmap(out_structure) do out_session
                V = out_session[:V]
                V = V[𝑡 = SpatiotemporalMotifs.INTERVAL]
                map(SpatiotemporalMotifs.trial_1onf,
                    eachslice(V, dims = (:Depth, :changetime)))
            end
        end
    end
    begin
        χ = map(oneoneff) do oneoneff_structure
            x = map(oneoneff_structure) do oneoneff_session
                getindex.(oneoneff_session, :χ)
            end
            ToolsArray(x, (SessionID(oursessions),))
        end
        χ = ToolsArray(χ, (Structure(structures),)) |> stack
        b = map(oneoneff) do oneoneff_structure
            x = map(oneoneff_structure) do oneoneff_session
                getindex.(oneoneff_session, :b)
            end
            ToolsArray(x, (SessionID(oursessions),))
        end
        b = ToolsArray(b, (Structure(structures),)) |> stack
        hitmiss = map(out) do structure
            x = map(structure) do session
                hitmiss = session[:trials].hit
            end
            ToolsArray(x, (SessionID(oursessions),))
        end
        hitmiss = ToolsArray(hitmiss, (Structure(structures),)) |> stack
    end

    begin # * Bin to common depths
        bins = intervals(0.0:0.1:0.9)

        Δχ = map(χ, hitmiss) do χ, hm
            local a = χ[:, hm]
            local b = χ[:, .!hm]
            d = median(a, dims = 2) .- median(b, dims = 2)
            y = dropdims(d, dims = 2)
            y = groupby(y, Depth => Bins(bins))
            y = mean.(y)
            set(y, Depth => mean.(lookup(y, Depth)))
        end |> stack

        Δχs = pmap(χ, hitmiss) do χ, hm
            z = map(1:shuffles) do _
                # * Shuffle the hit/miss labels
                idxs = randperm(length(hm))
                hm = hm[idxs]
                local a = χ[:, hm]
                local b = χ[:, .!hm]
                d = median(a, dims = 2) .- median(b, dims = 2)
                y = dropdims(d, dims = 2)
                y = groupby(y, Depth => Bins(bins))
                y = mean.(y)
                set(y, Depth => mean.(lookup(y, Depth)))
            end
            ToolsArray(z, (Dim{:repeats}(1:shuffles),)) |> stack
        end |> stack

        𝑝 = map(eachslice(Δχ, dims = (Depth, Structure)),
                eachslice(Δχs, dims = (Depth, Structure))) do x, y
            # * Perform the Mann-Whitney U test
            x = median(x)
            y = median(y, dims = SessionID)
            n_extreme = count(abs.(y) .>= abs(x))
            p = (n_extreme + 1) / (length(y) + 1)
            return p
        end
        𝑝[:] .= MultipleTesting.adjust(𝑝[:], BenjaminiHochberg())
    end

    begin # * Try and do a trial-by-trial classification using 1/f exponents
        seshχ = map(χ, hitmiss) do x, hm
            set(x, :changetime => Trial(hm))
        end

        # * Mean 1/f classifier
        bac_mean = pmap(eachslice(seshχ, dims = SessionID)) do x
            x = mean.(x, dims = Depth)
            x = ToolsArray(x, (Structure(structures),)) |> stack
            x = dropdims(x, dims = Depth)

            if mean(lookup(x, Trial)) < 0.9 # * Very low number of miss trials
                bac = classify_kfold(x; regcoef, k = folds, repeats = repeats,
                                     dim = Trial)

                bac_sur = map(1:nsur) do _
                    x_sur = set(x, Trial => lookup(x, Trial)[randperm(size(x, Trial))])
                    classify_kfold(x_sur; regcoef, k = folds, repeats = repeats,
                                   dim = Trial)
                end
            else
                return nothing
            end
            return bac, bac_sur
        end

        # * Layer-wise 1/f classifier
        bac = pmap(eachslice(seshχ, dims = SessionID)) do x
            hm = lookup(x |> first, Trial)
            @assert all(only(unique(lookup.(x, Trial))) == hm)
            x = vcat(parent.(x)...)
            x = ToolsArray(x, (Var(axes(x, 1)), Trial(hm)))

            if mean(lookup(x, Trial)) < 0.9 # * Very low number of miss trials
                bac = classify_kfold(x; regcoef, k = folds, repeats = repeats,
                                     dim = Trial)

                bac_sur = map(1:nsur) do _
                    x_sur = set(x, Trial => lookup(x, Trial)[randperm(size(x, Trial))])
                    classify_kfold(x_sur; regcoef, k = folds, repeats = repeats,
                                   dim = Trial)
                end
                return bac, bac_sur
            else
                return nothing
            end
        end

        W = pmap(eachslice(seshχ, dims = SessionID)) do x
            _x = x
            ds = intervals(0.0:0.1:0.9)
            x = map(x) do x
                    x = groupby(x, Depth => Bins(ds))
                    x = map(x) do x
                        if isempty(x)
                            return nothing
                        else
                            dropdims(mean(x, dims = Depth), dims = Depth)
                        end
                    end
                    y = first(x[.!isnothing.(x)])
                    y .= rand(length(y))# Dummy random values
                    x[isnothing.(x)] .= [y]
                    return x
                end |> stack |> stack
            x = set(x, Depth => mean.(lookup(x, Depth)))
            x = permutedims(x, (Structure, Depth, Trial))

            hm = lookup(x, Trial)
            X = hcat(collect.(Iterators.flatten.(eachslice(x, dims = Trial)))...)
            X = ToolsArray(X, (Var(axes(X, 1)), Trial(hm)))

            N, M = classifier(x; regcoef) # !!!
            W = projection(M)
            W = reshape(W, size(x)[1:2])
            return ToolsArray(W, dims(x)[1:2])
        end |> stack
        W .= W ./ maximum(abs.(W)) # ! Normalized
    end

    D = @strdict χ b hitmiss Δχ Δχs 𝑝 bac_mean bac W regcoef repeats folds nsur
    return D
end

begin # * Load data
    @unpack χ, b, hitmiss, Δχ, Δχs, 𝑝, bac_mean, bac, W = plot_data
    # @assert length(χ) == length(b)
end

begin # * Set up figure
    f = TwoPanel()
    gs = subdivide(f, 1, 4)
    mkpath(datadir("source_data", "figS1")) # ?
end

begin # * Plotting
    # * Filter for oursessions....
    ax = Axis(gs[1];
              ylabel = "Cortical depth (%)",
              xlabel = "1/f exponent (hit - miss)",
              ytickformat = depthticks, title = "1/f contrast", yreversed = true,
              limits = ((-0.07, 0.12), (0, 1)),
              xtickformat = terseticks)

    plotlayerints!(ax, layerints; axis = :y, flipside = true, newticks = false)
    SpatiotemporalMotifs.plot_layerwise!(ax, Δχ; scatter = true)

    # map(eachslice(Δχ, dims = Structure)) do Δχ_structure
    #     structure = only(refdims(Δχ_structure, Structure))
    #     μ, (σl, σh) = bootstrapmedian(Δχ_structure, dims = 2)

    #     _μ = upsample(μ, 5)
    #     _σl, _σh = upsample(σl, 5), upsample(σh, 5)

    #     band!(ax, Point2f.(collect(_σl), lookup(_μ, 1)),
    #           Point2f.(collect(_σh), lookup(_μ, 1)); label = structure,
    #           color = (structurecolormap[structure], 0.2))
    #     lines!(ax, collect(_μ), lookup(_μ, 1), label = structure,
    #            color = structurecolormap[structure])

    #     sigs = 𝑝[Structure = At(structure)] .< SpatiotemporalMotifs.PTHR

    #     scatter!(ax, collect(μ[sigs]), lookup(μ[sigs], 1),
    #              color = structurecolormap[structure], label = structure)

    #     scatter!(ax, collect(μ[.!sigs]), lookup(μ[.!sigs], 1),
    #              strokecolor = structurecolormap[structure], label = structure,
    #              color = :white, strokewidth = 4)
    # end
    display(f)

    begin # ? Save source data for panel a (1/f contrast)
        outdir = datadir("source_data", "figS1") # ?
        df = DataFrame() # ?
        for (i, s) in enumerate(structures) # ?
            sname = replace(s, "/" => "") # ?
            Δχ_s = Δχ[Structure = At(s)] # ?
            μ, (σl, σh) = bootstrapmedian(Δχ_s; dims = 2) # ?
            if i == 1 # ?
                df[!, "Cortical depth (%)"] = collect(lookup(μ, 1)) # ?
            end # ?
            df[!, "Median 1/f contrast $sname"] = collect(μ) # ?
            df[!, "CI lower $sname"] = collect(σl) # ?
            df[!, "CI upper $sname"] = collect(σh) # ?
        end # ?
        CSV.write(joinpath(outdir, "panel_a_1f_contrast.csv"), df) # ?
    end # ?
end

begin # * Plot classification accuracy
    ax = Axis(gs[2], ylabel = "Balanced accuracy",
              limits = ((0.5, 4.5), (nothing, 1.0)),
              xticks = ([1.5, 3.5], ["Regional", "Layer-wise"]),
              title = "1/f classification",
              ytickformat = terseticks)

    hlines!(ax, [0.5]; color = (:black, 0.5), linestyle = :dash, linewidth = 2)

    # * Mean 1/f accuracy
    x = first.(filter(!isnothing, bac_mean))
    xs = last.(filter(!isnothing, bac_mean)) |> Iterators.flatten |> collect
    vspan!(ax, 0.5, 2.5; alpha = 0.2, color = (california, 0.2))
    boxplot!(ax, ones(length(x)), x, color = (california, 0.7), whiskerlinewidth = 3)
    boxplot!(ax, ones(length(xs)) .+ 1, xs, color = (:gray, 0.7), whiskerlinewidth = 3)

    # * Layerwise 1/f accuracy
    x = first.(filter(!isnothing, bac))
    xs = last.(filter(!isnothing, bac)) |> Iterators.flatten |> collect
    vspan!(ax, 2.5, 4.5; color = (cucumber, 0.2))
    boxplot!(ax, ones(length(x)) .+ 2, x, color = (cucumber, 0.7), whiskerlinewidth = 3)
    boxplot!(ax, ones(length(xs)) .+ 3, xs, color = (:gray, 0.7), whiskerlinewidth = 3)

    display(f)

    begin # ? Save source data for panel b (classification accuracy)
        outdir = datadir("source_data", "figS1") # ?
        # ? Regional (mean 1/f) classifier
        regional_acc = first.(filter(!isnothing, bac_mean)) # ?
        regional_sur = last.(filter(!isnothing, bac_mean)) |> Iterators.flatten |> collect # ?
        # ? Layer-wise classifier
        layerwise_acc = first.(filter(!isnothing, bac)) # ?
        layerwise_sur = last.(filter(!isnothing, bac)) |> Iterators.flatten |> collect # ?
        # ? Build long-format DataFrame
        dfs = DataFrame[] # ?
        push!(dfs,
              DataFrame("Classifier" => fill("Regional", length(regional_acc)), # ?
                        "Type" => fill("Actual", length(regional_acc)), # ?
                        "Balanced accuracy" => regional_acc)) # ?
        push!(dfs,
              DataFrame("Classifier" => fill("Regional", length(regional_sur)), # ?
                        "Type" => fill("Surrogate", length(regional_sur)), # ?
                        "Balanced accuracy" => regional_sur)) # ?
        push!(dfs,
              DataFrame("Classifier" => fill("Layer-wise", length(layerwise_acc)), # ?
                        "Type" => fill("Actual", length(layerwise_acc)), # ?
                        "Balanced accuracy" => layerwise_acc)) # ?
        push!(dfs,
              DataFrame("Classifier" => fill("Layer-wise", length(layerwise_sur)), # ?
                        "Type" => fill("Surrogate", length(layerwise_sur)), # ?
                        "Balanced accuracy" => layerwise_sur)) # ?
        CSV.write(joinpath(outdir, "panel_b_classification_accuracy.csv"), vcat(dfs...)) # ?
    end # ?
end

begin # * Plot weights of lda classifier for full classification
    ax = Axis(gs[3], xlabel = "Weight", ylabel = "Cortical depth (%)",
              title = "Classifier weights", yreversed = true,
              limits = ((-0.35, 0.35), (0, 1)),
              xtickformat = terseticks)

    plotlayerints!(ax, layerints; axis = :y, flipside = true, newticks = false)
    SpatiotemporalMotifs.plot_layerwise!(ax, W; scatter = true)

    display(f)

    begin # ? Save source data for panel c (classifier weights)
        outdir = datadir("source_data", "figS1") # ?
        df = DataFrame() # ?
        for (i, s) in enumerate(structures) # ?
            sname = replace(s, "/" => "") # ?
            W_s = W[Structure = At(s)] # ?
            μ, (σl, σh) = bootstrapmedian(W_s; dims = 2) # ?
            if i == 1 # ?
                df[!, "Cortical depth (%)"] = collect(lookup(μ, 1)) # ?
            end # ?
            df[!, "Median weight $sname"] = collect(μ) # ?
            df[!, "CI lower $sname"] = collect(σl) # ?
            df[!, "CI upper $sname"] = collect(σh) # ?
        end # ?
        CSV.write(joinpath(outdir, "panel_c_classifier_weights.csv"), df) # ?
    end # ?
end

begin # * Save figure
    Legend(gs[4], ax, "Area"; merge = true) |> reverselegend!
    addlabels!(f, labelformat)
    wsave(plotdir("figS1", "figS1.pdf"), f)
    display(f)
end
