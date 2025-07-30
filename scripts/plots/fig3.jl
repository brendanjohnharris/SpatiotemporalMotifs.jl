#! /bin/bash
#=
exec julia +1.10.10 -t auto "${BASH_SOURCE[0]}" "$@"
=#
# ? Expected execution time: 30 mins
using DrWatson
@quickactivate "SpatiotemporalMotifs"

import TimeseriesTools: freqs
using FileIO
using Unitful
import DimensionalData: metadata
using MultivariateStats
using SpatiotemporalMotifs
using USydClusters
@preamble
set_theme!(foresight(:physics))
Random.seed!(32)

begin # * Parameters
    INTERVAL = SpatiotemporalMotifs.INTERVAL
    mainstructure = "VISl"
    maincolorrange = [-2.9, 2.9]
    csdcolorrange = [-1, 1]
    ylims = (0.05, 0.95)
    datafile = calcdir("plots", "theta_order_parameter.jld2")
    hyperfile = calcdir("plots", "hyperparameters", "theta_waves_task.jld2")

    regcoef = 0.5 # Set approximately after inspecting the hyperparameter search
    folds = 5
    repeats = 20
    path = "calculations"

    config = @strdict regcoef folds repeats path
end

if (isfile(hyperfile) && isfile(datafile)) || isfile(savepath("fig3", config, "jld2"))
else # * Use extra workers if we can
    if nprocs() == 1
        if SpatiotemporalMotifs.CLUSTER()
            ourprocs = USydClusters.Physics.addprocs(10; mem = 20, ncpus = 8,
                                                     project = projectdir(),
                                                     queue = "l40s")
        else
            addprocs(7)
        end
    end
    @everywhere using SpatiotemporalMotifs
    @everywhere SpatiotemporalMotifs.@preamble
end

plot_data, data_file = produce_or_load(config, calcdir("plots");
                                       filename = savepath("fig3")) do config
    @unpack regcoef, folds, repeats, path = config
    session_table = load(calcdir("posthoc_session_table.jld2"), "session_table")
    oursessions = session_table.ecephys_session_id
    # vars = [:csd, :k, :ω]
    vars = [:θ, :csd, :k, :ω]
    path = calcdir(path)
    mkpath(plotdir("fig3"))
    statsfile = plotdir("fig3", "$mainstructure.txt")
    close(open(statsfile, "w"))

    begin # * First, do the order parameters
        begin # * Flashes
            stimulus = "flash_250ms"
            Q = calcquality(path)[Structure = At(structures)]
            quality = mean(Q[stimulus = At(stimulus)])
            out = load_calculations(Q; stimulus, vars, path)

            begin # * Calculate a global order parameter and mean LFP at each time point
                out = map(out) do o
                    filter(o) do _o
                        _o[:sessionid] ∈ oursessions
                    end
                end
                Og = map(out) do o
                    O = map(o) do p
                        k = p[:k]
                        ω = p[:ω]
                        k[ustripall(ω) .< 0] .= NaN * unit(eltype(k))
                        ret = dropdims(nansafe(mean; dims = Depth)(sign.(k)), dims = Depth) # Ignore negative frequencies
                        set(ret, Dim{:changetime} => Trial) .|> Float32
                    end
                end
            end
            flashes = @strdict Og
        end

        begin # * Passive trials?
            stimulus = "Natural_Images_passive"
            Q = calcquality(path; require = ["θ"])[Structure = At(structures)]
            quality = mean(Q[stimulus = At(stimulus)])
            out = load_calculations(Q; stimulus, vars, path)

            session_table = load(calcdir("posthoc_session_table.jld2"), "session_table")
            oursessions = session_table.ecephys_session_id

            begin # * Calculate a global order parameter and mean LFP at each time point
                out = map(out) do o
                    filter(o) do _o
                        _o[:sessionid] ∈ oursessions
                    end
                end
                Og = map(out) do o
                    O = map(o) do p
                        k = p[:k]
                        ω = p[:ω]
                        k[ustripall(ω) .< 0] .= NaN * unit(eltype(k))
                        ret = dropdims(nansafe(mean; dims = Depth)(sign.(k)), dims = Depth) # Ignore negative frequencies
                        set(ret, Dim{:changetime} => Trial) .|> Float32
                    end
                end
            end
            passive = @strdict Og
        end
        begin # * Natural images
            stimulus = r"Natural_Images"
            Q = calcquality(path)[Structure = At(structures)]
            quality = mean(Q[stimulus = At(stimulus)])
            out = load_calculations(Q; stimulus, vars, path)

            map(out) do o # * Remove trials with a reaction time < 0.25s (so we can reliably divide trial periods into 'pre-reaction' and 'post-reaction' based on the timestamps. There are typically only a few of these per session.
                filter!(o) do _o
                    begin
                        reaction_times = _o[:trials].lick_latency
                        idxs = (reaction_times .> 0.25) .| isnan.(reaction_times) # NaN means a miss trial

                        @assert issorted(lookup(_o[:k], :changetime))
                        _o[:θ] = _o[:θ][:, :, idxs] # Remove trials with low reaction time
                        _o[:k] = _o[:k][:, :, idxs]
                        _o[:ω] = _o[:ω][:, :, idxs]
                        _o[:trials] = _o[:trials][idxs, :]
                    end
                    _o[:sessionid] ∈ oursessions
                end
            end
            begin
                Og = map(out) do o
                    O = map(o) do p
                        k = p[:k]
                        ω = p[:ω]
                        k[ustripall(ω) .< 0] .= NaN * unit(eltype(k))
                        # idxs = parselayernum.(metadata(k)[:layernames]) .∈ [[2, 3, 4]] # Just select middle layers
                        # k = k[:, parent(idxs), :]
                        ret = dropdims(nansafe(mean, dims = Depth)(sign.(k)), dims = Depth)
                        trials = p[:trials][1:size(k, :changetime), :]
                        trialtimes = trials.change_time_with_display_delay
                        @assert all(isapprox.(ustripall(lookup(k, :changetime)), trialtimes,
                                              atol = 0.05)) # These won't match exactly, because the data change times have been adjusted for rectification. This is OK.
                        set(ret, Dim{:changetime} => Trial(trials.hit)) .|> Float32
                    end
                end

                Og_h = [[o[:, lookup(o, Trial) .== true] for o in O] for O in Og]
                Og_m = [[o[:, lookup(o, Trial) .== false] for o in O] for O in Og]
            end
            layergroups = [1:5, [1, 2], [3], [4, 5]] # Superficial, middle (L4), deep
            Olfp = map(layergroups) do ls
                map(out) do o
                    map(o) do p
                        k = deepcopy(p[:θ])
                        idxs = parselayernum.(metadata(k)[:layernames]) .∈ [ls]
                        k = k[:, parent(idxs), :]
                        N = ZScore(k, dims = 1)
                        normalize!(k, N)
                        ret = dropdims(mean(k; dims = Depth), dims = Depth)
                        trials = p[:trials][1:size(k, :changetime), :]
                        trialtimes = trials.change_time_with_display_delay
                        @assert all(isapprox.(ustripall(lookup(k, :changetime)), trialtimes,
                                              atol = 0.05)) # These won't match exactly, because the data change times have been adjusted for rectification. This is OK.
                        set(ret, Dim{:changetime} => Trial(trials.hit))
                    end
                end
            end

            begin # * LDA input data and downsampling
                H = [getindex.(Og, s) for s in eachindex(Og[1])] # Order parameter
                Hlfp = map(Olfp) do O
                    [getindex.(O, s) for s in eachindex(O[1])]
                end

                # * Make sure temporal dims same size
                ts = intersect([intersect(lookup.(h, 𝑡)...) for h in H]...)
                H = [[_h[𝑡 = At(ts)] for _h in h] for h in H]
                H = stack.([Structure(structures)], H, dims = 3)
                H = permutedims.(H, [(1, 3, 2)])

                Hlfp = map(Hlfp) do H
                    H = [[_h[𝑡 = At(ts)] for _h in h] for h in H]
                    H = stack.([Structure(structures)], H, dims = 3)
                    H = permutedims.(H, [(1, 3, 2)])
                end

                cg = h -> coarsegrain(h; dims = 1, newdim = 4)
                function consist(H)
                    HH = pmap(H) do h
                        dropdims(nansafe(mean; dims = 4)((cg ∘ cg ∘ cg)(h)); dims = 4)
                    end
                    nanprop = maximum(sum.([isnan.(x) for x in HH]) ./
                                      length.([isnan.(x) for x in HH]))
                    @info "$(round(nanprop*100, sigdigits=2))% of order parameters are NaN"
                    @assert nanprop < 0.05

                    map(HH) do h
                        h[isnan.(h)] .= 0.0 # Ok. Could be better.
                    end
                    @assert maximum(sum.([isnan.(x) for x in HH])) .== 0
                    return HH
                end
                H = consist(H)
                Hlfp = consist.(Hlfp)
                GC.gc()
            end

            if isfile(hyperfile) && isfile(datafile)
                D = jldopen(datafile, "r") do D
                    map(keys(D)) do d
                        d => D[d]
                    end |> Dict
                end
            else # Run calculations; needs to be on a cluster
                if !isfile(hyperfile)
                    begin # * Get a rough estimate for a good hyperparameter. Currently on pre-offset data. This gives us ~0.5 as a good estimate
                        hfile = calcdir("plots", "hyperparameters", "theta_waves_task.jld2")
                        hyperr = pmap(SpatiotemporalMotifs.tuneclassifier, H)
                        tagsave(hfile, @strdict hyperr)
                        if isfile(hyperfile)
                            hyperr = load(hyperfile, "hyperr")
                            f = Figure()
                            ax = Axis(f[1, 1]; xlabel = "Regularization coefficient",
                                      ylabel = "Balanced accuracy",
                                      title = "Hyperparameter tuning")
                            scatter!(ax, first.(hyperr), last.(hyperr))
                            f
                            save(calcdir("plots", "hyperparameters",
                                         "theta_waves_task.pdf"), f)
                        end
                    end
                end

                begin # * Single-subject classifications, returning 5-fold balanced accuracy. Takes ages, about 1 hour
                    bac_pre = pmap(H) do h
                        h = h[𝑡 = -0.25u"s" .. 0.25u"s"]
                        bac = classify_kfold(h; regcoef, k = folds, repeats)
                    end

                    bac_post = pmap(H) do h
                        h = h[𝑡 = 0.25u"s" .. 0.75u"s"]
                        bac = classify_kfold(h; regcoef, k = folds, repeats)
                    end

                    bac_sur = pmap(H) do h
                        hs = h[𝑡 = -0.25u"s" .. 0.25u"s"]
                        idxs = randperm(size(h, Trial))
                        h = set(h, Trial => lookup(h, Trial)[idxs])
                        bac = classify_kfold(h; regcoef, k = folds, repeats)
                    end

                    bac_lfp_pre = map(Hlfp) do H
                        pmap(H) do h # * Mean LFP, pre-offset
                            h = h[𝑡 = -0.25u"s" .. 0.25u"s"]
                            bac = classify_kfold(h; regcoef, k = folds, repeats)
                        end
                    end
                    bac_lfp_post = map(Hlfp) do H
                        pmap(H) do h # * Mean LFP, pre-offset
                            h = h[𝑡 = 0.25u"s" .. 0.75u"s"]
                            bac = classify_kfold(h; regcoef, k = folds, repeats)
                        end
                    end

                    D = @strdict bac_pre bac_post bac_sur bac_lfp_pre bac_lfp_post regcoef folds repeats
                end

                begin # * Map of region-wise weightings
                    W = pmap(H) do h
                        h = h[𝑡 = -0.25u"s" .. 0.75u"s"] # !!!
                        N, M = classifier(h; regcoef) # !!!
                        w = projection(M)
                        w = reshape(w, size(h)[1:2])
                        return ToolsArray(w, dims(h)[1:2])
                    end
                    W = stack(SessionID(oursessions), W, dims = 3)
                    W = W ./ maximum(abs.(W)) # ! Normalized

                    Wlfp = map(Hlfp) do H
                        WL = pmap(H) do h
                            h = h[𝑡 = -0.25u"s" .. 0.75u"s"] # !!!
                            N, M = classifier(h; regcoef) # !!!
                            Wl = projection(M)
                            Wl = reshape(Wl, size(h)[1:2])
                            return ToolsArray(Wl, dims(h)[1:2])
                        end
                        WL = stack(SessionID(oursessions), WL, dims = 3)
                        WL = WL ./ maximum(abs.(WL)) # ! Normalized
                    end
                end
                push!(D, "W" => W .|> Float32)
                push!(D, "Wlfp" => Wlfp .|> Float32)
                tagsave(datafile, D)
            end

            Natural_Images = @strdict Og layergroups D
        end
        order_parameters = @strdict flashes Natural_Images passive
        out = [] # Clear for memory
        GC.gc()
    end

    begin # * Now do the unified wavenumbers
        stimulus = r"Natural_Images"
        uni = load_uni(; stimulus, vars)

        Natural_Images = map(uni) do u
            k = u[:k]
            ω = u[:ω]
            k = uconvert.(u"mm^-1", k)
            k[ustripall(ω) .< 0] .= NaN * unit(eltype(k)) # Mask out negative frequency periods
            @info "$(round(100*sum(ω .< 0u"Hz")/length(ω), sigdigits=2))% of data has a negative frequency"

            structure = metadata(k)[:structure]
            ints = u[:layerints]

            m_hit = nansafe(median, dims = Trial)(k[:, :, lookup(k, Trial) .== true])
            m_hit = dropdims(m_hit, dims = Trial)[𝑡(SpatiotemporalMotifs.INTERVAL)]

            m_miss = nansafe(median, dims = Trial)(k[:, :, lookup(k, Trial) .== false])
            m_miss = dropdims(m_miss, dims = Trial)[𝑡(SpatiotemporalMotifs.INTERVAL)]

            # * Calculate some stats of this figure
            if structure == mainstructure
                open(statsfile, "a+") do file
                    write(file, "## $mainstructure average wavenumbers\n")
                    write(file,
                          "Average wavenumber (median ± IQR) = $(only(nansafe(median)(k[:]))) ± $(nansafe(iqr)(k[:])|>only)\n")
                    write(file,
                          "Average wavenumber magnitude (median ± IQR) = $(only(nansafe(median)(abs.(k[:])))) ± $(nansafe(iqr)(abs.(k[:]))|>only)\n")
                end
            end

            # * Extract csd data
            csd = u[:csd]
            csd = uconvert.([u"cm^(-2)"], csd) # csd has dropped the volts dims but kept um dims
            csd_hit = nansafe(median, dims = Trial)(csd[:, :, lookup(csd, Trial) .== true])
            csd_hit = dropdims(csd_hit, dims = Trial)[𝑡(SpatiotemporalMotifs.INTERVAL)]
            csd_miss = nansafe(median, dims = Trial)(csd[:, :,
                                                         lookup(csd, Trial) .== false])
            csd_miss = dropdims(csd_miss, dims = Trial)[𝑡(SpatiotemporalMotifs.INTERVAL)]

            return (@strdict m_hit m_miss structure ints csd_hit csd_miss)
        end

        # ? Flashes
        stimulus = "flash_250ms"
        flash_uni = load_uni(; stimulus, vars)
        flashes = map(flash_uni) do u
            k = u[:k]
            ω = u[:ω]
            k = uconvert.(u"mm^-1", k)
            k[ustripall(ω) .< 0] .= NaN * unit(eltype(k)) # Mask out negative frequency periods
            structure = metadata(k)[:structure]

            m_flashes = nansafe(median, dims = (Trial, SessionID))(k)
            m_flashes = dropdims(m_flashes, dims = (Trial, SessionID))[𝑡(SpatiotemporalMotifs.INTERVAL)]
            ints = u[:layerints]

            csd = u[:csd]
            csd = uconvert.([u"cm^(-2)"], csd) # k has dropped the volts dims but kept um dims
            csd_flashes = nansafe(median, dims = (Trial, SessionID))(csd)
            csd_flashes = dropdims(csd_flashes, dims = (Trial, SessionID))[𝑡(SpatiotemporalMotifs.INTERVAL)]

            return (@strdict m_flashes csd_flashes structure ints)
        end
        wavenumbers = @strdict Natural_Images flashes
    end
    return (@strdict wavenumbers order_parameters)
end

# begin # * Current source density
#     csd = getindex.(uni, :csd)[2][:, :, :]
#     csd = dropdims(mean(csd, dims = Trial); dims = Trial)
#     csd = upsample(ustripall(csd), 5, Depth) # For VISl
#     # csd = set(csd, csd[:, end:-1:1])
#     f = Figure()
#     ax = Axis(f[1, 1]; yreversed = true)
#     plotlayermap!(ax, csd)
#     plotlayerints!(ax, uni[2][:layerints])
#     # heatmap(decompose(csd)...; colormap = binarysunset, axis = (; yreversed = true))
#     f
# end

begin # * Set up main figure
    fullfig = SixPanel()
    mf = fullfig[1, :]
    mfs = subdivide(mf, 1, 4)
    fig = fullfig[2:3, :]
    gs = subdivide(fig, 2, 2)
end

begin # * Wavenumbers
    @info "Plotting wavenumbers"
    begin # * Supplemental material: average wavenumbers in each region
        f = Figure(size = (1440, 720) .* 1.25)
        fgs = subdivide(f, 3, 2)

        map(plot_data["wavenumbers"]["Natural_Images"], fgs) do data, fg
            @unpack m_hit, m_miss, structure, ints = data

            # * Hit
            ax = Axis(fg[1, 1], yreversed = true)
            ax.limits = (nothing, ylims)
            ax.title = structure * ": hit"
            p = plotlayermap!(ax, m_hit, ints; arrows = true,
                              colorrange = maincolorrange) |>
                first
            if fg === fgs[5] || fg === fgs[6]
                ax.xlabel = "Time (s)"
            end

            if structure == mainstructure # * Plot into main figure
                ax = Axis(mfs[1], yreversed = true)
                ax.xlabel = "Time (s)"
                ax.limits = (nothing, ylims)
                ax.title = structure * ": hit"
                p = plotlayermap!(ax, m_hit, ints; arrows = true,
                                  colorrange = maincolorrange) |>
                    first
                c = Colorbar(mfs[end]; colorrange = maincolorrange,
                             colormap = defaultcolormap,
                             highclip = defaultcolormap[end], lowclip = defaultcolormap[1])
                c.label = "θ wavenumber ($(unit(eltype(m_hit))))"
            end

            # * Miss
            ax = Axis(fg[1, 2], yreversed = true)
            ax.limits = (nothing, ylims)
            ax.title = structure * ": miss"
            # ax.yticklabelsvisible = false
            p = plotlayermap!(ax, m_miss, ints; arrows = true,
                              colorrange = maincolorrange) |>
                first
            if fg === fgs[5] || fg === fgs[6]
                ax.xlabel = "Time (s)"
            end

            if structure == mainstructure # * Plot into main figure
                ax = Axis(mfs[2], yreversed = true)
                ax.xlabel = "Time (s)"
                ax.limits = (nothing, ylims)
                ax.title = structure * ": miss"
                p = plotlayermap!(ax, m_miss, ints; arrows = true,
                                  colorrange = maincolorrange) |>
                    first
            end
        end
        addlabels!(fgs[:], f, labelformat)
        display(f)
    end

    begin # * Wavenumber for flashes
        map(plot_data["wavenumbers"]["flashes"], fgs) do data, fg
            @unpack m_flashes, structure, ints = data

            ax = Axis(fg[1, 3], yreversed = true)
            ax.limits = (nothing, ylims)
            ax.title = structure * ": flashes"
            p = plotlayermap!(ax, m_flashes, ints; arrows = true,
                              colorrange = maincolorrange) |>
                first
            c = Colorbar(fg[1, 4]; colorrange = maincolorrange,
                         colormap = defaultcolormap,
                         highclip = defaultcolormap[end], lowclip = defaultcolormap[1])
            c.label = "θ wavenumber ($(unit(eltype(m_flashes))))"
            if structure == mainstructure
                colorrange = maincolorrange
                ax = Axis(mfs[3], yreversed = true, xlabel = "Time (s)")
                ax.limits = (nothing, ylims)
                ax.title = structure * ": flashes"
                p = plotlayermap!(ax, m_flashes, ints; arrows = true,
                                  colorrange = maincolorrange) |>
                    first
            end
            if fg === fgs[5] || fg === fgs[6]
                ax.xlabel = "Time (s)"
            end
        end
        colgap!(f.layout, 1, Relative(0.05))
        addlabels!(fgs[:], f, labelformat)
        display(f)
        # wsave(plotdir("fig3", "supplemental_wavenumber_flash.pdf"), f)
        wsave(plotdir("fig3", "supplemental_wavenumber.pdf"), f)
    end
end

begin # * Supplemental material: csd in each region
    @info "Plotting CSDs"
    f = Figure(size = (1440, 720) .* 1.25)
    fgs = subdivide(f, 3, 2)

    map(plot_data["wavenumbers"]["Natural_Images"], fgs) do data, fg
        @unpack csd_hit, csd_miss, ints, structure = data

        # * Hit
        ax = Axis(fg[1, 1], yreversed = true)
        ax.limits = (nothing, ylims)
        ax.title = structure * ": hit"

        # colorrange = maximum(abs.(ustripall(m))) * [-1, 1]
        p = plotlayermap!(ax, csd_hit, ints; arrows = false, colorrange = csdcolorrange,
                          colormap = defaultcolormap) |> first
        if fg === fgs[5] || fg === fgs[6]
            ax.xlabel = "Time (s)"
        end

        # if structure == mainstructure # * Plot into main figure
        #     ax = Axis(mfs[1], yreversed = true)
        #     ax.xlabel = "Time (s)"
        #     ax.limits = (nothing, ylims)
        #     ax.title = structure * ": hit"
        #     m = nansafe(median, dims = Trial)(k[:, :, lookup(k, Trial) .== true])
        #     m = dropdims(m, dims = Trial)[𝑡(SpatiotemporalMotifs.INTERVAL)]
        #     ints = uni[i][:layerints]
        #     p = plotlayermap!(ax, m, ints; arrows = true, colorrange) |> first
        #     c = Colorbar(mfs[end]; colorrange = maincolorrange, colormap = defaultcolormap,
        #                  highclip = defaultcolormap[end], lowclip = defaultcolormap[1])
        #     c.label = "CSD ($(unit(eltype(k))))"
        # end

        # * Miss
        ax = Axis(fg[1, 2], yreversed = true)
        ax.limits = (nothing, ylims)
        ax.title = structure * ": miss"
        p = plotlayermap!(ax, csd_miss, ints; arrows = false, colorrange = csdcolorrange,
                          colormap = defaultcolormap) |> first
        if fg === fgs[5] || fg === fgs[6]
            ax.xlabel = "Time (s)"
        end

        # if structure == mainstructure # * Plot into main figure
        #     ax = Axis(mfs[2], yreversed = true)
        #     ax.xlabel = "Time (s)"
        #     ax.limits = (nothing, ylims)
        #     ax.title = structure * ": miss"
        #     m = nansafe(median; dims = Trial)(k[:, :, lookup(k, Trial) .== false])
        #     m = dropdims(m, dims = Trial)[𝑡(SpatiotemporalMotifs.INTERVAL)]
        #     ints = uni[i][:layerints]
        #     p = plotlayermap!(ax, m, ints; arrows = true, colorrange = maincolorrange) |>
        #         first
        # end
    end
    begin # * CSD for flashes
        map(plot_data["wavenumbers"]["flashes"], fgs) do data, fg
            @unpack csd_flashes, structure, ints = data

            ax = Axis(fg[1, 3], yreversed = true)
            ax.limits = (nothing, ylims)
            ax.title = structure * ": flashes"
            # colorrange = maximum(abs.(ustripall(m))) * [-1, 1]
            p = plotlayermap!(ax, csd_flashes, ints; arrows = false,
                              colorrange = csdcolorrange,
                              colormap = defaultcolormap) |> first
            c = Colorbar(fg[1, 4]; colorrange = csdcolorrange,
                         colormap = defaultcolormap,
                         highclip = defaultcolormap[end], lowclip = defaultcolormap[1])
            c.label = "CSD ($(u"V"*unit(eltype(csd_flashes))))" # Correct for dropped units by multiplying with V

            if fg === fgs[5] || fg === fgs[6]
                ax.xlabel = "Time (s)"
            end
        end
    end
    addlabels!(fgs[:], f, labelformat)
    colgap!(f.layout, 1, Relative(0.05))
    display(f)
    wsave(plotdir("fig3", "supplemental_csd.pdf"), f)
end

begin # * Order parameters
    @info "Plotting order parameters"
    function orderparameter(Og)
        out = map(Og) do O
            o = dropdims.(nansafe(mean; dims = Trial).(O), dims = Trial)
            sessionids = [metadata(_o)[:sessionid] for _o in O]
            o = [_o[(end - minimum(length.(o)) + 1):end] for _o in o]
            stack(SessionID(sessionids), o)
        end
        ts = lookup.(out, 𝑡)
        minterval = maximum(minimum.(ts)) .. minimum(maximum.(ts))
        out = [o[𝑡 = minterval] for o in out]
    end
    begin # * Flashes
        @unpack Og = plot_data["order_parameters"]["flashes"]
        begin # * Plot the mean order parameter across time
            Ō = orderparameter(Og)

            ax = Axis(gs[2]; xlabel = "Time (s)",
                      ylabel = rich("Mean order parameter ",
                                    rich("R", subscript("θ"), font = "Times Italic")),
                      title = "Order parameter during flashes",
                      xautolimitmargin = (0, 0), xminorticksvisible = true,
                      xminorticks = IntervalsBetween(5), yminorticksvisible = true,
                      yminorticks = IntervalsBetween(5))
            hlines!(ax, [0]; color = (:black, 0.5), linestyle = :dash, linewidth = 2)
            vlines!(ax, [0, 0.25]; color = (:black, 0.5), linestyle = :dash, linewidth = 2)
            for (i, O) in reverse(collect(enumerate(Ō)))
                structure = metadata(O)[:structure]
                O = O[𝑡(SpatiotemporalMotifs.INTERVAL)]
                μ = dropdims(mean(O, dims = SessionID), dims = SessionID)
                σ = dropdims(std(O, dims = SessionID), dims = SessionID)
                σ = σ ./ 2
                bargs = [times(μ), μ .- σ, μ .+ σ] .|> ustripall .|> collect
                band!(ax, bargs..., color = (structurecolors[i], 0.3), label = structure)
                lines!(ax, times(μ) |> ustripall, μ |> ustripall |> collect,
                       color = (structurecolors[i], 0.7),
                       label = structure)
            end
            l = axislegend(ax, position = :lt, nbanks = 2, framevisible = true,
                           labelsize = 12,
                           merge = true)
            reverselegend!(l)
        end
    end
    begin # * Task stimulus (natural images)
        @unpack Og, layergroups, D = plot_data["order_parameters"]["Natural_Images"]
        @unpack bac_pre, bac_post, bac_sur, bac_lfp_pre, bac_lfp_post, regcoef, folds, repeats, W, Wlfp = D

        Og_h = [[o[:, lookup(o, Trial) .== true] for o in O] for O in Og]
        Og_m = [[o[:, lookup(o, Trial) .== false] for o in O] for O in Og]

        begin # * Plot the mean order parameter across time
            Ō = orderparameter(Og)
            Ō_h = orderparameter(Og_h)
            Ō_m = orderparameter(Og_m)

            ax = Axis(gs[1]; xlabel = "Time (s)",
                      ylabel = rich("Mean order parameter ",
                                    rich("R", subscript("θ"), font = "Times Italic")),
                      title = "Order parameter during task",
                      xautolimitmargin = (0, 0), xminorticksvisible = true,
                      xminorticks = IntervalsBetween(5), yminorticksvisible = true,
                      yminorticks = IntervalsBetween(5))
            hlines!(ax, [0]; color = (:black, 0.5), linestyle = :dash, linewidth = 2)
            vlines!(ax, [0, 0.25]; color = (:black, 0.5), linestyle = :dash, linewidth = 2)
            for (i, O) in reverse(collect(enumerate(Ō)))
                structure = metadata(O)[:structure]
                O = O[𝑡(SpatiotemporalMotifs.INTERVAL)]
                μ = dropdims(mean(O, dims = SessionID), dims = SessionID)
                σ = dropdims(std(O, dims = SessionID), dims = SessionID)
                σ = σ ./ 2
                bargs = [times(μ), μ .- σ, μ .+ σ] .|> ustripall .|> collect
                band!(ax, bargs..., color = (structurecolors[i], 0.3), label = structure)
                lines!(ax, times(μ) |> ustripall, μ |> ustripall |> collect,
                       color = (structurecolors[i], 0.7),
                       label = structure)
            end
            l = axislegend(ax, position = :lt, nbanks = 2, framevisible = true,
                           labelsize = 12,
                           merge = true)
            reverselegend!(l)
        end
        begin # * Order parameter correlation to hierarchy
            N = 1e6
            statsfile = plotdir("fig3", "theta_order_parameter.txt")
            close(open(statsfile, "w")) # Create the file or clear it

            x = getindex.([SpatiotemporalMotifs.hierarchy_scores], structures)
            _y = stack(Structure(structures), Ō)
            _y = stack(Depth([0]), [_y])

            for T in [0.0u"s" .. 0.15u"s", 0.25u"s" .. 0.4u"s"]
                y = _y[𝑡 = T]
                y = dropdims(mean(y; dims = 𝑡); dims = 𝑡)
                open(statsfile, "a+") do file
                    write(file, "\n# $T\n")
                end

                # * Group level
                μ, σ, 𝑝 = hierarchicalkendall(x, y, :group; N) .|> first
                open(statsfile, "a+") do file
                    write(file, "\n## Group level")
                    write(file, "\nmedian τ = $μ")
                    write(file, "\n95\\% conf. = $σ")
                    write(file, "\n𝑝 = $𝑝")
                    write(file, "\nN = $N")
                    write(file, "\n")
                end

                # * Individual level
                μ, σ, 𝑝 = hierarchicalkendall(x, y, :individual; N) .|> first
                open(statsfile, "a+") do file
                    write(file, "\n## Individual level")
                    write(file, "\nmedian τ = $μ")
                    write(file, "\nIQR = $(σ[2] - σ[1])")
                    write(file, "\n𝑝 = $𝑝")
                    write(file, "\n")
                end
            end
        end

        begin # * Plot classification performance
            ax = Axis(fig[2, 1];
                      xticks = ([2, 5, 7], ["Post-offset", "Pre-offset", "Null"]),
                      ylabel = "Balanced accuracy", title = "Hit/miss classification",
                      limits = ((0.5, 7.5), (0.35, 0.95)), xminorticksvisible = false)
            boxargs = (; width = 0.75, strokewidth = 5, whiskerwidth = 0.2,
                       strokecolor = (:gray, 0.0), whiskerlinewidth = 6) # !!!! Show outliers??

            vspan!(ax, 0.5, 3.5; color = (california, 0.2))
            vspan!(ax, 3.5, 6.5; color = (cucumber, 0.2))
            vspan!(ax, 6.5, 7.5; color = (:gray, 0.2))
            hlines!(ax, [0.5]; color = (:black, 0.5), linestyle = :dash, linewidth = 2)

            begin # * Write statistics
                open(statsfile, "a+") do file
                    write(file, "\n# Classification performance\n")
                    write(file, "\n## Order parameter pre-offset\n")
                    write(file, "\nMedian BAc = $(median(bac_pre))")
                    write(file, "\nIQR = $(iqr(bac_pre))")
                    𝑝 = HypothesisTests.pvalue(HypothesisTests.MannWhitneyUTest(bac_pre,
                                                                                bac_sur);
                                               tail = :right)
                    write(file, "\nU-test, right-sided 𝑝 to sur= $𝑝")
                    write(file, "\n")

                    write(file, "\n## Order parameter post-offset\n")
                    write(file, "\nMedian BAc = $(median(bac_post))")
                    write(file, "\nIQR = $(iqr(bac_post))")
                    𝑝 = HypothesisTests.pvalue(HypothesisTests.MannWhitneyUTest(bac_pre,
                                                                                bac_post);
                                               tail = :both)
                    write(file, "\nU-test, two-sided 𝑝 to pre. = $𝑝")
                    write(file, "\n")

                    write(file, "\n## Mean LFP pre-offset \n")
                    write(file, "\nMedian BAc = $(median(bac_lfp_pre[1]))")
                    write(file, "\nIQR = $(iqr(bac_lfp_pre[1]))")
                    𝑝 = HypothesisTests.pvalue(HypothesisTests.MannWhitneyUTest(bac_pre,
                                                                                bac_lfp_pre[1]);
                                               tail = :right)
                    write(file, "\nU-test, two-sided 𝑝 to pre. = $𝑝")
                    write(file, "\n")

                    write(file, "\n## Order parameter pre-offset null\n")
                    write(file, "\nMedian BAc = $(median(bac_sur))")
                    write(file, "\nIQR = $(iqr(bac_sur))")
                    write(file, "\n")
                end
            end

            boxplot!(ax, fill(1, length(bac_post)), bac_post; boxargs...,
                     color = cornflowerblue,
                     label = "Order parameter")
            boxplot!(ax, fill(2, length(bac_lfp_post[1])), bac_lfp_post[1]; boxargs...,
                     color = juliapurple, label = "Mean LFP")
            boxplot!(ax,
                     vcat([fill(3 + i, length(bac_lfp_post[1])) for i in [-0.3, 0, 0.3]]...),
                     reverse(vcat(bac_lfp_post[2:end]...)); boxargs..., width = 0.3,
                     color = crimson, label = "Compartmental LFP")
            text!(ax, 3 .+ [-0.3, 0, 0.3], [0.9, 0.9, 0.9]; text = reverse(["S", "M", "D"]),
                  align = (:center, :center))

            boxplot!(ax, fill(4, length(bac_pre)), bac_pre; boxargs...,
                     color = cornflowerblue,
                     label = "Order parameter")
            boxplot!(ax, fill(5, length(bac_lfp_pre[1])), bac_lfp_pre[1]; boxargs...,
                     color = juliapurple, label = "Mean LFP")
            boxplot!(ax,
                     vcat([fill(6 + i, length(bac_lfp_pre[1])) for i in [-0.3, 0, 0.3]]...),
                     reverse(vcat(bac_lfp_pre[2:end]...)); boxargs..., width = 0.3,
                     color = crimson,
                     label = "Compartmental LFP")
            text!(ax, 6 .+ [-0.3, 0, 0.3], [0.8, 0.8, 0.8]; text = reverse(["S", "M", "D"]),
                  align = (:center, :center))

            boxplot!(ax, fill(7, length(bac_sur)), bac_sur; color = :gray, boxargs...)
            axislegend(ax; merge = true, labelsize = 10)
        end
        begin # * Plot region-wise weightings
            ax = Axis(fig[2, 2]; xlabel = "Time (s)", ylabel = "Normalized LDA weight",
                      title = "Regional classification weights")
            vlines!(ax, [0, 0.25], color = (:black, 0.2), linestyle = :dash)
            hlines!(ax, [0], color = (:black, 0.5), linewidth = 2)
            for structure in (collect(lookup(W, Structure)))
                ws = W[Structure = At(structure)]
                ws = upsample(ws, 5, 1)
                ts = ustripall(lookup(ws, 𝑡))
                μ = mean(ws, dims = SessionID) |> vec
                σ = std(ws, dims = SessionID) ./ sqrt(size(ws, SessionID)) |> vec |>
                    ustripall
                band!(ax, ts, μ - σ, μ + σ; color = (structurecolormap[structure], 0.3),
                      label = structure)
                lines!(ax, ts, μ; color = structurecolormap[structure], label = structure)
            end
            tightlimits!(ax)
            l = axislegend(ax, position = :lt, nbanks = 2, framevisible = true,
                           labelsize = 12,
                           merge = true)
            # reverselegend!(l)
        end
    end

    begin # * Save
        addlabels!(fullfig, labelformat)
        rowsize!(fullfig.layout, 1, Relative(0.25))
        wsave(plotdir("fig3", "theta_order_parameter.pdf"), fullfig)
    end
end

begin # * Supplementary figure: contrast between active and passive stimulus
    begin # * Active rder parameter
        @unpack Og = plot_data["order_parameters"]["Natural_Images"]
        Ō_active = orderparameter(Og)

        @unpack Og = plot_data["order_parameters"]["passive"]
        Ō_passive = orderparameter(Og)

        begin # * Plot the mean order parameter across time
            ax = Axis(gs[2]; xlabel = "Time (s)",
                      ylabel = rich("Mean order parameter ",
                                    rich("R", subscript("θ"), font = "Times Italic")),
                      title = "Order parameter during flashes",
                      xautolimitmargin = (0, 0), xminorticksvisible = true,
                      xminorticks = IntervalsBetween(5), yminorticksvisible = true,
                      yminorticks = IntervalsBetween(5))
            hlines!(ax, [0]; color = (:black, 0.5), linestyle = :dash, linewidth = 2)
            vlines!(ax, [0, 0.25]; color = (:black, 0.5), linestyle = :dash, linewidth = 2)
            for (i, O) in reverse(collect(enumerate(Ō)))
                structure = metadata(O)[:structure]
                O = O[𝑡(SpatiotemporalMotifs.INTERVAL)]
                μ = dropdims(mean(O, dims = SessionID), dims = SessionID)
                σ = dropdims(std(O, dims = SessionID), dims = SessionID)
                σ = σ ./ 2
                bargs = [times(μ), μ .- σ, μ .+ σ] .|> ustripall .|> collect
                band!(ax, bargs..., color = (structurecolors[i], 0.3), label = structure)
                lines!(ax, times(μ) |> ustripall, μ |> ustripall |> collect,
                       color = (structurecolors[i], 0.7),
                       label = structure)
            end
            l = axislegend(ax, position = :lt, nbanks = 2, framevisible = true,
                           labelsize = 12,
                           merge = true)
            reverselegend!(l)
        end
    end
end
