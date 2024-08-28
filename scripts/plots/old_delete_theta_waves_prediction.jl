#! /bin/bash
#=
exec julia -t auto "${BASH_SOURCE[0]}" "$@"
=#
using DrWatson
@quickactivate "SpatiotemporalMotifs"

import DimensionalData: metadata
using Optim
using SpatiotemporalMotifs
@preamble
import TimeseriesTools: freqs
set_theme!(foresight(:physics))
Random.seed!(42)

stimulus = r"Natural_Images"
vars = [:k]

session_table = load(datadir("session_table.jld2"), "session_table")
oursessions = session_table.ecephys_session_id

path = datadir("calculations")
Q = calcquality(path)[structure = At(structures)]
quality = mean(Q[stimulus = At(stimulus)])

config = @strdict stimulus vars
data, file = produce_or_load(produce_out(Q), config, datadir(); filename = savepath,
                             prefix = "out")
out = data["out"]

# * Filter out poor-performing mice
performance = load_performance(; path)
newsessions = performance.sessionid[performance.mean_dprime .> 1]
out = [filter(x -> x.sessionid ∈ newsessions, o) for o in out]
GC.gc()

begin # * Calculate a global order parameter at each time point
    Og = map(out) do o
        O = map(o) do p
            k = p[:k]
            ret = dropdims(mean(sign.(k), dims = Depth), dims = Depth)
            trials = p[:trials]
            trialtimes = trials.change_time_with_display_delay
            @assert all(isapprox.(ustripall(lookup(k, :changetime)), trialtimes,
                                  atol = 0.05)) # These won't match exactly, because the data change times have been adjusted for rectification. This is OK.
            ret = set(ret, Dim{:changetime} => Trial)
            set(ret, DimensionalData.format(Trial(trials.hit), lookup(ret, Trial)))
        end
    end

    Og_h = [[o[:, lookup(o, Trial) .== true] for o in O] for O in Og]
    Og_m = [[o[:, lookup(o, Trial) .== false] for o in O] for O in Og]
    sessionids = [o[:sessionid] for o in out[1]]
end

begin # * LDA input data and downsampling
    H = [getindex.(Og, s) for s in eachindex(Og[1])] # Grouped by subject

    # * Make sure temporal dims same size
    ts = intersect([intersect(lookup.(h, 𝑡)...) for h in H]...)
    H = [[_h[Ti = At(ts)] for _h in h] for h in H]
    H = stack.([Dim{:structure}(structures)], H, dims = 3)
    H = permutedims.(H, [(1, 3, 2)])

    H = [h[1:10:end, :, :] for h in H]
end

begin # * What sort of hyperparameter? Currently on including-offset data
    function tuneclassifier(H; r0 = 0.5, kwargs...)
        r0 = log10(r0)
        objective = r -> -classify_kfold(H; regcoef = exp10(only(r)), kwargs...)
        o = optimize(objective, [r0], ParticleSwarm(), Optim.Options(; iterations = 10))
        return exp10(only(o.minimizer)), -o.minimum
    end
    hfile = datadir("hyperparameters", "theta_waves_prediction.jld2")
    if haskey(ENV, "JULIA_DISTRIBUTED")
        procs = addprocs(10; ncpus = 8, mem = 45, walltime = 96)
        hyperr = pmap(tuneclassifier, H)
        tagsave(hfile, @strdict hyperr)
    end
    if isfile(hfile)
        hyperr = load(hfile, "hyperr")
        f = Figure()
        ax = Axis(f[1, 1]; xlabel = "Regularization coefficient",
                  ylabel = "Balanced accuracy",
                  title = "Hyperparameter tuning")
        scatter!(ax, first.(hyperr), last.(hyperr))
        f
        save(datadir("hyperparameters", "theta_waves_prediction.pdf"), f)
    end
end

begin # * Single-subject classifications
    bac_pre = progressmap(H; description = "Pre-stim") do h
        h = h[Ti = -0.25u"s" .. 0.25u"s"]
        bac = classify_kfold(h; regcoef = 0.5)
    end

    bac_post = progressmap(H; description = "Post-stim") do h
        h = h[Ti = 0.25u"s" .. 0.75u"s"]
        bac = classify_kfold(h; regcoef = 0.5)
    end

    bac_sur = progressmap(H; description = "Surrogate") do h
        h = h[Ti = -0.25u"s" .. 0.25u"s"]
        idxs = randperm(size(h, Trial))
        h = set(h, Trial => lookup(h, Trial)[idxs])
        bac = classify_kfold(h; regcoef = 0.5)
    end
end
begin
    f = Figure()
    ax = Axis(f[1, 1]; xlabel = "Balanced accuracy", ylabel = "Density",
              title = "Balanced accuracy distribution")
    ziggurat!(ax, bac_sur; bins = 10, label = "Null", normalizaton = :pdf, color = crimson)
    ziggurat!(ax, bac_post; bins = 10, label = "Post-stimulus", normalizaton = :pdf,
              color = cucumber)
    ziggurat!(ax, bac_pre; bins = 10, label = "Pre-stimulus", normalizaton = :pdf,
              color = cornflowerblue)
    axislegend(ax, position = :lt)
    tightlimits!(ax)
    f
end

begin # * Map of layer-wise weightings
    W = map(H) do h
        h = h[Ti = -0.25u"s" .. 0.75u"s"] # !!!
        N, M = classifier(h; regcoef = 0.5) # !!!
        W = projection(M)
        W = reshape(W, size(h)[1:2])
        return ToolsArray(W, dims(h)[1:2])
    end
    W = stack(SessionID(sessionids), W, dims = 3)
    W = W ./ maximum(abs.(W))
end
begin
    f = Figure()
    ax = Axis(f[1, 1]; xlabel = "Time (s)", ylabel = "Normalized LDA weight",
              title = "Regional classification weights")
    vlines!(ax, [0, 0.25], color = (:black, 0.2), linestyle = :dash)
    hlines!(ax, [0], color = (:black, 0.5), linewidth = 2)
    for structure in lookup(W, :structure)
        ws = W[structure = At(structure)]
        ws = upsample(ws, 5, 1)
        ts = ustripall(lookup(ws, 𝑡))
        μ = mean(ws, dims = SessionID) |> vec
        σ = std(ws, dims = SessionID) ./ sqrt(size(ws, SessionID)) |> ustripall |> vec
        band!(ax, ts, μ - σ, μ + σ; color = (structurecolormap[structure], 0.3))
        lines!(ax, ts, μ; color = structurecolormap[structure], label = structure)
    end
    tightlimits!(ax)
    axislegend(ax, position = :lt)
    f
end

begin
    begin # * Compare cross-validated balanced accuracy to dprime
        dprimes = [o[:performance_metrics]["mean_dprime_engaged"] for o in out[1]]
        r = round(corspearman(dprimes, bac); sigdigits = 2)
        f = Figure()
        ax = Axis(f[1, 1]; xlabel = "Task performance (d-prime engaged)",
                  ylabel = "Classification accuracy", title = "ρ =  $r")
        scatter!(ax, dprimes, bac)
        f
    end
    begin # * The most correlated performance metric?
        metrics = out[1][1][:performance_metrics] |> keys |> collect
        ρ = map(metrics) do m
            x = [o[:performance_metrics][m] for o in out[1]]
            corspearman(x, bac)
        end
        abs_ρ = abs.(ρ)
        D = DataFrame()
        @pack! D = metrics, ρ, abs_ρ
        sort!(D, [:abs_ρ], rev = true)
        D
    end
    begin # * And to reaction time
        reactiontime = [o[:trials].lick_latency for o in out[1]]
        reactiontime = map(x -> iqr(filter(!isnan, x)), reactiontime)
        r = round(corspearman(reactiontime, bac); sigdigits = 2)
        f = Figure()
        ax = Axis(f[1, 1]; xlabel = "Iqr of reaction time (s. from onset)",
                  ylabel = "Classification accuracy", title = "ρ =  $r")
        scatter!(ax, reactiontime, bac)
        f
    end
end
begin # * Now try regressing against reaction time... nothing
    # * We have reaction time for each trial. We:
    # 1. Take the LDA model trained on all hit/miss trials
    # 2. Calculate the 1D LDA projection for each hit trial
    # 3. Correlate the LDA projection with the response latency (for each hit trial)

    scores = map(H) do h
        h = h[Ti = -0.25u"s" .. 0.25u"s"]
        N, M = classifier(h; regcoef = 0.5f0)
        _h = reshape(h, prod(size(h)[1:2]), size(h)[3])
        _h = _h[:, lookup(h, Trial) .== true]
        predict(M, normalize(_h, N))[:]
    end

    reactions = map(out[1]) do o
        trials = o[:trials]
        trials.lick_latency[trials.hit]
    end

    ρs = map(scores, reactions) do s, r
        corspearman(s, r)
    end

    scores_sur = map(H) do h
        h = h[Ti = -0.25u"s" .. 0.25u"s"]
        idxs = randperm(size(h, Trial))
        h = set(h, Trial => lookup(h, Trial)[idxs])
        N, M = classifier(h; regcoef = 0.5f0)
        _h = reshape(h, prod(size(h)[1:2]), size(h)[3])
        _h = _h[:, lookup(h, Trial) .== true]
        predict(M, normalize(_h, N))[:]
    end
    ρs_sur = map(scores_sur, reactions) do s, r
        corspearman(s, r)
    end

    hist(ρs)
    hist!(ρs_sur)
    current_figure()
end # ...nothing. Now see `theta_waves_reactiontime.jl`
