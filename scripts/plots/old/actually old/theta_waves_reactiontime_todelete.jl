#! /bin/bash
#=
exec julia -t auto "${BASH_SOURCE[0]}" "$@"
=#
using DrWatson
@quickactivate "SpatiotemporalMotifs"

import DimensionalData: metadata
using SpatiotemporalMotifs
import MLJ
import MLJScikitLearnInterface
import MLJLinearModels
using RobustModels
@preamble
import TimeseriesTools: freqs
set_theme!(foresight(:physics))
Random.seed!(42)

stimulus = r"Natural_Images"
vars = [:k]

session_table = load(datadir("session_table.jld2"), "session_table")
oursessions = session_table.ecephys_session_id

path = datadir("calculations")
Q = calcquality(path)[Structure = At(structures)]
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

begin # * Regression input data and downsampling
    H = [getindex.(Og, s) for s in eachindex(Og[1])] # Grouped by subject
    reaction_times = [o[:trials].lick_latency for o in first(out)]

    # * Make sure temporal dims same size
    ts = intersect([intersect(lookup.(h, 𝑡)...) for h in H]...)
    H = [[_h[Ti = At(ts)] for _h in h] for h in H]
    H = stack.([Structure(structures)], H, dims = 3)
    H = permutedims.(H, [(1, 3, 2)])

    H = [h[1:15:end, :, :] for h in H] # Downsample so we have more observations than features
end

begin # * Regression as classification?
    bac = pmap(eachindex(H)) do i
        labels = deepcopy(reaction_times[i])
        labels[isnan.(labels)] .= NaN
        idxs = .!isnan.(labels) .& (labels .< 1) .* (labels .> 0.2)
        labels = labels .> nansafe(median)(labels)
        h = H[i][Ti = -0.25u"s" .. 0.25u"s"]
        labels = labels[idxs]
        h = h[:, :, idxs]
        h = set(h, Trial => labels)
        classify_kfold(h; regcoef = 1, k = 5)
    end
    bac_sur = pmap(eachindex(H)) do i
        labels = deepcopy(reaction_times[i])
        labels[isnan.(labels)] .= NaN
        idxs = .!isnan.(labels) .& (labels .< 1.5) .* (labels .> 0.1)
        labels = labels .> nansafe(median)(labels)
        h = H[i][Ti = -0.25u"s" .. 0.25u"s"]
        labels = labels[idxs]
        labels = labels[randperm(length(labels))]
        h = h[:, :, idxs]
        h = set(h, Trial => labels)
        classify_kfold(h; regcoef = 0.5, k = 5)
    end

    ziggurat(bac)
    ziggurat!(bac_sur)
    current_figure()
end
begin # * K-fold validation with Huber (outlier-robust) regression
    regcoef = 1e5
    bac = pmap(eachindex(H)) do i
        regress_kfold(H[i][Ti = -0.0u"s" .. 0.25u"s"], reaction_times[i]; regcoef,
                      k = 5)
    end
    bac_sur = pmap(eachindex(H)) do i
        regress_kfold(H[i][Ti = -0.0u"s" .. 0.25u"s"],
                      reaction_times[i][randperm(length(reaction_times[i]))]; regcoef,
                      k = 5)
    end
    ziggurat(bac)
    ziggurat!(bac_sur)
    current_figure()
end

# begin # * Test MLJ models
#     i = 11
#     # rhos = pmap(eachindex(H)) do i
#     targets = deepcopy(reaction_times[i])
#     h = H[i][Ti = -0.25u"s" .. 0.25u"s"]

#     h = reshape(h, (prod(size(h)[1:2]), size(h, Trial))) # Flattened data, _ × trial

#     N = fit(Normalization.MixedZScore, h; dims = 2)
#     h = normalize(h, N)
#     idxs = .!isnan.(targets) .& (targets .< 1) .* (targets .> 0.2)
#     h = h[:, idxs]
#     targets = targets[idxs]

#     HuberRegressor = MLJ.@load HuberRegressor pkg=MLJScikitLearnInterface
#     h = DataFrame(h', Symbol.(1:size(h, 1)))
#     mach = MLJ.machine(HuberRegressor(; max_iter = 10000, alpha = 1e6), h, targets)
#     MLJ.fit!(mach)
#     yhat = MLJ.predict(mach, h)
#     rho = cor(targets, yhat)
#     # end
# end

# ziggurat(rhos)

# begin # robust test
#     m = fit(RobustLinearModel, randn(100, 20), randn(100), MEstimator{HuberLoss}())
# end
