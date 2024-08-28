using Random
using MultivariateStats
using HypothesisTests
using RobustModels
using Distributed
import Bootstrap
# import LIBSVM

if !haskey(ENV, "DRWATSON_STOREPATCH")
    ENV["DRWATSON_STOREPATCH"] = "true"
end

function _preamble()
    quote
        using Statistics
        using StatsBase
        using MultivariateStats
        using Random
        using LinearAlgebra
        using Unitful
        using FileIO
        import AllenNeuropixels as AN
        using ComplexityMeasures
        using Clustering
        using CairoMakie
        using DataFrames
        using DimensionalData
        using DrWatson
        using DSP
        using Foresight
        using IntervalSets
        using Normalization
        using Distributed
        using JLD2
        using JSON
        using HypothesisTests
        using MultipleTesting
        using TimeseriesFeatures
        using TimeseriesTools
        using TimeseriesTools.Unitful
        import TimeseriesTools.TimeSeries
        using UnPack
        import Foresight.clip
        import CairoMakie.save
        import DimensionalData: metadata
        using Term
        using Term.Progress
        import AllenNeuropixels: Chan, Unit, Depth, Log𝑓
    end
end
macro preamble()
    _preamble()
end
@preamble
const THETA = (5, 10)
const GAMMA = (30, 100)
const INTERVAL = -0.25u"s" .. 0.75u"s"
const structures = ["VISp", "VISl", "VISrl", "VISal", "VISpm", "VISam"]

DimensionalData.@dim Sesh ToolsDim "Session"
DimensionalData.@dim Trial ToolsDim "Trial"
export Sesh, Trial

function calcquality(dirname; suffix = "jld2", connector = connector)
    files = readdir(dirname)
    ps = []
    for f in files
        f = joinpath(dirname, f)
        _, parameters, _suffix = parse_savename(f; connector)
        if _suffix == suffix
            try
                fl = jldopen(f)
                if haskey(fl, "error")
                    continue
                else
                    push!(ps, parameters)
                end
            catch e
                @warn e
                continue
            end
        end
    end
    ks = keys.(ps) |> collect
    vs = values.(ps) .|> collect
    vs = stack(vs)
    dims = unique(stack(ks))
    uvs = unique.(eachrow(vs))
    si = findfirst(dims .== ["stimulus"])
    if !isnothing(si)
        vs = replace(vs, "Natural_Images" => r"Natural_Images")
        uvs[si] = replace(uvs[si], "Natural_Images" => r"Natural_Images")
    end
    ddims = Tuple([Dim{Symbol(d)}(s) for (s, d) in zip(uvs, dims)])

    Q = Array{Bool}(undef, length.(uvs)...)
    Q = ToolsArray(Q, ddims)
    Q .= false
    for v in eachcol(vs)
        ds = [Dim{Symbol(d)}(At(s)) for (s, d) in zip(v, dims)]
        Q[ds...] = true
    end
    if any(isa.(DimensionalData.dims(Q), (Dim{:structure},))) &&
       all(lookup(Q, :structure) .∈ [structures])
        Q = Q[structure = At(structures)] # * Sort to global structures order
    end
    return Q
end

function commondepths(depths)
    # Find a Range of depths that best approximates the given collection of depths
    # * Get the smallest possible bound
    I = ceil(median(minimum.(depths)), sigdigits = 2) ..
        floor(median(maximum.(depths)), sigdigits = 2)
    depths = map(depths) do d
        d[d .∈ (I,)]
    end
    N = minimum(length.(depths)) # * We never want to sample a depth twice
    guess = range(start = I.left, stop = I.right,
                  length = 20)
end

function parselayernum(layername)
    m = match(r"\d+", layername)
    if m !== nothing
        m = parse(Int, m.match)
    else
        m = 0
    end
    if m > 2
        m -= 1 # Merge layer 2 and 3
    end
    return m
end

function layernum2name(num)
    if num < 2
        return string(num)
    elseif num == 2
        return "2/3"
    else
        return string(num + 1)
    end
end

function isbad(outfile; retry_errors = true)
    if !isfile(outfile)
        return true
    end
    if retry_errors
        f = jldopen(outfile)
        ind = haskey(f, "error")
        close(f)
        if ind
            return true
        else
            return false
        end
    else
        @info "File $outfile exists but contains an error. It will not be overwritten as $retry_errors is `false`."
        return false
    end
end

function on_error(e)
    @warn e
    return false
end

function val_to_string(v)
    if v isa Regex
        return v.pattern
    else
        return string(v)
    end
end

allowedtypes = (Real, String, Regex, Symbol, TimeType, Vector, Tuple)

symextrema(x) = (m = maximum(abs.(extrema(x))); (-m, m))

function savepath(D::Dict, ext = "", args...)
    filename = savename(D, ext; connector, val_to_string, allowedtypes)
    return joinpath(args..., filename)
end
function savepath(prefix::String, D::Dict, ext = "", args...)
    filename = savename(prefix, D, ext; connector, val_to_string, allowedtypes)
    return joinpath(args..., filename)
end

function unique_inverse(A::AbstractArray)
    out = Array{eltype(A)}(undef, 0)
    out_idx = Array{Vector{Int}}(undef, 0)
    seen = Dict{eltype(A), Int}()
    for (idx, x) in enumerate(A)
        if !in(x, keys(seen))
            seen[x] = length(seen) + 1
            push!(out, x)
            push!(out_idx, Int[])
        end
        push!(out_idx[seen[x]], idx)
    end
    out, out_idx
end

function crossvalidate(n::Integer, k, rng::AbstractRNG = Random.default_rng())
    idxs = 1:n
    idxs = shuffle(rng, idxs)
    N = Int(round(n / k))
    test = map(1:k) do i
        if i == k
            view(idxs, (1 + (i - 1) * N):n)
        else
            view(idxs, (1 + (i - 1) * N):min(i * N, n))
        end
    end
    train = setdiff.([idxs], test)
    return train, test
end

"""
Stratified k-fold cross validation
"""
function crossvalidate(classlabels::AbstractVector, k,
                       rng::AbstractRNG = Random.default_rng())
    labels = unique(classlabels)
    idxs = findall.(map(x -> x .== classlabels, unique(labels)))
    @assert minimum(length.(idxs)) ≥ k # Otherwise one or more folds contain only one class
    idxs = map(x -> shuffle!(rng, x), idxs)
    coeffs = length.(idxs) ./ k
    test = map(1:k) do s
        out = map(idxs, coeffs) do i, c
            a, b = round.(Integer, [s - 1, s] .* c)
            view(i, (a + 1):b)
        end
        vcat(out...)
    end
    train = setdiff.([eachindex(classlabels)], test)
    return train, test
end

function classification_metrics(y::Vector{<:Bool}, ŷ::Vector{<:Bool})
    @assert eltype(y) == eltype(ŷ) == Bool
    Np = sum(y)
    Nn = sum(.!y)
    tp = y .& ŷ |> sum
    fn = y .& .!ŷ |> sum
    tn = .!y .& .!ŷ |> sum
    fp = .!y .& ŷ |> sum
    sensitivity = tp / (tp + fn)
    specificity = tn / (tn + fp)
    acc = (tp + tn) / (tp + tn + fp + fn)
    bacc = (sensitivity + specificity) / 2
    D = Dict{Symbol, Any}()
    @pack! D = Np, Nn, tp, fn, tn, fp, sensitivity, specificity, acc, bacc
    return DataFrame(D)
end
balanced_accuracy(args...) = classification_metrics(args...).bacc |> only
accuracy(args...) = classification_metrics(args...).acc |> only

function classifier(h, labels; regcoef = 0.5)
    N = fit(Normalization.ZScore, h; dims = 2)
    h = normalize(h, N)

    # * LDA
    M = fit(MulticlassLDA, collect(h), collect(labels);
            regcoef = convert(eltype(h), regcoef))
    ŷ = (ptiredict(M, h) .> 0) |> vec |> collect # Predicted output classes

    # # * SVM
    # M = LIBSVM.svmtrain(h, labels, cost = Float64(regcoef))
    # ŷ = LIBSVM.svmpredict(M, h) |> first

    if cor(labels, ŷ) < 0
        M.proj = .-M.proj
    end
    # if cor(labels, ŷ) < 0
    #     @warn "Correlation is negative: $(cor(labels, ŷ))"
    # end
    return N, M
end

function classifier(H; dim = :trial, regcoef = 0.1)
    labels = lookup(H, dim)
    @assert eltype(labels) == Bool
    negdims = setdiff(1:ndims(H), [dimnum(H, dim)])
    negsize = size(H)[negdims]
    h = reshape(H, (prod(negsize), size(H, dim))) # Flattened data, _ × trial
    classifier(h, labels; regcoef)
end

function classify_kfold(H, rng::AbstractRNG = Random.default_rng(); k = 5, dim = :trial,
                        repeats = 10,
                        kwargs...)
    labels = lookup(H, dim)
    @assert eltype(labels) == Bool
    negdims = setdiff(1:ndims(H), [dimnum(H, dim)])
    negsize = size(H)[negdims]
    h = reshape(H, (prod(negsize), size(H, dim))) # Flattened data, _ × trial

    bacs = map(1:repeats) do _
        trains, tests = crossvalidate(labels, k, rng)

        bac = map(trains, tests) do train, test
            try
                N, M = classifier(h[:, train], labels[train]; kwargs...)
                y = labels[test] |> collect
                ŷ = predict(M, normalize(h[:, test], N)) .> 0 # Predicted output classes
                # ŷ = LIBSVM.svmpredict(M, normalize(h[:, test], N)) |> first # Predicted output classes
                ŷ = ŷ |> vec |> collect
                bac = balanced_accuracy(y, ŷ)
                return bac
            catch e
                @warn e
                return NaN
            end
        end
        return mean(bac)
    end
    return mean(bacs)
end

# mutable struct Regressor{T}
#     a::Vector{T}
#     b::T
# end
# function (M::Regressor{T})(x::AbstractMatrix{<:T})::Vector{T} where {T}
#     return M.a' * x .+ M.b |> vec |> collect
# end

# function regressor(h::AbstractMatrix, targets; regcoef = 0.1)
#     N = fit(Normalization.MixedZScore, h; dims = 2)
#     h = normalize(h, N)
#     # a..., b = ridge(h', convert.(eltype(h), targets), regcoef)
#     idxs = 0.2 .< targets .< 1 # Remove outliers
#     # idxs = trues(length(targets))

#     # # m = fit(LassoPath, Float64.(h[:, idxs]'), Float64.(targets[idxs]))
#     # m = fit(RobustLinearModel, Float64.(h[:, idxs]'), Float64.(targets[idxs]),
#     #         MEstimator{HuberLoss}(); σ0 = 2)

#     HuberRegressor = MLJ.@load HuberRegressor verbosity=0 pkg=MLJScikitLearnInterface
#     h = DataFrame(h[:, idxs]', Symbol.(1:size(h, 1)))
#     mach = MLJ.machine(HuberRegressor(; max_iter = 10000, alpha = regcoef), h,
#                        targets[idxs])
#     MLJ.fit!(mach; verbosity = 0)

#     # M = Regressor(a, b)
#     M = x -> MLJ.predict(mach, DataFrame(x', Symbol.(1:size(x, 1))))
#     return N, M
# end

# function regressor(H::AbstractArray{T, 3}, targets; dim = :trial, regcoef = 0.1) where {T}
#     negdims = setdiff(1:ndims(H), [dimnum(H, dim)])
#     negsize = size(H)[negdims]
#     h = reshape(H, (prod(negsize), size(H, dim))) # Flattened data, _ × trial
#     regressor(h, targets; regcoef)
# end

# function regress_kfold(H, reaction_times, rng::AbstractRNG = Random.default_rng(); k = 5,
#                        dim = :trial,
#                        kwargs...)
#     labels = lookup(H, dim)
#     @assert eltype(labels) == Bool
#     negdims = setdiff(1:ndims(H), [dimnum(H, dim)])
#     negsize = size(H)[negdims]
#     h = reshape(H, (prod(negsize), size(H, dim))) # Flattened data, _ × trial

#     h = h[:, .!isnan.(reaction_times)]
#     reaction_times = reaction_times[.!isnan.(reaction_times)]

#     trains, tests = crossvalidate(length(reaction_times), k, rng)

#     ρ = map(trains, tests) do train, test
#         try
#             N, M = regressor(h[:, train], reaction_times[train]; kwargs...)
#             y = reaction_times[test] |> collect
#             ŷ = M(normalize(h[:, test], N))[:]
#             ρ = corspearman(y, ŷ)
#             return ρ
#         catch e
#             @warn e
#             return NaN
#         end
#     end
#     return mean(ρ) # Use fisher z transform
# end

function ppc(x::AbstractVector{T})::T where {T} # Eq. 14 of Vinck 2010
    isempty(x) && return NaN
    N = length(x)
    Δ = zeros(N - 1)
    Threads.@threads for i in 1:(N - 1)
        δ = @views x[i] .- x[(i + 1):end]
        Δ[i] = sum(cos.(δ))
    end
    return (2 / (N * (N - 1))) * sum(Δ)
end
function ppc(ϕ::UnivariateTimeSeries{T}, spikes::AbstractVector)::NTuple{3, T} where {T}
    spikes = spikes[spikes .∈ [Interval(ϕ)]]
    isempty(spikes) && return (NaN, NaN, NaN)
    phis = ϕ[𝑡(Near(spikes))] |> parent
    γ = ppc(phis)
    𝑝 = isempty(phis) ? 1.0 : HypothesisTests.pvalue(RayleighTest(phis))
    p = phis |> resultant |> angle
    return (γ, p, 𝑝)
end
function ppc(ϕ::AbstractVector{<:UnivariateTimeSeries}, spikes::AbstractVector)
    γs = ppc.(ϕ, [spikes])
    return first.(γs), getindex.(γs, 2), last.(γs)
end

function initialize_spc_dataframe!(spikes, T)
    if !("trial_pairwise_phase_consistency" ∈ names(spikes))
        @debug "Initializing missing columns"
        spikes[!, :trial_pairwise_phase_consistency] = [Vector{T}()
                                                        for _ in 1:size(spikes, 1)]
    end
    if !("trial_pairwise_phase_consistency_pvalue" ∈ names(spikes))
        @debug "Initializing missing columns"
        spikes[!, :trial_pairwise_phase_consistency_pvalue] = [Vector{T}()
                                                               for _ in 1:size(spikes,
                                                                               1)]
    end
    if !("trial_pairwise_phase_consistency_angle" ∈ names(spikes))
        @debug "Initializing missing columns"
        spikes[!, :trial_pairwise_phase_consistency_angle] = [Vector{T}()
                                                              for _ in 1:size(spikes,
                                                                              1)]
    end
    if !("pairwise_phase_consistency" ∈ names(spikes))
        @debug "Initializing missing columns"
        spikes[!, :pairwise_phase_consistency] .= NaN
    end
    if !("pairwise_phase_consistency_pvalue" ∈ names(spikes))
        @debug "Initializing missing columns"
        spikes[!, :pairwise_phase_consistency_pvalue] .= NaN
    end
    if !("pairwise_phase_consistency_angle" ∈ names(spikes))
        @debug "Initializing missing columns"
        spikes[!, :pairwise_phase_consistency_angle] .= NaN
    end
end

function spc!(spikes::AbstractDataFrame, ϕ::AbstractTimeSeries; pbar = nothing) # Mutate the dataframe
    T = eltype(lookup(ϕ, 𝑡) |> ustripall)
    initialize_spc_dataframe!(spikes, T)

    probeid = metadata(ϕ)[:probeid]

    !isnothing(pbar) &&
        (job = addjob!(pbar; N = size(spikes, 1), transient = true, description = "Unit"))
    for unit in eachrow(spikes)
        if unit.probe_id == probeid
            unitid = unit.ecephys_unit_id
            _ϕ = ϕ[Depth(Near(unit.streamlinedepth))]
            spiketimes = unit.spiketimes

            _ϕ = map(eachslice(_ϕ, dims = 2)) do x # Individual trial
                x = set(x, 𝑡 => lookup(x, 𝑡) .+ refdims(x, :changetime))
            end
            γ_trial, p_trial, 𝑝_trial = ppc(_ϕ, spiketimes)
            spikes.trial_pairwise_phase_consistency[spikes.ecephys_unit_id .== unitid] .= [γ_trial]
            spikes.trial_pairwise_phase_consistency_pvalue[spikes.ecephys_unit_id .== unitid] .= [𝑝_trial]
            spikes.trial_pairwise_phase_consistency_angle[spikes.ecephys_unit_id .== unitid] .= [p_trial]

            idxs = [any(s .∈ Interval.(_ϕ)) for s in spiketimes]
            spiketimes = spiketimes[idxs]
            _ϕ = cat(_ϕ..., dims = 𝑡(vcat(lookup.(_ϕ, 𝑡)...)))
            γ, p, 𝑝 = ppc(_ϕ, spiketimes)
            spikes.pairwise_phase_consistency[spikes.ecephys_unit_id .== unitid] .= γ
            spikes.pairwise_phase_consistency_pvalue[spikes.ecephys_unit_id .== unitid] .= 𝑝
            spikes.pairwise_phase_consistency_angle[spikes.ecephys_unit_id .== unitid] .= p
        end
        !isnothing(pbar) && update!(job)
    end
    !isnothing(pbar) && stop!(job)
end
# * Assumes each LFP comes from one session, from one structure, and has has the correct metadata
function spc!(spikes::AbstractDataFrame, ϕ::AbstractVector{<:AbstractTimeSeries};
              job = nothing)
    T = eltype(lookup(first(ϕ), 𝑡) |> ustripall)
    initialize_spc_dataframe!(spikes, T)

    Threads.@threads for _ϕ in ϕ
        idxs = spikes.probe_id .== metadata(_ϕ)[:probeid]
        _spikes = @views spikes[idxs, :] # * Select structure
        spc!(_spikes, ustripall(_ϕ))
        isnothing(job) || update!(job)
    end
    isnothing(job) || stop!(job)
end

function sac(r::UnivariateTimeSeries{T}, spikes::AbstractVector)::T where {T}
    spikes = spikes[spikes .∈ [Interval(r)]]
    isempty(spikes) && return NaN
    return r[𝑡(Near(spikes))] |> parent |> mean
end
sac(r::AbstractVector{<:UnivariateTimeSeries}, spikes::AbstractVector) = sac.(r, [spikes])

function initialize_sac_dataframe!(spikes, T)
    if !("trial_spike_amplitude_coupling" ∈ names(spikes))
        @debug "Initializing missing columns"
        spikes[!, :trial_spike_amplitude_coupling] = [Vector{T}()
                                                      for _ in 1:size(spikes, 1)]
    end
    if !("spike_amplitude_coupling" ∈ names(spikes))
        @debug "Initializing missing columns"
        spikes[!, :spike_amplitude_coupling] .= NaN
    end
end

function sac!(spikes::AbstractDataFrame, r::AbstractTimeSeries; pbar = nothing,
              normfunc = HalfZScore) # At the level of each probe. We normalize here
    T = eltype(lookup(r, 𝑡) |> ustripall)
    initialize_sac_dataframe!(spikes, T)

    probeid = metadata(r)[:probeid]

    !isnothing(pbar) &&
        (job = addjob!(pbar; N = size(spikes, 1), transient = true, description = "Unit"))
    for unit in eachrow(spikes)
        if unit.probe_id == probeid
            unitid = unit.ecephys_unit_id
            _r = r[Depth(Near(unit.streamlinedepth))]
            if !isnothing(normfunc) # * Normalize over trials AND time
                N = fit(normfunc, _r) # * Crucial; normalization makes this correlation-like
                _r = normalize(_r, N)
            end
            spiketimes = unit.spiketimes

            _r = map(eachslice(_r, dims = 2)) do x # Individual trial
                x = set(x, 𝑡 => lookup(x, 𝑡) .+ refdims(x, :changetime))
            end
            γ_trial = sac(_r, spiketimes)
            spikes.trial_spike_amplitude_coupling[spikes.ecephys_unit_id .== unitid] .= [γ_trial]

            idxs = [any(s .∈ Interval.(_r)) for s in spiketimes]
            spiketimes = spiketimes[idxs]
            _r = cat(_r..., dims = 𝑡(vcat(lookup.(_r, 𝑡)...)))
            γ = sac(_r, spiketimes)
            spikes.spike_amplitude_coupling[spikes.ecephys_unit_id .== unitid] .= γ
        end
        !isnothing(pbar) && update!(job)
    end
    !isnothing(pbar) && stop!(job)
end

# * Assumes each LFP comes from one session, from one structure, and has has the correct
#   metadata. Each element is one subject.
function sac!(spikes::AbstractDataFrame, r::AbstractVector{<:AbstractTimeSeries};
              job = nothing) # At the level of structures
    T = eltype(lookup(first(r), 𝑡) |> ustripall)
    initialize_sac_dataframe!(spikes, T)
    Threads.@threads for _r in r
        idxs = spikes.probe_id .== metadata(_r)[:probeid]
        _spikes = @views spikes[idxs, :] # * Select structure
        sac!(_spikes, ustripall(_r))
        isnothing(job) || update!(job)
    end
    isnothing(job) || stop!(job)
end

struct HistBins
    bints::AbstractVector{<:AbstractInterval}
    idxs::AbstractVector{AbstractVector{<:Integer}}
end
bincenters(B::HistBins) = mean.(B.bints)
(B::HistBins)(x) = ToolsArray(getindex.([x], B.idxs), (Dim{:bin}(bincenters(B)),))

function HistBins(x; bins = StatsBase.histrange(x, 10))
    if bins isa Integer
        bins = StatsBase.histrange(x, bins)
    end

    bints = map(eachindex(bins)[1:(end - 1)]) do i
        if i == length(bins) - 1
            Interval{:closed, :closed}(bins[i], bins[i + 1])
        else
            Interval{:closed, :open}(bins[i], bins[i + 1])
        end
    end
    idxs = map(bin -> findall(x .∈ [bin]), bints)
    return HistBins(bints, idxs)
end

function pac(ϕ::AbstractVector, r::AbstractVector; kwargs...)
    ϕ = mod2pi.(ϕ .+ pi) .- pi
    ModulationIndices.tort2010(ϕ, r; kwargs...)
end
function pac(ϕ::AbstractToolsArray, r::AbstractToolsArray; dims, kwargs...)
    out = similar(first(eachslice(ϕ; dims = dims))) # Template
    dims = dimnum.([ϕ], dims)
    @assert (length(dims) == (ndims(ϕ) - 1)) || length(dims) == 1
    if length(dims) > 1
        negdim = setdiff(collect(1:ndims(ϕ)), dims) |> only
        ns = size(ϕ, negdim)
        ϕ = eachslice(parent(ϕ); dims = negdim)
        ϕ = hcat([p[:] for p in ϕ]...)
        r = eachslice(parent(r); dims = negdim)
        r = hcat([a[:] for a in r]...)
        @assert ndims(ϕ) == 2
        @assert size(ϕ, 2) == ns
        dims = 1
    end
    dims = setdiff(collect(1:ndims(ϕ)), dims)
    dims = Tuple(dims)
    progressmap(eachindex(out), eachslice(ϕ; dims), eachslice(r; dims);
                parallel = true) do i, ϕ, r
        out[i] = pac(collect(ϕ), collect(r); kwargs...)
    end
    return out
end
_confidence(x, z = 1.96) = z * std(x) ./ sqrt(length(x)) # * 1.96 for 95% confidence interval
confidence(x, args...; dims = 1:ndims(x)) = mapslices(x -> _confidence(x, args...), x; dims)
function quartiles(X::AbstractArray; dims = 1)
    q1 = mapslices(x -> quantile(x, 0.25), X; dims)
    q2 = mapslices(x -> quantile(x, 0.5), X; dims)
    q3 = mapslices(x -> quantile(x, 0.75), X; dims)
end

function bootstrapaverage(average, x::AbstractVector{T}; confint = 0.95,
                          N = 10000)::Tuple{T, Tuple{T, T}} where {T}
    sum(!isnan, x) < 5 && return (NaN, (NaN, NaN))

    # * Estimate a sampling distribution of the average
    x = filter(!isnan, x)
    b = Bootstrap.bootstrap(nansafe(average), x, Bootstrap.BalancedSampling(N))
    μ, σ... = only(Bootstrap.confint(b, Bootstrap.BCaConfInt(confint)))
    return μ, σ
end

function bootstrapaverage(average, X::AbstractArray; dims = 1, kwargs...)
    ds = [i == dims ? 1 : Colon() for i in 1:ndims(X)]
    μ = similar(X[ds...])
    σl = similar(μ)
    σh = similar(μ)
    negdims = filter(!=(dims), 1:ndims(X)) |> Tuple
    Threads.@threads for (i, x) in collect(enumerate(eachslice(X; dims = negdims)))
        μ[i], (σl[i], σh[i]) = bootstrapaverage(average, x; kwargs...)
    end
    return μ, (σl, σh)
end
function bootstrapaverage(average, X::AbstractToolsArray; dims = 1, kwargs...)
    dims = dimnum(X, dims)
    ds = [i == dims ? 1 : Colon() for i in 1:ndims(X)]
    μ = similar(X[ds...])
    σl = similar(μ)
    σh = similar(μ)
    negdims = filter(!=(dims), 1:ndims(X)) |> Tuple
    Threads.@threads for (i, x) in collect(enumerate(eachslice(X; dims = negdims)))
        μ[i], (σl[i], σh[i]) = bootstrapaverage(average, parent(x); kwargs...)
    end
    return μ, (σl, σh)
end
bootstrapmedian(args...; kwargs...) = bootstrapaverage(median, args...; kwargs...)
