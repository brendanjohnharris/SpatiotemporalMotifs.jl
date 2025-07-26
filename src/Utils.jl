using Random
using MultivariateStats
using HypothesisTests
using RobustModels
using Distributed
import Bootstrap
using MultipleTesting
using Optim
using IntervalSets
using USydClusters

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
        import SpatiotemporalMotifs.calcdir
        using TerminalLoggers, Logging
        global_logger(TerminalLogger(right_justify = 200))
        using ProgressLogging
    end
end
macro preamble()
    _preamble()
end

function _DIR()
    dir = ""
    if (first(THETA()) != 3) || (last(THETA()) != 10)
        dir *= "&THETA=$(THETA())"
    end
    if (first(GAMMA()) != 30) || (last(GAMMA()) != 100)
        dir *= "&GAMMA=$(GAMMA())"
    end
    return dir
end

calcdir(args...; kwargs...) = projectdir("data" * _DIR(), args...; kwargs...)
plotdir(args...; kwargs...) = projectdir("plots" * _DIR(), args...; kwargs...)
export calcdir, plotdir

THETA() = haskey(ENV, "SM_THETA") ? eval(Meta.parse(ENV["SM_THETA"])) : (3, 10)
GAMMA() = haskey(ENV, "SM_GAMMA") ? eval(Meta.parse(ENV["SM_GAMMA"])) : (30, 100)
const INTERVAL = -0.25u"s" .. 0.75u"s"
const structures = ["VISp", "VISl", "VISrl", "VISal", "VISpm", "VISam"]
const PTHR = 1e-2
@preamble

DimensionalData.@dim SessionID ToolsDim "SessionID"
DimensionalData.@dim Trial ToolsDim "Trial"
DimensionalData.@dim Structure ToolsDim "Structure"
const Freq = TimeseriesTools.𝑓
export SessionID, Trial, Structure, Freq

check_keys(D, required_keys) = all(haskey.([D], required_keys))
function check_calc_keys(D)
    required_keys = [
        "trials",
        "channels",
        "streamlinedepths",
        "units",
        "spiketimes",
        "layerinfo",
        "pass_θ",
        "pass_γ",
        "performance_metrics"
    ]
    return check_keys(D, required_keys)
end

function tmap(f, d::DimensionalData.Dimension; kwargs...)
    ToolsArray(map(f, parent(d); kwargs...), (d,))
end

"""
Check the quality of the collected calculations files e.g. `out&stimulus=...`
"""
function calcquality(; path = calcdir(), fallbackpath = calcdir("calculations"),
                     prefix = "out", suffix = "jld2",
                     connector = connector)
    _Q = calcquality(fallbackpath; suffix, connector, require = false)
    if !all(_Q)
        error("Some calculation files are missing or have errors. Please rerun `scripts/calculations/calculations.jl`.") |>
        throw
    end

    files = readdir(path)
    files = filter(Base.Fix1(==, "." * suffix) ∘ last ∘ Base.splitext, files)
    files = filter(Base.Fix2(==, prefix * connector) ∘
                   Base.Fix2(getindex, 1:length(prefix * connector)), files)

    vars = [
        "pass_θ",
        "pass_γ",
        "performance_metrics",
        "trials",
        "streamlinedepths",
        "layernames",
        "spiketimes"
    ]

    Q = map(files) do f
        stimulus = parse_savename(f; connector)[2]["stimulus"]
        if stimulus == "Natural_Images"
            stimulus = r"Natural_Images"
        end
        f = joinpath(path, f)
        Q = jldopen(f, "r"; iotype = IOStream) do fl
            grp = fl.root_group
            areas = keys(grp)
            areas = areas[areas .!== ["_types"]]
            tmap(Structure(areas)) do area
                sessions = keys(grp[area])
                tmap(SessionID(Meta.parse.(sessions))) do session
                    if contains(stimulus |> string, "nochange")
                        subvars = sort([:csd])
                    else
                        subvars = sort([:V, :csd, :θ, :ϕ, :r, :k, :ω])
                    end
                    # * Check all subvars exist
                    subkeys = keys(grp[area][string(session)])
                    q = setdiff(vars, subkeys)
                    if !contains(stimulus |> string, "omission") && !isempty(q)
                        @debug "Missing variables in `$f`: $q"
                        return false
                    end
                    q = setdiff(string.(subvars), subkeys)
                    if !isempty(q)
                        @debug "Missing variables in `$f`: $q"
                        return false
                    end
                    return true
                end
            end |> stack
        end |> stack
        return Q, stimulus
    end
    stimuli = last.(Q)
    Q = first.(Q)
    badstimuli = [any(.!q) for q in Q]
    if any(badstimuli)
        # map(findall(badstimuli)) do stimulus
        #     collect_calculations(_Q; path = fallbackpath,
        #                          stimulus,
        #                          rewrite = false,
        #                          outpath = path)
        # end
        # Q = calcquality(; path, fallbackpath,
        #                 prefix, suffix,
        #                 connector)
        @warn """Some collected calculations are missing or have errors for stimuli: `$(join(stimuli[badstimuli], "; "))`"""
    end
    return ToolsArray(Q, (Dim{:stimulus}(stimuli),)) |> stack
end

"""
Check the quality of a a calculations directory e.g. data/calculations/`
"""
function calcquality(dirname; suffix = "jld2", connector = connector, require = true)
    if !(require isa Bool)
        _require = require
        require = true
    else
        _require = []
    end
    @info "Checking quality of calculations in `$dirname`"
    files = readdir(dirname)
    threadlog = Threads.Atomic{Int}(0)
    threadmax = length(files)
    lk = Threads.ReentrantLock()
    ps = []
    @withprogress name="Checking quality" begin
        Threads.@threads for f in files
            f = joinpath(dirname, f)
            _, parameters, _suffix = parse_savename(f; connector)
            if _suffix == suffix
                try
                    if require
                        cannotload = jldopen(f, "r"; iotype = IOStream) do fl
                            # fl["performance_metrics"] # Can load
                            haskey(fl, "error") || !check_calc_keys(fl) ||
                                !check_keys(fl, string.(_require))
                        end
                    else
                        cannotload = false
                    end
                    if cannotload
                        continue
                    end
                catch e
                    @warn e
                    continue
                end
                lock(lk) do
                    push!(ps, parameters) # Only add to list of good files if file exists and has no error
                    if require
                        if Threads.threadid() ∈ 1:2
                            Threads.atomic_add!(threadlog, Threads.nthreads() ÷ 2)
                            @logprogress threadlog[] / threadmax
                        end
                    end
                end
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

    function map2dims(d)
        if d == "sessionid"
            return SessionID
        elseif d == "trial"
            return Trial
        elseif d == "structure"
            return Structure
        else
            return Dim{Symbol(d)}
        end
    end
    ddims = Tuple([map2dims(d)(s) for (s, d) in zip(uvs, dims)])

    Q = falses(length.(uvs)...)
    Q = ToolsArray(Q, ddims)
    for v in eachcol(vs)
        ds = [map2dims(d)(At(s)) for (s, d) in zip(v, dims)]
        Q[ds...] = true
    end
    if any(isa.(DimensionalData.dims(Q), (Structure,))) &&
       all(lookup(Q, Structure) .∈ [structures])
        Q = Q[Structure = At(structures)] # * Sort to global structures order
    end
    @info "Mean quality: $(mean(Q))"
    return Q
end

function submit_calculations(exprs; queue = ``, mem = 50)
    exprs = deepcopy(exprs)
    if length(exprs) > 2
        N3 = (length(exprs) ÷ 3)
        shuffle!(exprs) # ? Shuffle so restarted calcs are more even
        if isempty(queue)
            USydClusters.Physics.runscripts(exprs[1:N3]; ncpus = 8, mem,
                                            walltime = 8,
                                            project = projectdir(), exeflags = `+1.10.10`,
                                            queue = `h100`)
            USydClusters.Physics.runscripts(exprs[(N3 + 1):(N3 * 2)]; ncpus = 8, mem,
                                            walltime = 8,
                                            project = projectdir(), exeflags = `+1.10.10`,
                                            queue = `l40s`)
            USydClusters.Physics.runscripts(exprs[(2 * N3 + 1):end]; ncpus = 8, mem,
                                            walltime = 8,
                                            project = projectdir(), exeflags = `+1.10.10`,
                                            queue = `taiji`)
        else
            USydClusters.Physics.runscripts(exprs; ncpus = 8, mem,
                                            walltime = 8,
                                            project = projectdir(), exeflags = `+1.10.10`,
                                            queue = queue)
        end
    else
        USydClusters.Physics.runscript.(exprs; ncpus = 8, mem,
                                        walltime = 8,
                                        project = projectdir(), exeflags = `+1.10.10`,
                                        queue = `h100`)
    end
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

function coordinates(X::AbstractDimArray)
    cs = Iterators.product(lookup(X)...)
    cs = [getindex.(cs, i) for i in 1:ndims(X)]
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

function isbad(outfile; retry_errors = true, check_other_file = false)
    if !isfile(outfile)
        return true
    end
    if isnothing(check_other_file)
        check_other_file = true
    else
        check_other_file = isfile(check_other_file)
    end
    if retry_errors
        try
            jldopen(outfile, "r") do f
                ind = haskey(f, "error")
                if ind
                    if contains(f["error"], "Region error")
                        return false
                    else
                        return check_other_file
                    end
                else
                    return false
                end
            end
        catch e
            @warn e
            @warn "Could not open $outfile"
            return true
        end
    else
        @info "File `$outfile` exists but contains an error. It will not be overwritten as $retry_errors is `false`."
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

function savepath(D::Union{Dict, NamedTuple}, ext = "", args...)
    filename = savename(D, ext; connector, val_to_string, allowedtypes)
    return joinpath(args..., filename)
end
function savepath(prefix::String, D::Union{Dict, NamedTuple}, ext = "", args...)
    filename = savename(prefix, D, ext; connector, val_to_string, allowedtypes)
    return joinpath(args..., filename)
end
savepath(prefix::String) = (args...) -> savepath(prefix, args...)

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
    ŷ = (predict(M, h) .> 0) |> vec |> collect # Predicted output classes

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

function classifier(H; dim = Trial, regcoef = 0.1)
    labels = lookup(H, dim)
    @assert eltype(labels) == Bool
    negdims = setdiff(1:ndims(H), [dimnum(H, dim)])
    negsize = size(H)[negdims]
    h = reshape(H, (prod(negsize), size(H, dim))) # Flattened data, _ × trial
    classifier(h, labels; regcoef)
end

function classify_kfold(H, rng::AbstractRNG = Random.default_rng(); k = 5, dim = Trial,
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

function tuneclassifier(H; r0 = 0.5, repeats = 10, k = 5, kwargs...)
    r0 = log10(r0)
    objective = r -> -classify_kfold(H[𝑡 = -0.25u"s" .. 0.25u"s"]; # Just pre-offset
                                     regcoef = exp10(only(r)), repeats, k,
                                     kwargs...)
    o = optimize(objective, [r0], ParticleSwarm(),
                 Optim.Options(; iterations = 10))
    return exp10(only(o.minimizer)), -o.minimum
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

# function regressor(H::AbstractArray{T, 3}, targets; dim = Trial, regcoef = 0.1) where {T}
#     negdims = setdiff(1:ndims(H), [dimnum(H, dim)])
#     negsize = size(H)[negdims]
#     h = reshape(H, (prod(negsize), size(H, dim))) # Flattened data, _ × trial
#     regressor(h, targets; regcoef)
# end

# function regress_kfold(H, reaction_times, rng::AbstractRNG = Random.default_rng(); k = 5,
#                        dim = Trial,
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
    spikes = spikes[spikes .∈ [Interval(ϕ)]] # * Important: we only want spikes that are within the time range of the ϕ
    isempty(spikes) && return (NaN, NaN, NaN)
    phis = ϕ[𝑡(Near(spikes))] |> parent
    idxs = .!isnan.(phis)
    sum(idxs) == 0 && return (NaN, NaN, NaN)
    phis = phis[idxs]
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
    trial_ppcs = [:trial_pairwise_phase_consistency,
        :trial_pairwise_phase_consistency_pvalue,
        :trial_pairwise_phase_consistency_angle]

    map(trial_ppcs) do col
        if !(string(col) ∈ names(spikes))
            spikes[!, col] = [Vector{T}() for _ in 1:size(spikes, 1)]
        end
    end

    ppcs = [:pairwise_phase_consistency,
        :pairwise_phase_consistency_pvalue,
        :pairwise_phase_consistency_angle,
        :onset_pairwise_phase_consistency,
        :onset_pairwise_phase_consistency_pvalue,
        :onset_pairwise_phase_consistency_angle,
        :offset_pairwise_phase_consistency,
        :offset_pairwise_phase_consistency_pvalue,
        :offset_pairwise_phase_consistency_angle,
        :hit_onset_pairwise_phase_consistency,
        :hit_onset_pairwise_phase_consistency_pvalue,
        :hit_onset_pairwise_phase_consistency_angle,
        :hit_offset_pairwise_phase_consistency,
        :hit_offset_pairwise_phase_consistency_pvalue,
        :hit_offset_pairwise_phase_consistency_angle,
        :miss_onset_pairwise_phase_consistency,
        :miss_onset_pairwise_phase_consistency_pvalue,
        :miss_onset_pairwise_phase_consistency_angle,
        :miss_offset_pairwise_phase_consistency,
        :miss_offset_pairwise_phase_consistency_pvalue,
        :miss_offset_pairwise_phase_consistency_angle
    ]

    map(ppcs) do col
        if !(string(col) ∈ names(spikes))
            spikes[!, col] .= NaN
        end
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
            ϕd = ϕ[Depth(Near(unit.streamlinedepth))]
            spiketimes = unit.spiketimes
            thisunit = spikes.ecephys_unit_id .== unitid
            hitmiss = spikes[spikes.ecephys_unit_id .== unitid, :hitmiss] |> unique |> only

            begin # * Trial PPC
                _ϕ = map(eachslice(ϕd, dims = 2)) do x # Individual trial
                    toffset = refdims(x, :changetime)
                    @assert !isempty(toffset)
                    x = set(x, 𝑡 => lookup(x, 𝑡) .+ toffset)
                end
                γ_trial, p_trial, 𝑝_trial = ppc(_ϕ, spiketimes)
                spikes.trial_pairwise_phase_consistency[thisunit] .= [γ_trial]
                spikes.trial_pairwise_phase_consistency_pvalue[thisunit] .= [𝑝_trial]
                spikes.trial_pairwise_phase_consistency_angle[thisunit] .= [p_trial]
            end
            begin # * PPC for onset periods
                Δt = 0 .. 0.25 # s Relative to changetime
                @assert Δt ⊆ Interval(ϕd)
                __ϕ = map(eachslice(ϕd[𝑡 = Δt], dims = 2)) do x # Individual trial
                    toffset = refdims(x, :changetime)
                    @assert !isempty(toffset)
                    x = set(x, 𝑡 => lookup(x, 𝑡) .+ toffset)
                end

                _idxs = [any(s .∈ Interval.(__ϕ)) for s in spiketimes]
                _spiketimes = spiketimes[_idxs]

                # * First, hit/miss
                ___ϕ = deepcopy(__ϕ)[hitmiss]
                ___ϕ = cat(___ϕ..., dims = 𝑡(vcat(lookup.(___ϕ, 𝑡)...)))
                γ, p, 𝑝 = ppc(___ϕ, _spiketimes)
                spikes.hit_onset_pairwise_phase_consistency[thisunit] .= γ
                spikes.hit_onset_pairwise_phase_consistency_pvalue[thisunit] .= 𝑝
                spikes.hit_onset_pairwise_phase_consistency_angle[thisunit] .= p

                ___ϕ = deepcopy(__ϕ)[.!hitmiss]
                ___ϕ = cat(___ϕ..., dims = 𝑡(vcat(lookup.(___ϕ, 𝑡)...)))
                γ, p, 𝑝 = ppc(___ϕ, _spiketimes)
                spikes.miss_onset_pairwise_phase_consistency[thisunit] .= γ
                spikes.miss_onset_pairwise_phase_consistency_pvalue[thisunit] .= 𝑝
                spikes.miss_onset_pairwise_phase_consistency_angle[thisunit] .= p

                # * Then all trials
                __ϕ = cat(__ϕ..., dims = 𝑡(vcat(lookup.(__ϕ, 𝑡)...)))
                γ, p, 𝑝 = ppc(__ϕ, _spiketimes)
                spikes.onset_pairwise_phase_consistency[thisunit] .= γ
                spikes.onset_pairwise_phase_consistency_pvalue[thisunit] .= 𝑝
                spikes.onset_pairwise_phase_consistency_angle[thisunit] .= p
            end
            begin # * PPC for offset periods
                Δt = 0.25 .. 0.5 # s Relative to changetime
                @assert Δt ⊆ Interval(ϕd)

                __ϕ = map(eachslice(ϕd[𝑡 = Δt], dims = 2)) do x # Individual trial
                    toffset = refdims(x, :changetime)
                    @assert !isempty(toffset)
                    x = set(x, 𝑡 => lookup(x, 𝑡) .+ toffset)
                end

                # * First hit/miss
                ___ϕ = deepcopy(__ϕ)[hitmiss]
                ___ϕ = cat(___ϕ..., dims = 𝑡(vcat(lookup.(___ϕ, 𝑡)...)))
                γ, p, 𝑝 = ppc(___ϕ, _spiketimes)
                spikes.hit_offset_pairwise_phase_consistency[thisunit] .= γ
                spikes.hit_offset_pairwise_phase_consistency_pvalue[thisunit] .= 𝑝
                spikes.hit_offset_pairwise_phase_consistency_angle[thisunit] .= p

                ___ϕ = deepcopy(__ϕ)[.!hitmiss]
                ___ϕ = cat(___ϕ..., dims = 𝑡(vcat(lookup.(___ϕ, 𝑡)...)))
                γ, p, 𝑝 = ppc(___ϕ, _spiketimes)
                spikes.miss_offset_pairwise_phase_consistency[thisunit] .= γ
                spikes.miss_offset_pairwise_phase_consistency_pvalue[thisunit] .= 𝑝
                spikes.miss_offset_pairwise_phase_consistency_angle[thisunit] .= p

                # * Then all trials
                _idxs = [any(s .∈ Interval.(__ϕ)) for s in spiketimes]
                _spiketimes = spiketimes[_idxs]
                __ϕ = cat(__ϕ..., dims = 𝑡(vcat(lookup.(__ϕ, 𝑡)...)))
                γ, p, 𝑝 = ppc(__ϕ, _spiketimes)

                spikes.offset_pairwise_phase_consistency[thisunit] .= γ
                spikes.offset_pairwise_phase_consistency_pvalue[thisunit] .= 𝑝
                spikes.offset_pairwise_phase_consistency_angle[thisunit] .= p
            end
            # * Cat all trials
            idxs = [any(s .∈ Interval.(_ϕ)) for s in spiketimes]
            spiketimes = spiketimes[idxs]
            _ϕ = cat(_ϕ..., dims = 𝑡(vcat(lookup.(_ϕ, 𝑡)...)))
            γ, p, 𝑝 = ppc(_ϕ, spiketimes)
            spikes.pairwise_phase_consistency[thisunit] .= γ
            spikes.pairwise_phase_consistency_pvalue[thisunit] .= 𝑝
            spikes.pairwise_phase_consistency_angle[thisunit] .= p
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
    x = r[𝑡(Near(spikes))] |> parent
    return mean(x[.!isnan.(x)])
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
            _r = r[Depth(Near(unit.streamlinedepth))] # * We have selected a single depth, so...
            if !isnothing(normfunc) # * ... normalize over trials AND time
                N = fit(normfunc, _r) # * Crucial; normalization makes this correlation-like. dim 1 is time
                _r = normalize(_r, N)
            end
            spiketimes = unit.spiketimes

            _r = map(eachslice(_r, dims = 2)) do x # Individual trial
                toffset = refdims(x, :changetime)
                @assert !isempty(toffset)
                x = set(x, 𝑡 => lookup(x, 𝑡) .+ toffset)
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
    sum(!isnan, x; init = 0) < 5 && return (NaN, (NaN, NaN))

    # * Estimate a sampling distribution of the average
    x = filter(!isnan, x)
    b = Bootstrap.bootstrap(nansafe(average), x, Bootstrap.BalancedSampling(N))
    μ, σ... = only(Bootstrap.confint(b, Bootstrap.BCaConfInt(confint)))
    return μ, σ
end

function bootstrapaverage(average, X::AbstractArray; dims = 1, kwargs...)
    if length(dims) > 1
        error("Only one dimension can be specified")
    end
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
    if length(dims) > 1
        error("Only one dimension can be specified")
    end
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
bootstrapmean(args...; kwargs...) = bootstrapaverage(mean, args...; kwargs...)

"""
hierarchicalkendall(x::AbstractVector{<:Real}, y::AbstractDimArray,
                         mode = :group; kwargs...)
Compute the Kendall's tau between a vector of hierarchy scores `x` and a 3D array of
features `y`. `y` should have dimensions `Depth`, `SessionID`, and `Structure`.
mode can be either `:group` or `:individual`. In the `:group` mode, the Kendall's tau
"""
function hierarchicalkendall(x::AbstractVector{<:Real}, y::AbstractDimArray,
                             mode = Val(:group); kwargs...)
    hasdim(y, Depth) ||
        throw(ArgumentError("Argument 2 should have a Depth dimension"))
    hasdim(y, SessionID) ||
        throw(ArgumentError("Argument 2 should have a SessionID dimension"))
    hasdim(y, Structure) ||
        throw(ArgumentError("Argument 2 should have a Structure dimension"))
    hierarchicalkendall(x, y, mode; kwargs...)
end
function hierarchicalkendall(x::AbstractVector{<:Real}, y::AbstractDimArray,
                             mode::Symbol; kwargs...)
    hierarchicalkendall(x, y, Val(mode); kwargs...)
end
function pairedkendall(xy)
    x = first.(xy)
    y = last.(xy)
    notnan = .!(isnan.(y) .| isnan.(x))
    corkendall(x[notnan], y[notnan])
end
function _hierarchicalkendall(xx, yy; N = 10000, confint = 0.95)
    b = Bootstrap.bootstrap(pairedkendall, collect(zip(xx, yy)),
                            Bootstrap.BalancedSampling(N))
    μ, σ... = only(Bootstrap.confint(b, Bootstrap.BCaConfInt(confint)))
    nsesh = hasdim(yy, SessionID) ? size(yy, SessionID) : 1
    μsur = map(1:N) do _
        idxs = stack(randperm(size(yy, Structure)) for _ in 1:nsesh)
        idxs = yy isa AbstractMatrix ? idxs' : idxs[:]
        ys = view(yy, idxs) # A different shuffle for each session
        @assert size(xx) == size(ys)
        pairedkendall(collect(zip(xx, ys)))
    end
    𝑝 = mean(abs.(μ) .< abs.(μsur))
    return μ, σ, 𝑝
end

function mediankendallpvalue(x::AbstractVector, Y::AbstractMatrix; N = 10000) # Take the correlation of each column
    τ = map(eachslice(collect(Y), dims = 2)) do y
        notnan = .!isnan.(y)
        corkendall(x[notnan], y[notnan])
    end
    τsur = asyncmap(1:N) do _
        idxs = randperm(length(x))
        taus = map(eachslice(collect(Y), dims = 2)) do y
            notnan = .!isnan.(y)
            corspearman(x[idxs][notnan], y[notnan])
        end
    end
    τsur = vcat(τsur...)
    𝑝 = MannWhitneyUTest(τ[:], τsur[:]) |> pvalue
    return median(τ), 𝑝
end

function hierarchicalkendall(x::AbstractVector{<:Real}, y::AbstractDimArray,
                             ::Val{:group}; kwargs...)
    y = permutedims(y, (Depth, SessionID, Structure)) # Check right orientation
    xx = repeat(x', size(y, SessionID), 1) # Sessionid × Structure
    ms = asyncmap(eachslice(y, dims = Depth)) do yy
        μ, σ, 𝑝 = _hierarchicalkendall(xx, yy; kwargs...)
    end
    𝑝 = last.(ms)
    𝑝 = set(𝑝, MultipleTesting.adjust(collect(𝑝), BenjaminiHochberg()))
    return first.(ms), getindex.(ms, 2), 𝑝
end
function hierarchicalkendall(x::AbstractVector{<:Real}, y::AbstractDimArray,
                             ::Val{:individual}; N)
    y = permutedims(y, (Depth, SessionID, Structure)) # Check right orientation
    ms = asyncmap(eachslice(y, dims = Depth)) do yy
        mnms = map(eachslice(yy, dims = SessionID)) do yyy # Individual 6-vector
            μ = corkendall(x, yyy)
            s = map(1:N) do _
                idxs = randperm(length(yyy))
                corkendall(x, yyy[idxs]) # A shuffle
            end
            return μ, s
        end
        μ = first.(mnms)
        s = last.(mnms)
        σ = (percentile(μ, 25), percentile(μ, 75))
        s = collect(Iterators.flatten(s))
        𝑝 = MannWhitneyUTest(s, μ) |> pvalue
        μ = median(μ)
        return μ, σ, 𝑝
    end
    𝑝 = last.(ms)
    𝑝 = set(𝑝, MultipleTesting.adjust(collect(𝑝), BenjaminiHochberg()))
    return first.(ms), getindex.(ms, 2), 𝑝
end
