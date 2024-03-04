module SpatiotemporalMotifs
import SpatiotemporalMotifs as SM
using DimensionalData
using PythonCall
using Conda
using StatsBase
using JLD2
using CairoMakie
using Foresight
set_theme!(foresight(:physics))
using Statistics
export plotdir

function _preamble()
    quote
        file = @__FILE__
        display(file)
        using Statistics
        using StatsBase
        using MultivariateStats
        using Random
        using LinearAlgebra
        import AllenNeuropixels as AN
        using AperiodicSurrogates
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
    end
end
macro preamble()
    _preamble()
end

plotdir(args...) = projectdir("plots", args...)
export plotdir

const structures = ["VISp", "VISl", "VISrl", "VISal", "VISpm", "VISam"]
const connector = "-"
const layers = ["1", "2/3", "4", "5", "6"]
const layercolors = cgrad(:inferno)[range(start = 0, stop = 1, length = length(layers) + 1)][1:length(layers)]
colors = [Foresight.cornflowerblue,
          Foresight.crimson,
          Foresight.cucumber,
          Foresight.california,
          Foresight.juliapurple]
const structurecolors = [Makie.RGB(0.5, 0.5, 0.5),
                         colors...]

if !haskey(ENV, "DRWATSON_STOREPATCH")
    ENV["DRWATSON_STOREPATCH"] = "true"
end

using IntervalSets
using Normalization
using Distributed
using Foresight
using TimeseriesTools.Unitful
using Peaks
using DimensionalData
using JLD2
using GeneralizedPhase
using TimeseriesTools
using ModulationIndices
using DataFrames
using ComplexityMeasures
using JSON
using DrWatson
import AllenNeuropixels as AN
import AllenNeuropixelsBase as ANB

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
            catch
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
    ddims = Tuple([Dim{Symbol(d)}(s) for (s, d) in zip(uvs, dims)])
    Q = Array{Bool}(undef, length.(uvs)...)
    Q = DimArray(Q, ddims)
    Q .= false
    for v in eachcol(vs)
        ds = [Dim{Symbol(d)}(At(s)) for (s, d) in zip(v, dims)]
        Q[ds...] = true
    end
    return Q
end

function commondepths(depths)
    # Find a Range of depths that best approximates the given collection of depths
    # * Get the smallest possible bound
    I = ceil(maximum(minimum.(depths)), sigdigits = 2) ..
        floor(minimum(maximum.(depths)), sigdigits = 2)
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

function send_powerspectra(sessionid; outpath = datadir("PowerSpectra"),
                           plotpath = plotdir("PowerSpectra", "full"),
                           rewrite = false)
    params = (;
              sessionid,
              epoch = :longest,
              pass = (1, 300))

    stimuli = ["spontaneous", "flash_250ms"] #, r"Natural_Images"]
    structures = SM.structures

    ssession = []

    for stimulus in stimuli
        for structure in structures
            @info "Loading $(stimulus) LFP in $(structure) for session $(params[:sessionid])"
            # if stimulus == "flashes" && partition
            #     @info "Partitioning $stimulus"
            #     _params = (; params..., stimulus, structure)
            #     LFP = AN.extracttheta(_params)u"V"
            # else
            _params = (; params..., stimulus, structure)
            if stimulus == r"Natural_Images"
                _params = (; _params..., epoch = (:longest, :active))
            end
            outfile = joinpath(outpath,
                               savename(Dict("sessionid" => params[:sessionid],
                                             "stimulus" => string(stimulus),
                                             "structure" => structure), "jld2",
                                        connector = "-"))
            @info outfile
            function isbad()
                if !isfile(outfile)
                    return true
                end
                f = jldopen(outfile)
                ind = haskey(f, "error")
                close(f)
                if ind
                    return true
                else
                    return false
                end
            end
            if !rewrite && !isbad()
                @info "Already calculated: $(stimulus), $(structure), $(params[:sessionid])"
                continue
            end
            if isempty(ssession) # Only initialize session if we have to
                session = AN.Session(params[:sessionid])
                push!(ssession, session)
            else
                session = ssession[1]
            end
            @info "Calculating $(stimulus) LFP in $(structure) for session $(params[:sessionid])"
            # try
            LFP = AN.formatlfp(session; tol = 4, _params...)u"V"
            # end
            LFP = set(LFP, Ti => Ti((times(LFP))u"s"))
            channels = lookup(LFP, :channel)
            # N = fit(UnitPower, LFP, dims=1)
            # normalize!(LFP, N)
            S = powerspectrum(LFP, 0.1; padding = 10000)

            S = S[Freq(params[:pass][1] * u"Hz" .. params[:pass][2] * u"Hz")]
            depths = AN.getchanneldepths(session, LFP; method = :probe)
            S = set(S, Dim{:channel} => Dim{:depth}(depths))

            # * Find peaks in the average power spectrum
            s = mean(S, dims = Dim{:depth})[:, 1]
            pks, vals = findmaxima(s, 10)
            pks, proms = peakproms(pks, s)
            promidxs = (proms ./ vals .> 0.25) |> collect
            maxs = maximum(S, dims = Dim{:depth})[:, 1]
            pks = pks[promidxs]
            pks = TimeseriesTools.freqs(s)[pks]
            vals = s[Freq(At(pks))]

            f = Figure()
            colorrange = extrema(depths)
            ax = Axis(f[1, 1]; xscale = log10, yscale = log10, xtickformat = "{:.1f}",
                      limits = (params[:pass], (nothing, nothing)), xgridvisible = true,
                      ygridvisible = true, topspinevisible = true,
                      title = "$(_params[:structure]), $(_params[:stimulus])",
                      xminorticksvisible = true, yminorticksvisible = true,
                      xminorgridvisible = true, yminorgridvisible = true,
                      xminorgridstyle = :dash)
            p = traces!(ax, S[2:end, :]; colormap = cgrad(sunset, alpha = 0.4),
                        linewidth = 3, colorrange)
            scatter!(ax, ustrip.(pks), collect(ustrip.(vals)), color = :black,
                     markersize = 10, marker = :dtriangle)
            text!(ax, ustrip.(pks), collect(ustrip.(vals));
                  text = string.(round.(eltype(pks), pks, digits = 1)),
                  align = (:center, :bottom), color = :black, rotation = 0,
                  fontsize = 12,
                  offset = (0, 3))

            c = Colorbar(f[1, 2]; label = "Channel depth (μm)", colorrange,
                         colormap = sunset)
            # rowsize!(f.layout, 1, Relative(0.8))
            mkpath(joinpath(plotpath, "$(_params[:sessionid])"))
            wsave(joinpath(plotpath, "$(_params[:sessionid])",
                           "$(_params[:stimulus])_$(_params[:structure]).pdf"), f)

            # * Format data for saving
            streamlinedepths = AN.getchanneldepths(session, LFP; method = :streamlines)
            layerinfo = AN.Plots._layerplot(session, channels)
            D = Dict(DimensionalData.metadata(LFP))
            @pack! D = channels, streamlinedepths, layerinfo
            S = rebuild(S; metadata = D)
            tagsave(outfile, @strdict S)
            # catch e
            #     @warn e
            #     tagsave(outfile, Dict("error" => sprint(showerror, e)))
            # end
            LFP = []
            S = []
            s = []
            f = []
            depths = []
            GC.gc()
        end
    end
end

function send_calculations(D;
                           pass_θ = (4, 10) .* u"Hz",
                           pass_γ = (40, 100) .* u"Hz",
                           ΔT = -0.10u"s" .. 0.6u"s",
                           doupsample = 0,
                           rewrite = false,
                           outpath)
    @unpack sessionid, structure, stimulus = D
    filename = joinpath(outpath, savename(D, "jld2"), connector = "-")
    if !rewrite && isfile(filename)
        @info "Calculations already complete for $(sessionid), $(structure), $(stimulus)"
        return GC.gc()
    end
    session = AN.Session(sessionid)

    probeid = AN.getprobe(session, structure)

    if stimulus == r"Natural_Images"
        LFP = AN.formatlfp(session; probeid, structure, stimulus, rectify = false,
                           epoch = (:longest, :active))
        trials = AN.getchangetrials(session)
        starttimes = trials.change_time_with_display_delay
    else
        LFP = AN.formatlfp(session; probeid, structure, stimulus, rectify = false,
                           epoch = :first)
        trials = AN.stimulusintervals(session, stimulus)
        starttimes = trials.start_time # For flashes
    end
    tmap = Timeseries(times(LFP), zeros(size(LFP, 1)))
    tmap[Ti(Near(starttimes))] .= starttimes
    tmap = rectify(tmap, dims = Ti, tol = 4)
    dev = tmap[findlast(tmap .> 1)] - times(tmap)[findlast(tmap .> 1)] # * Gonna have to fix this
    LFP = rectify(LFP, dims = Ti)
    @assert times(LFP) == times(tmap)
    starttimes = times(tmap[tmap .> 1]) |> collect # * The change times transformed to rectified time coordinates

    channels = lookup(LFP, :channel)
    depths = AN.getchanneldepths(session, LFP; method = :probe)
    streamlinedepths = AN.getchanneldepths(session, LFP; method = :streamlines)
    layerinfo = AN.Plots._layerplot(session, channels)
    LFP = set(LFP, Dim{:channel} => Dim{:depth}(depths))

    LFP = set(LFP, Ti => times(LFP) .* u"s")
    LFP = set(LFP, Dim{:depth} => lookup(LFP, :depth) .* u"μm")
    starttimes = starttimes .* u"s"
    tmap = set(tmap .* u"s", Ti => times(tmap) .* u"s")
    LFP = rectify(LFP, dims = Dim{:depth})

    # N = fit(RobustZScore, LFP, dims=1)
    # LFP = normalize(LFP, N)

    θ = bandpass(LFP, pass_θ)
    doupsample > 0 && (θ = upsample(θ, doupsample, Dim{:depth}))

    γ = bandpass(LFP, pass_γ)
    doupsample > 0 && (γ = upsample(γ, doupsample, Dim{:depth}))

    a = hilbert(θ)
    aᵧ = hilbert(γ)

    # r = abs.(aᵧ)
    # ∂r∂t = centralderiv(r, dims=Ti)
    # ∂r∂x = centralderiv(r, dims=Dim{:depth})
    # v̂ᵧ = ∂r∂x ./ ∂r∂t

    ϕ = angle.(a) # _generalized_phase(ustrip(θ)) #
    ϕᵧ = angle.(aᵧ) # _generalized_phase(ustrip(γ)) #

    ω = centralderiv(ϕ, dims = Ti, grad = phasegrad) # Angular Frequency
    k = -centralderiv(ϕ, dims = Dim{:depth}, grad = phasegrad) # Wavenumber
    # We have a minus sign because the hilbert transform uses the convention of ϕ = ωt (for univariate; phase increases with time), which implies ϕ = ωt - kx (for multivariate).
    v = ω ./ k # Phase velocity of theta

    ωᵧ = centralderiv(ϕᵧ, dims = Ti, grad = phasegrad) # Frequency
    kᵧ = -centralderiv(ϕᵧ, dims = Dim{:depth}, grad = phasegrad) # Wavenumber
    ∂ωᵧ = centralderiv(ωᵧ; dims = Ti)
    ∂kᵧ = centralderiv(kᵧ; dims = Ti)
    vₚ = ωᵧ ./ kᵧ # Phase velocity
    vᵧ = ∂ωᵧ ./ ∂kᵧ # Group velocity

    function alignmatchcat(x)
        x = align(x, starttimes, ΔT)[2:(end - 2)]
        x = matchdim(x; dims = 1)
        x = stack(Dim{:changetime}(times(x)), x)
    end
    x = alignmatchcat(θ)
    y = alignmatchcat(γ)
    a, ϕ, ϕᵧ, ω, k, v, aᵧ, ωᵧ, kᵧ, vₚ, vᵧ = alignmatchcat.([a, ϕ, ϕᵧ, ω, k, v, aᵧ, ωᵧ, kᵧ,
                                                            vₚ, vᵧ])
    out = Dict("channels" => channels,
               "trials" => trials[2:(end - 2), :],
               "streamlinedepths" => streamlinedepths,
               "layerinfo" => layerinfo,
               "pass_θ" => pass_θ,
               "pass_γ" => pass_γ,
               "x" => x,
               "y" => y,
               "a" => a,
               "ϕ" => ϕ,
               "ϕᵧ" => ϕᵧ,
               "ω" => ω,
               "k" => k,
               "v" => v,
               "aᵧ" => aᵧ,
               # "r" => r,
               # "∂r∂t" => ∂r∂t,
               # "∂r∂x" => ∂r∂x,
               # "v̂ᵧ" => v̂ᵧ,
               "ωᵧ" => ωᵧ,
               "kᵧ" => kᵧ,
               "vₚ" => vₚ,
               "vᵧ" => vᵧ)
    @tagsave filename out

    x = []
    y = []
    a = []
    ϕ = []
    ϕᵧ = []
    ω = []
    k = []
    v = []
    aᵧ = []
    r = []
    ∂r∂t = []
    ∂r∂x = []
    v̂ᵧ = []
    ωᵧ = []
    kᵧ = []
    vₚ = []
    vᵧ = []
    GC.gc()
end

function send_SaveBehaviorResponse(sessionid;
                                   structures = ["VISp", "VISl", "VISrl", "VISal", "VISpm",
                                                 "VISam"], kwargs...)
    session = AN.Session(sessionid)
    probestructures = AN.getprobestructures(session, structures)
    for (probeid, structure) in probestructures
        for stimulus in ["flash_250ms", r"Natural_Images"]
            @info "Saving $(sessionid) $(structure) $(stimulus)"
            D = Dict()
            @pack! D = sessionid, structure, stimulus
            send_SaveBehaviorResponse(sessionid, structure, stimulus; kwargs...)
        end
    end
end

end
