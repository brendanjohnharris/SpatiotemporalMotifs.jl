
function send_powerspectra(sessionid; outpath = datadir("power_spectra"),
                           plotpath = plotdir("power_spectra", "full"),
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
            # if stimulus == "flashes" && partition @info "Partitioning $stimulus" _params =
            #     (; params..., stimulus, structure) LFP = AN.extracttheta(_params)u"V" else
            _params = (; params..., stimulus, structure)
            if stimulus == r"Natural_Images"
                _params = (; _params..., epoch = (:longest, :active))
            end
            outfile = savepath(Dict("sessionid" => params[:sessionid],
                                    "stimulus" => string(stimulus),
                                    "structure" => structure), "jld2", outpath)
            @info outfile

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
            LFP = AN.formatlfp(session; tol = 3, _params...)u"V"
            # end
            LFP = set(LFP, Ti => Ti((times(LFP))u"s"))
            channels = lookup(LFP, :channel)
            # N = fit(UnitPower, LFP, dims=1) normalize!(LFP, N)
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
            # catch e @warn e tagsave(outfile, Dict("error" => sprint(showerror, e))) end
            LFP = []
            S = []
            s = []
            f = []
            depths = []
            GC.gc()
        end
    end
end

function _calculations(LFP; pass_θ, pass_γ, ΔT, doupsample, starttimes)
    θ = bandpass(LFP, pass_θ)
    doupsample > 0 && (θ = upsample(θ, doupsample, Dim{:depth}))

    γ = bandpass(LFP, pass_γ)
    doupsample > 0 && (γ = upsample(γ, doupsample, Dim{:depth}))

    a = hilbert(θ)
    aᵧ = hilbert(γ)

    # r = abs.(aᵧ)

    ϕ = angle.(a) # _generalized_phase(ustrip(θ)) #
    ϕᵧ = angle.(aᵧ) # _generalized_phase(ustrip(γ)) #

    ω = centralderiv(ϕ, dims = Ti, grad = phasegrad) # Angular Frequency

    k = -centralderiv(ϕ, dims = Dim{:depth}, grad = phasegrad) # Wavenumber
    # We have a minus sign because the hilbert transform uses the convention of ϕ = ωt (for
    # univariate; phase increases with time), which implies ϕ = ωt - kx (for multivariate).
    v = ω ./ k # Phase velocity of theta

    ωᵧ = centralderiv(ϕᵧ, dims = Ti, grad = phasegrad) # Frequency
    kᵧ = -centralderiv(ϕᵧ, dims = Dim{:depth}, grad = phasegrad) # Wavenumber
    # ∂ωᵧ = centralderiv(ωᵧ; dims = Ti)
    # ∂kᵧ = centralderiv(kᵧ; dims = Ti)
    # vₚ = ωᵧ ./ kᵧ # Phase velocity
    # vᵧ = ∂ωᵧ ./ ∂kᵧ # Group velocity

    function alignmatchcat(x)
        x = align(x, starttimes, ΔT)[2:(end - 2)]
        x = matchdim(x; dims = 1)
        x = stack(Dim{:changetime}(times(x)), x)
    end
    # a, ϕ, ϕᵧ, ω, k, v, aᵧ, ωᵧ, kᵧ, vₚ, vᵧ = alignmatchcat.([a, ϕ, ϕᵧ, ω, k, v, aᵧ, ωᵧ, kᵧ,
    #                                                         vₚ, vᵧ])
    V, x, y, a, ϕ, ϕᵧ, ω, k, v, aᵧ, ωᵧ, kᵧ = alignmatchcat.([LFP, θ, γ, a, ϕ, ϕᵧ, ω, k,
                                                             v,
                                                             aᵧ, ωᵧ,
                                                             kᵧ])
    r = abs.(aᵧ) .* unit(eltype(V))

    return V, x, y, a, ϕ, ϕᵧ, ω, k, v, r, ωᵧ, kᵧ
end

function send_calculations(D::Dict, session = AN.Session(D[:sessionid]);
                           pass_θ = (4, 10) .* u"Hz",
                           pass_γ = (40, 100) .* u"Hz",
                           ΔT = -0.5u"s" .. 0.75u"s",
                           doupsample = 0,
                           rewrite = false,
                           outpath)
    @unpack sessionid, structure, stimulus = D
    outfile = savepath(D, "jld2", outpath)

    @info outfile
    if !rewrite && isfile(outfile)
        @info "Calculations already complete for $(sessionid), $(structure), $(stimulus)"
        return GC.gc()
    end

    probeid = AN.getprobe(session, structure)
    performance_metrics = AN.getperformancemetrics(session)

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
    _tmap = Timeseries(times(LFP), zeros(size(LFP, 1)))
    tmap = deepcopy(_tmap)
    tmap[Ti(Near(starttimes))] .= starttimes
    tmap = rectify(tmap, dims = Ti, tol = 3)
    dev = tmap[findlast(tmap .> 1)] - times(tmap)[findlast(tmap .> 1)] # * Gonna have to fix this
    LFP = rectify(LFP, dims = Ti, tol = 3)
    @assert times(LFP) == times(tmap)
    starttimes = times(tmap[tmap .> 1]) |> collect # * The change times transformed to rectified time coordinates
    spiketimes = AN.getspiketimes(session, structure)
    units = AN.getunitmetrics(session)
    units = units[units.ecephys_unit_id .∈ [keys(spiketimes)], :]
    spiketimes = map(collect(spiketimes)) do (u, ts) # * Rectify the spike times
        tmap = deepcopy(_tmap)
        tmap[Ti(Near(ts))] .= ts
        tmap = rectify(tmap, dims = Ti, tol = 3)
        ts = times(tmap[tmap .> 1]) |> collect
        return u => ts
    end |> Dict

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

    V, x, y, a, ϕ, ϕᵧ, ω, k, v, r, ωᵧ, kᵧ = _calculations(LFP; pass_θ, pass_γ, ΔT,
                                                          doupsample, starttimes)

    out = Dict("channels" => channels,
               "trials" => trials[2:(end - 2), :],
               "streamlinedepths" => streamlinedepths,
               "layerinfo" => layerinfo,
               "pass_θ" => pass_θ,
               "pass_γ" => pass_γ,
               "V" => V,
               "x" => x,
               "y" => y,
               "a" => a,
               "ϕ" => ϕ,
               #    "ϕᵧ" => ϕᵧ,
               "ω" => ω,
               "k" => k,
               "v" => v,
               "r" => r,
               "ωᵧ" => ωᵧ,
               "kᵧ" => kᵧ,
               #    "vₚ" => vₚ,
               #    "vᵧ" => vᵧ,
               "units" => units,
               "spiketimes" => spiketimes,
               "performance_metrics" => performance_metrics)
    @tagsave outfile out

    out = V = x = y = a = ϕ = ϕᵧ = ω = k = v = r = ωᵧ = kᵧ = spiketimes = []
    GC.gc()
    return true
end

function send_calculations(sessionid;
                           structures = ["VISp", "VISl", "VISrl", "VISal", "VISpm",
                                         "VISam"], kwargs...)
    session = AN.Session(sessionid)
    probestructures = AN.getprobestructures(session, structures)
    for (probeid, structure) in probestructures
        for stimulus in ["flash_250ms", r"Natural_Images"]
            @info "Saving $(sessionid) $(structure) $(stimulus)"
            D = @dict sessionid structure stimulus
            send_calculations(D, session; kwargs...)
        end
    end
    return true
end

function test_calculations(args...; kwargs...)
    while true
        X = randn(20000, 20000) # Requires at least 3Gb memory
        sleep(120)
        @info "Poll: $(ENV["PBS_JOBID"])"
    end
end

function load_performance(; path = datadir("calculations"), stimulus = r"Natural_Images")
    Q = calcquality(path)[structure = At(structures)]
    out = map([lookup(Q, :structure) |> first]) do structure
        out = map(lookup(Q, :sessionid)) do sessionid
            if Q[sessionid = At(sessionid), structure = At(structure)] == 0
                return nothing
            end
            filename = savepath((@strdict sessionid structure stimulus), "jld2", path)
            f = jldopen(filename, "r")
            @unpack performance_metrics = f
            close(f)
            return performance_metrics
        end
        out = filter(!isnothing, out)
    end
    performance = only(out)

    newsessions = performance |> DataFrame
    newsessions.sessionid = lookup(Q, :sessionid) .|> Int
    return newsessions
end

function load_calculations(Q; path = datadir("calculations"), stimulus, vars = [:x, :k])
    out = map(lookup(Q, :structure)) do structure
        out = map(lookup(Q, :sessionid)) do sessionid
            if Q[sessionid = At(sessionid), structure = At(structure)] == 0
                return nothing
            end
            filename = savepath((@strdict sessionid structure stimulus), "jld2", path)
            f = jldopen(filename, "r")
            @unpack streamlinedepths, layerinfo, pass_γ, pass_θ, performance_metrics, spiketimes, trials = f
            ovars = Dict()
            for v in vars
                push!(ovars, v => f[string(v)])
            end
            close(f)

            # Depth and layer info
            layernames = DimArray(layerinfo[1],
                                  Dim{:depth}(lookup(ovars[first(vars)], Dim{:depth})))
            layernums = DimArray(layerinfo[3],
                                 Dim{:depth}(lookup(ovars[first(vars)], Dim{:depth})))

            # Remove poor quality depth estimates
            idxs = 1:size(ovars[first(vars)], 2)
            try
                while !issorted(streamlinedepths)
                    idxs = idxs[1:(end - 1)]
                    streamlinedepths = streamlinedepths[idxs]
                end
            catch e
                return nothing
            end
            layernames = layernames[idxs]

            ovars = map(collect(ovars)) do (v, k)
                k = k[:, idxs, :]

                # We know the depths are sorted from above
                k = set(k, Dim{:depth}(streamlinedepths))
                k = set(k,
                        Dim{:depth} => DimensionalData.Irregular(extrema(streamlinedepths)))

                @assert issorted(lookup(k, Dim{:depth}))
                push!(k.metadata.val, :layernames => layernames)
                return v => k
            end
            return (; streamlinedepths, layernames, pass_γ, pass_θ, trials, sessionid,
                    performance_metrics, spiketimes, ovars...)
        end
        out = filter(!isnothing, out)
        out = filter(x -> maximum(x[:streamlinedepths]) > 0.90, out) # Remove sessions that don't have data to a reasonable depth
    end
    if length(unique(length.(out))) != 1
        @warn "Not all structures returned the same sessions; filtering to common sessions"
        commonsessions = [[o[:sessionid] for o in O]
                          for O in out]
        commonsessions = intersect(commonsessions...)
        out = map(out) do O
            filter(o -> (o[:sessionid] in commonsessions), O)
        end
    end
    return out
end

function unify_calculations(out; vars = [:x, :k])
    stimuli = [[DimensionalData.metadata(o[first(vars)])[:stimulus] for o in out[i]]
               for i in eachindex(out)]
    stimulus = only(unique(vcat(stimuli...)))
    uni = map(out) do o
        sessionids = getindex.(o, :sessionid)
        streamlinedepths = getindex.(o, :streamlinedepths)
        unidepths = commondepths(streamlinedepths)
        layernames = map(o) do p
            l = p[:layernames]
            l = set(l, Dim{:depth} => p[:streamlinedepths])
            if last(l) ∈ ["or", "scwm"] # Sometime anomalies at the boundary; we fall back on the channel structure labels, the ground truth, for confidence that this is still a cortical channel
                l[end] = l[end - 1]
            end
            l[depth = Near(unidepths)]
        end
        layernums = [DimArray(parselayernum.(parent(p)), dims(p)) for p in layernames]

        # * We want to get a unified array, plus (potentially overlapping) intervals for
        # * each layer
        layerints = map(unique(vcat(parent.(layernums)...))) do l
            depths = map(layernums) do ls
                if sum(ls .== l) == 0 # !!! BROKEN!!!
                    m = mean([maximum(lookup(ls[ls .== l - 1], :depth)),
                              minimum(lookup(ls[ls .== l + 1], :depth))])
                    [m, m]
                else
                    collect(extrema(lookup(ls, :depth)[ls .== l]))
                end
            end
            Interval(extrema(vcat(depths...))...)
        end
        ovars = Dict()
        for v in vars
            k = map(o) do p
                k = p[v]
                k = k[depth = Near(unidepths)]
                k = set(k, Dim{:depth} => Dim{:depth}(unidepths))
                if stimulus == r"Natural_Images"
                    d = Dim{:trial}(p[:trials].hit[1:size(k,
                                                          :changetime)])
                else
                    d = Dim{:trial}(1:size(k, :changetime))
                end
                k = set(k, Dim{:changetime} => Dim{:trial})
                set(k, DimensionalData.format(d, lookup(k, :trial)))
            end
            mints = lookup.(k, Ti)
            _, minti = findmin(length.(mints))
            mints = mints[minti]
            k = map(k) do p # Match to smallest time interval. Should only differ by a sample or so
                p[Ti(At(mints))]
            end
            if stimulus == r"Natural_Images"
                k = cat(k...; dims = Dim{:trial}(vcat(lookup.(k, :trial)...)))
            else
                k = stack(Dim{:sessionid}(sessionids), k)
            end
            push!(ovars, v => k)
        end

        return (; unidepths, layerints, layernames, layernums, ovars...)
    end
    return uni
end
