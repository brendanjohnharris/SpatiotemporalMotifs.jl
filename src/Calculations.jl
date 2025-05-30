using DataFrames
using JLD2, CodecZstd
using TimeseriesTools

function powerspectra_quality(sessionid, stimulus, structure;
                              outpath = datadir("power_spectra"),
                              plotpath = datadir("plots", "power_spectra_plots"),
                              rewrite = false, retry_errors = true)
    plotfile = joinpath(plotpath, "$(sessionid)",
                        "$(stimulus)_$(structure)_pac.pdf")
    outfile = savepath(Dict("sessionid" => sessionid,
                            "stimulus" => stimulus,
                            "structure" => structure), "jld2", outpath)
    badis = rewrite || isbad(outfile; retry_errors, check_other_file = plotfile)
    return !badis
end

function send_powerspectra(sessionid, stimulus, structure;
                           outpath = datadir("power_spectra"),
                           plotpath = datadir("plots", "power_spectra_plots"),
                           rewrite = false, retry_errors = true)
    params = (;
              sessionid,
              epoch = :longest,
              pass = (1, 300)) # To match typical range, e.g., https://www.biorxiv.org/content/10.1101/2021.07.28.454235v2.full

    ssession = []
    GC.safepoint()
    @info "Loading $(stimulus) LFP in $(structure) for session $(params[:sessionid])"
    # if stimulus == "flashes" && partition @info "Partitioning $stimulus" _params =
    #     (; params..., stimulus, structure) LFP = AN.extracttheta(_params)u"V" else
    _params = (; params..., stimulus, structure)
    if stimulus == r"Natural_Images"
        _params = (; _params..., epoch = (:longest, :active))
    end
    outfile = savepath(Dict("sessionid" => params[:sessionid],
                            "stimulus" => stimulus,
                            "structure" => structure), "jld2", outpath)
    @info outfile

    plotfile = joinpath(plotpath, "$(_params[:sessionid])",
                        "$(_params[:stimulus])_$(_params[:structure]).pdf")

    if !rewrite && !isbad(outfile; retry_errors, check_other_file = plotfile)
        @info "Already calculated: $(stimulus), $(structure), $(params[:sessionid])"
        return
    end
    if isempty(ssession) # Only initialize session if we have to
        session = AN.Session(params[:sessionid])
        push!(ssession, session)
    else
        session = ssession[1]
    end
    @info "Calculating $(stimulus) LFP in $(structure) for session $(params[:sessionid])"
    try
        probestructures = AN.getprobestructures(session)
        probestructures = unique(vcat(values(probestructures)...))
        if structure ∉ probestructures
            str = "Region error: structure $(structure) not found in $(params[:sessionid])"
            @warn str
            tagsave(outfile, Dict("error" => str))
            return
        end

        LFP = AN.formatlfp(session; tol = 3, _params...)u"V"

        LFP = set(LFP, 𝑡 => 𝑡((times(LFP))u"s"))
        channels = lookup(LFP, Chan)
        # N = fit(UnitPower, LFP, dims=1) normalize!(LFP, N)
        S = powerspectrum(LFP, 0.1; padding = 10000)
        S = S[Freq(params[:pass][1] * u"Hz" .. params[:pass][2] * u"Hz")]
        depths = AN.getchanneldepths(session, LFP; method = :probe)
        S = set(S, Chan => Depth(depths))

        # * Find peaks in the average power spectrum
        s = mean(S, dims = Depth)[:, 1]
        pks, vals = findmaxima(s, 10)
        pks, proms = peakproms(pks, s)
        promidxs = (proms ./ vals .> 0.25) |> collect
        maxs = maximum(S, dims = Depth)[:, 1]
        pks = pks[promidxs]
        pks = TimeseriesTools.freqs(s)[pks]
        vals = s[Freq(At(pks))]

        begin
            f = Figure()
            colorrange = extrema(depths)
            ax = Axis(f[1, 1]; xscale = log10, yscale = log10,
                      xtickformat = "{:.1f}",
                      limits = (params[:pass], (nothing, nothing)),
                      xgridvisible = true,
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
            wsave(plotfile, f)
            @info "Saved plot to $plotfile"
        end

        # * Calculate a comodulogram for these time series
        c = x -> comodulogram(ustripall(x); fs = ustripall(samplingrate(LFP)),
                              fₚ = 3:0.25:12,
                              fₐ = 25:1:125, dp = 2, da = 20)
        N = min(360000, size(LFP, 1))
        C = progressmap(c, eachslice(LFP[1:N, :], dims = 2); parallel = true)
        C = stack(dims(LFP[1:N, :], Chan), ToolsArray.(C), dims = 3)
        sLFP = mapslices(x -> surrogate(ustripall(parent(x)), IAAFT()), LFP[1:N, :],
                         dims = 1)
        sC = progressmap(c, eachslice(sLFP, dims = 2); parallel = true)
        sC = stack(dims(LFP[1:N, :], Chan), ToolsArray.(sC), dims = 3)

        begin
            f = Figure()
            ax = Axis(f[1, 1]; xlabel = "Phase frequency (Hz)",
                      ylabel = "Phase frequency (HZ)")
            mc = dropdims(mean(C, dims = 3); dims = 3)
            colorrange = (0, maximum(mc))
            p = heatmap!(ax, lookup(C, 1), lookup(C, 2), collect(ustripall(mc));
                         colorrange,
                         colormap = seethrough(reverse(sunrise)), rasterize = 5)
            c = Colorbar(f[1, 2], p; label = "Mean modulation strength")
            mkpath(joinpath(plotpath, "$(_params[:sessionid])"))
            plotfile = joinpath(plotpath, "$(_params[:sessionid])",
                                "$(_params[:stimulus])_$(_params[:structure])_pac.pdf")
            wsave(plotfile, f)
            @info "Saved plot to $plotfile"
        end

        # * Format data for saving
        streamlinedepths = AN.getchanneldepths(session, LFP; method = :streamlines)
        layerinfo = AN.Plots._layerplot(session, channels)

        begin # * Spontaneous order parameter. Need to rectify for this
            LFP = set(LFP, Chan => Depth(depths))
            origts = deepcopy(times(LFP))
            LFP = rectify(LFP; dims = Depth)
            θ = bandpass(LFP, SpatiotemporalMotifs.THETA .* u"Hz")
            ϕ = analyticphase(θ)
            ω = centralderiv(ϕ, dims = 𝑡, grad = phasegrad)

            k = -centralderiv(ϕ, dims = Depth, grad = phasegrad)

            k[ω .< 0u"Hz"] .= NaN * unit(eltype(k))
            R = dropdims(nansafe(mean, dims = Depth)(sign.(k)); dims = Depth)

            # * Surrogates
            idxs = randperm(size(ϕ, Depth))
            ϕs = set(ϕ, ϕ[:, idxs]) # Spatially shuffle channels
            ωs = centralderiv(ϕs, dims = 𝑡, grad = phasegrad)
            ks = -centralderiv(ϕs, dims = Depth, grad = phasegrad)
            ks[ωs .< 0u"Hz"] .= NaN * unit(eltype(ks))
            sR = dropdims(nansafe(mean, dims = Depth)(sign.(k)); dims = Depth)
        end

        begin # * Spontaneous spike-LFP coupling. We set the LFP time indices back to their original, non-rectified values
            unitdepths = produce_unitdepths(session)
            unitdepths[!, :stimulus] .= stimulus

            spiketimes = AN.getspiketimes(session, structure)
            units = AN.getunitmetrics(session)
            units = units[units.ecephys_unit_id .∈ [keys(spiketimes)], :]
            unitdepths = unitdepths[unitdepths.ecephys_unit_id .∈ [units.id], :]

            γ = bandpass(LFP, SpatiotemporalMotifs.GAMMA .* u"Hz")
            r = abs.(hilbert(γ))
            r[ω .< 0u"Hz"] .= NaN * unit(eltype(r))
            ϕ[ω .< 0u"Hz"] .= NaN * unit(eltype(ϕ))

            LFP = set(LFP, 𝑡 => origts)
            ϕ = set(ϕ, 𝑡 => origts)
            r = set(r, 𝑡 => origts)
            @assert dims(r, 1) isa DimensionalData.TimeDim
            r = ustripall(r) # Remove units, not used from here on
            N = fit(nansafe(HalfZScore), r; dims = 1) # * Crucial; normalization makes this correlation-like. Dimension 1 is time (check)
            normalize!(r, N)

            unitdepths[:, :spc] .= NaN
            unitdepths[:, :spc_angle] .= NaN
            unitdepths[:, :spc_pvalue] .= NaN
            unitdepths[:, :sac] .= NaN

            for u in 1:size(unitdepths, 1)
                unitid = unitdepths[u, :].ecephys_unit_id
                depth = unitdepths[u, :].probedepth
                _ϕ = ϕ[Depth(Near(depth))]
                _r = r[Depth(Near(depth))]
                spikes = spiketimes[unitid]

                pγ, p, 𝑝 = ppc(ustripall(_ϕ), spikes)
                rγ = sac(_r, spikes)
                unitdepths[u, :].spc = pγ
                unitdepths[u, :].spc_angle = p
                unitdepths[u, :].spc_pvalue = 𝑝
                unitdepths[u, :].sac = rγ
            end
        end

        D = Dict(DimensionalData.metadata(LFP))
        @pack! D = channels, streamlinedepths, layerinfo
        S = rebuild(S; metadata = D)
        outD = Dict("S" => S .|> Float32,
                    "C" => C .|> Float32,
                    "sC" => sC .|> Float32,
                    "R" => R .|> Float32,
                    "sR" => sR .|> Float32,
                    "unitdepths" => unitdepths)
        tagsave(outfile, outD)

        GC.safepoint()
        GC.gc()
    catch e
        GC.safepoint()
        GC.gc()
        @warn e
        tagsave(outfile, Dict("error" => sprint(showerror, e)))
    end
    @info "Finished calculations"
    GC.safepoint()
    GC.gc()
end

function _calculations(session::AN.AbstractSession, structure, stimulus)
    probeid = AN.getprobe(session, structure)

    if (stimulus == r"Natural_Images") || (stimulus == "Natural_Images_omission")
        _stimulus = r"Natural_Images"
        LFP = AN.formatlfp(session; probeid, structure, stimulus = _stimulus,
                           rectify = false, epoch = (:longest, :active))
        if stimulus == r"Natural_Images"
            trials = AN.getchangetrials(session)
            starttimes = trials.change_time_with_display_delay
        elseif stimulus == "Natural_Images_omission"
            stimuli = AN.getstimuli(session)
            stimuli = stimuli[.!ismissing.(stimuli.omitted), :]
            stimuli = stimuli[identity.(stimuli.omitted), :] # Just the omission trials
            stimuli = stimuli[stimuli.active, :] # * just the active-epoch omissions
            starttimes = stimuli.start_time
            trials = stimuli
        end
    elseif stimulus == "Natural_Images_passive"
        _stimulus = r"Natural_Images"
        LFP = AN.formatlfp(session; probeid, structure, stimulus = _stimulus,
                           rectify = false,
                           epoch = (:longest, :active => ByRow(==(false))))
        stimuli = AN.getstimuli(session)
        stimuli = stimuli[stimuli.start_time .∈ [Interval(LFP)], :]
        stimuli = stimuli[.!ismissing.(stimuli.is_change), :]
        stimuli = stimuli[identity.(stimuli.is_change), :]
        stimuli = stimuli[.!stimuli.active, :] # * just the passive changes
        starttimes = stimuli.start_time
        trials = stimuli
    else
        LFP = AN.formatlfp(session; probeid, structure, stimulus, rectify = false,
                           epoch = :first)
        trials = AN.stimulusintervals(session, stimulus)
        starttimes = trials.start_time # For flashes
    end
    _tmap = Timeseries(times(LFP), zeros(size(LFP, 1)))
    tmap = deepcopy(_tmap)
    tmap[𝑡(Near(starttimes))] .= starttimes
    tmap = rectify(tmap, dims = 𝑡, tol = 3)
    dev = tmap[findlast(tmap .> 1)] - times(tmap)[findlast(tmap .> 1)]
    LFP = rectify(LFP, dims = 𝑡, tol = 3)
    @assert times(LFP) == times(tmap)
    starttimes = times(tmap[tmap .> 1]) |> collect # * The change times transformed to rectified time coordinates
    spiketimes = AN.getspiketimes(session, structure)
    units = AN.getunitmetrics(session)
    units = units[units.ecephys_unit_id .∈ [keys(spiketimes)], :]
    spiketimes = map(collect(spiketimes)) do (u, ts) # * Rectify the spike times
        tmap = deepcopy(_tmap)
        tmap[𝑡(Near(ts))] .= ts
        tmap = rectify(tmap, dims = 𝑡, tol = 3)
        ts = times(tmap[tmap .> 1]) |> collect
        return u => Float32.(ts)
    end |> Dict

    channels = lookup(LFP, Chan)
    depths = AN.getchanneldepths(session, LFP; method = :probe)
    streamlinedepths = AN.getchanneldepths(session, LFP; method = :streamlines)
    LFP = set(LFP, Chan => Depth(depths))

    LFP = set(LFP, 𝑡 => times(LFP) .* u"s")
    LFP = set(LFP, Depth => lookup(LFP, Depth) .* u"μm")
    starttimes = starttimes .* u"s"
    tmap = set(tmap .* u"s", 𝑡 => times(tmap) .* u"s")
    LFP = rectify(LFP, dims = Depth)

    output = @strdict LFP trials starttimes channels depths streamlinedepths spiketimes units

    return output
end

function aligned_calculations(LFP; pass_θ, pass_γ, ΔT, doupsample, starttimes)
    θ = bandpass(LFP, pass_θ)
    doupsample > 0 && (θ = upsample(θ, doupsample, Depth))

    γ = bandpass(LFP, pass_γ)
    doupsample > 0 && (γ = upsample(γ, doupsample, Depth))

    a = hilbert(θ)
    aᵧ = hilbert(γ)

    # r = abs.(aᵧ)

    ϕ = angle.(a) # _generalized_phase(ustrip(θ)) #
    ϕᵧ = angle.(aᵧ) # _generalized_phase(ustrip(γ)) #

    ω = centralderiv(ϕ, dims = 𝑡, grad = phasegrad) # Angular Frequency

    k = -centralderiv(ϕ, dims = Depth, grad = phasegrad) # Wavenumber
    # We have a minus sign because the hilbert transform uses the convention of ϕ = ωt (for
    # univariate; phase increases with time for positive frequencies), which implies ϕ = ωt - kx (for multivariate).
    v = ω ./ k # Phase velocity of theta

    ωᵧ = centralderiv(ϕᵧ, dims = 𝑡, grad = phasegrad) # Frequency
    kᵧ = -centralderiv(ϕᵧ, dims = Depth, grad = phasegrad) # Wavenumber
    # ∂ωᵧ = centralderiv(ωᵧ; dims = 𝑡)
    # ∂kᵧ = centralderiv(kᵧ; dims = 𝑡)
    # vₚ = ωᵧ ./ kᵧ # Phase velocity
    # vᵧ = ∂ωᵧ ./ ∂kᵧ # Group velocity

    csd = centralderiv(centralderiv(LFP, dims = Depth); dims = Depth)

    function alignmatchcat(x)
        x = align(x, starttimes, ΔT)[2:(end - 2)]
        x = matchdim(x; dims = 1)
        x = stack(Dim{:changetime}(times(x)), x)
    end

    output = @strdict LFP csd θ γ a ϕ ϕᵧ ω k v aᵧ ωᵧ kᵧ
    for k in keys(output)
        output[k] = alignmatchcat(output[k])
    end
    output["r"] = abs.(output["aᵧ"]) .* unit(eltype(output["V"]))

    for k in keys(output)
        if k == "a"
            output[k] = map(ComplexF32, output[k])
        else
            output[k] = map(Float32, output[k])
        end
    end

    return output
end

function send_calculations(D::Dict, session = AN.Session(D[:sessionid]);
                           pass_θ = THETA .* u"Hz",
                           pass_γ = GAMMA .* u"Hz",
                           ΔT = -0.5u"s" .. 0.75u"s",
                           doupsample = 0,
                           rewrite = false,
                           outpath)
    @unpack sessionid, structure, stimulus = D
    outfile = savepath(D, "jld2", outpath)

    @info outfile
    if !rewrite && isfile(outfile)
        fl = jldopen(outfile, "r")
        if !haskey(fl, "error")
            @info "Calculations already complete for $(sessionid), $(structure), $(stimulus)"
            close(fl)
            return GC.gc()
        else
            close(fl)
        end
    end

    performance_metrics = AN.getperformancemetrics(session)

    outputs = _calculations(session, structure, stimulus)

    layerinfo = AN.Plots._layerplot(session, outputs["channels"])

    aligned_outputs = aligned_calculations(outputs["LFP"]; pass_θ, pass_γ, ΔT,
                                           doupsample,
                                           starttimes = outputs["starttimes"])

    @unpack aligned_outputs

    out = Dict("channels" => channels,
               "trials" => trials[2:(end - 2), :],
               "streamlinedepths" => streamlinedepths,
               "layerinfo" => layerinfo,
               "pass_θ" => pass_θ,
               "pass_γ" => pass_γ,
               "V" => V,
               "csd" => csd .|> Float32,
               "x" => x .|> Float32,
               "y" => y .|> Float32,
               "a" => a .|> ComplexF32,
               "ϕ" => ϕ .|> Float32,
               #    "ϕᵧ" => ϕᵧ,
               "ω" => ω .|> Float32,
               "k" => k .|> Float32,
               "v" => v .|> Float32,
               "r" => r .|> Float32,
               "ωᵧ" => ωᵧ .|> Float32,
               "kᵧ" => kᵧ .|> Float32,
               #    "vₚ" => vₚ,
               #    "vᵧ" => vᵧ,
               "units" => units,
               "spiketimes" => spiketimes,
               "performance_metrics" => performance_metrics)
    @tagsave outfile out compress=ZstdFrameCompressor()

    out = V = csd = x = y = a = ϕ = ϕᵧ = ω = k = v = r = ωᵧ = kᵧ = spiketimes = units = []
    GC.gc()
    return filenames
end

function send_calculations(sessionid;
                           structures = ["VISp", "VISl", "VISrl",
                               "VISal", "VISpm", "VISam"],
                           stimuli = ["flash_250ms",
                               r"Natural_Images",
                               "Natural_Images_passive"],
                           outpath = datadir("calculations"),
                           kwargs...)
    session = AN.Session(sessionid)
    probestructures = AN.getprobestructures(session, structures)
    for (probeid, structure) in probestructures
        for stimulus in stimuli
            @info "Saving $(sessionid) $(structure) $(stimulus)"
            D = @dict sessionid structure stimulus
            try
                send_calculations(D, session; outpath, kwargs...)
            catch e
                outfile = savepath(D, "jld2", outpath)
                @warn "Error in $(outfile). Deleting file."
                rm(outfile)
                throw(e)
            end
        end
    end
    return true
end

function send_thalamus_calculations(D::Dict, session = AN.Session(D[:sessionid]);
                                    pass_θ = THETA .* u"Hz",
                                    pass_γ = GAMMA .* u"Hz",
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

    performance_metrics = AN.getperformancemetrics(session)
    output = _calculations(session, structure, stimulus)
    @unpack output

    θ = bandpass(LFP, pass_θ)
    doupsample > 0 && (θ = upsample(θ, doupsample, Depth))

    γ = bandpass(LFP, pass_γ)
    doupsample > 0 && (γ = upsample(γ, doupsample, Depth))

    a = hilbert(θ)
    aᵧ = hilbert(γ)

    ϕ = angle.(a)
    ϕᵧ = angle.(aᵧ)

    csd = centralderiv(centralderiv(LFP, dims = Depth); dims = Depth)

    function alignmatchcat(x)
        x = align(x, starttimes, ΔT)[2:(end - 2)]
        x = matchdim(x; dims = 1)
        x = stack(Dim{:changetime}(times(x)), x)
    end
    V, csd, x, y, a, ϕ, ϕᵧ, aᵧ = alignmatchcat.([LFP, csd, θ, γ, a, ϕ, ϕᵧ, aᵧ])
    r = abs.(aᵧ) .* unit(eltype(V))

    out = Dict("channels" => channels,
               "trials" => trials[2:(end - 2), :],
               "pass_θ" => pass_θ,
               "pass_γ" => pass_γ,
               "V" => V .|> Float32,
               "csd" => csd,
               "x" => x .|> Float32,
               "y" => y .|> Float32,
               "a" => a .|> ComplexF32,
               "ϕ" => ϕ .|> Float32,
               "r" => r .|> Float32,
               "units" => units,
               "spiketimes" => spiketimes,
               "performance_metrics" => performance_metrics)
    @tagsave outfile out compress=ZstdFrameCompressor()

    out = V = csd = x = y = a = ϕ = ϕᵧ = r = spiketimes = []
    GC.gc()
    return true
end

function send_thalamus_calculations(sessionid;
                                    structures = ["LGd-co", "LGd-sh", "LGd"], kwargs...)
    session = AN.Session(sessionid)
    probestructures = AN.getprobestructures(session, structures)
    for (probeid, structure) in probestructures
        for stimulus in ["flash_250ms", r"Natural_Images"]
            @info "Saving $(sessionid) $(structure) $(stimulus)"
            D = @dict sessionid structure stimulus
            try
                send_thalamus_calculations(D, session; kwargs...)
            catch e
                @warn e
            end
        end
    end
    return true
end

function test_calculations(args...; N = 10, kwargs...)
    i = 0
    while i < N
        X = randn(10000, 10000) # Requires at least some amount of memory (1GB?)
        sleep(12)
        i += 1
        @info "Poll: $(ENV["PBS_JOBID"])"
    end
end

function load_performance(; path = datadir("calculations"), stimulus = r"Natural_Images")
    Q = calcquality(path)[Structure = At(structures)]
    out = map([lookup(Q, Structure) |> first]) do structure
        out = map(lookup(Q, SessionID)) do sessionid
            if Q[SessionID = At(sessionid), Structure = At(structure)] == 0
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
    newsessions.sessionid = lookup(Q, SessionID) .|> Int
    return newsessions
end

function _collect_calculations(outfile; sessionid, structure, stimulus, path, subvars)
    filename = savepath((@strdict sessionid structure stimulus), "jld2", path)
    f = jldopen(filename, "r")
    begin
        # @unpack streamlinedepths, layerinfo, pass_γ, pass_θ, performance_metrics, spiketimes, trials = f
        streamlinedepths = f["streamlinedepths"]
        layerinfo = f["layerinfo"]
        pass_γ = f["pass_γ"]
        pass_θ = f["pass_θ"]
        spiketimes = f["spiketimes"]

        performance_metrics = nothing
        try
            performance_metrics = f["performance_metrics"]
        catch
            @warn "Performance metrics not found in $filename"
        end

        trials = nothing
        try
            trials = f["trials"]
        catch
            @warn "Trials table not found in $filename"
        end
    end
    ovars = Dict()
    for v in subvars
        push!(ovars, v => f[string(v)])
    end
    close(f)

    # Depth and layer info
    layernames = ToolsArray(layerinfo[1], (Depth(layerinfo[2]),))
    layernums = ToolsArray(layerinfo[3], (Depth(layerinfo[2]),))

    # Remove poor quality depth estimates
    idxs = 1:length(streamlinedepths)
    # try
    while !issorted(streamlinedepths) # Make sure depths are increasing only
        idxs = idxs[1:(end - 1)]
        streamlinedepths = streamlinedepths[idxs]
        @assert length(streamlinedepths) > 15
    end
    # catch e
    #     return nothing
    # end
    layernames = layernames[idxs]

    ovars = map(collect(ovars)) do (v, k)
        k = k[:, idxs, :] .|> Float32

        # We know the depths are sorted from above
        k = set(k, Depth(streamlinedepths))
        k = set(k,
                Depth => DimensionalData.Irregular(extrema(streamlinedepths)))

        @assert issorted(lookup(k, Depth))
        push!(k.metadata.val, :layernames => layernames)
        return v => k
    end

    D = @strdict streamlinedepths layernames pass_γ pass_θ trials sessionid performance_metrics spiketimes
    for (k, v) in pairs(D)
        outfile[string(structure) * "/" * string(sessionid) * "/" * k] = v
    end
    for (k, v) in ovars
        outfile[string(structure) * "/" * string(sessionid) * "/" * string(k)] = v .|>
                                                                                 Float32
    end
    return nothing
end

function collect_calculations(Q; path = datadir("calculations"), stimulus, rewrite = false)
    outfilepath = savepath("out", Dict("stimulus" => stimulus), "jld2", datadir())
    subvars = sort([:V, :csd, :x, :ϕ, :r, :k, :ω])
    if !isfile(outfilepath) || rewrite
        jldopen(outfilepath, "w") do outfile
            out = map(lookup(Q, Structure)) do structure
                @info "Collecting data for structure $(structure)"
                map(lookup(Q, SessionID)) do sessionid
                    if Q[SessionID = At(sessionid), Structure = At(structure),
                         stimulus = At(stimulus)] == 0
                        return nothing
                    end
                    _collect_calculations(outfile; sessionid, structure, stimulus, path,
                                          subvars)
                end
            end
        end
        # out = filter(!isnothing, out)
        # out = filter(x -> maximum(x[:streamlinedepths]) > 0.90, out) # Remove sessions
        # that don't have data to a reasonable depth --> Moved to posthoc session filter
        # if length(unique(length.(out))) != 1
        #     error("Not all structures returned the same sessions")
        #     # commonsessions = [[o[:sessionid] for o in O]
        #     #                   for O in out]
        #     # commonsessions = intersect(commonsessions...)
        #     # out = map(out) do O
        #     #     filter(o -> (o[:sessionid] in commonsessions), O)
        #     # end
        # end
    end
    if rewrite == false
        @info "Already calculated: $(stimulus). Checking quality"
        jldopen(outfilepath, "a+") do outfile
            sessionids = [keys(outfile[s]) for s in structures]
            if length(unique(sessionids)) > 1
                @warn "Not all structures returned the same sessions. Attempting to repair"
            end
            sessionids = unique(vcat(sessionids...))
            for sessionid in sessionids
                for structure in structures
                    if !haskey(outfile[structure], sessionid)
                        @warn "Missing $(structure) $(sessionid)"
                        sessionid = tryparse(Int, sessionid)
                        if Q[SessionID = At(sessionid), Structure = At(structure),
                             stimulus = At(stimulus)] == 0
                            return nothing
                        end
                        _collect_calculations(outfile; sessionid, structure, stimulus, path,
                                              subvars)
                    end
                end
            end
        end
    end
    return outfilepath
end

function load_calculations(Q; vars = sort([:V, :csd, :x, :ϕ, :r, :k, :ω]), stimulus,
                           kwargs...)
    commonkeys = [
        :streamlinedepths, :layernames, :pass_γ, :pass_θ, :trials, :sessionid,
        :performance_metrics, :spiketimes
    ]
    vars = vcat(vars, commonkeys)
    outfilepath = collect_calculations(Q; stimulus, kwargs...)
    jldopen(outfilepath, "r") do outfile
        @assert sort(keys(outfile)) == sort(structures)
        sessionids = unique([keys(outfile[s]) for s in structures]) |> only
        out = progressmap(structures) do structure
            progressmap(sessionids) do sessionid
                D = Dict()
                for v in vars
                    D[v] = outfile[structure][sessionid][string(v)]
                end
                return D
            end
        end
        return out
    end
end

function unify_calculations(Q; stimulus, vars = sort([:V, :csd, :x, :ϕ, :r, :k, :ω]),
                            rewrite = false, kwargs...)
    unifilepath = savepath("uni", Dict("stimulus" => stimulus), "jld2", datadir())
    # * Filter to posthoc sessions
    if !isfile(unifilepath) || rewrite
        @info "Loading calculations"
        out = load_calculations(Q; vars, stimulus, kwargs...)
        @info "Loaded calculations, unifying"
        session_table = load(datadir("posthoc_session_table.jld2"), "session_table")
        oursessions = session_table.ecephys_session_id
        out = map(out) do O
            filter(o -> (o[:sessionid] in oursessions), O)
        end
        @assert all(length.(out) .== length(oursessions))

        stimuli = [[DimensionalData.metadata(o[first(vars)])[:stimulus] for o in out[i]]
                   for i in eachindex(out)]
        stimulus = only(unique(vcat(stimuli...)))
        jldopen(unifilepath, "w") do outfile
            progressmap(enumerate(out)) do (si, o)
                @info "Collecting data for structure $(structures[si])"
                sessionids = getindex.(o, :sessionid)
                streamlinedepths = getindex.(o, :streamlinedepths)
                unidepths = commondepths(streamlinedepths)
                layernames = map(o) do p
                    l = p[:layernames]
                    l = set(l, Depth => p[:streamlinedepths])
                    if last(l) ∈ ["or", "scwm", "cing"] # Sometime anomalies at the boundary; we fall back on the channel structure labels, the ground truth, for confidence that this is still a cortical channel
                        l[end] = l[end - 1]
                    end
                    l[Depth = Near(unidepths)]
                end
                layernums = [ToolsArray(parselayernum.(parent(p)), dims(p))
                             for p in layernames]

                # * We want to get a unified array, plus (potentially overlapping) intervals for
                # * each layer
                layerints = map(unique(vcat(parent.(layernums)...))) do l
                    depths = map(enumerate(layernums)) do (sid, ls)
                        if sum(ls .== l) == 0
                            ma = sum(ls .== l - 1) == 0 ? 1.0 :
                                 maximum(lookup(ls[ls .== l - 1], Depth))
                            mi = sum(ls .== l + 1) == 0 ? 0.0 :
                                 minimum(lookup(ls[ls .== l + 1], Depth))
                            m = mean([ma, mi])
                            @info "Session $(sessionids[sid]) has no layer $(l)"
                            [m, m]
                        else
                            collect(extrema(lookup(ls, Depth)[ls .== l]))
                        end
                    end
                    Interval(extrema(vcat(depths...))...)
                end

                ovars = Dict()
                for (k, v) in pairs(@strdict oursessions unidepths layerints layernames layernums)
                    push!(ovars, k => v)
                end

                for v in vars
                    k = map(o) do p
                        k = p[v]
                        k = k[Depth = Near(unidepths)]
                        k = set(k, Depth => Depth(unidepths))
                        if stimulus == r"Natural_Images"
                            d = Trial(p[:trials].hit[1:size(k,
                                                            :changetime)])
                        else
                            d = Trial(1:size(k, :changetime))
                        end
                        k = set(k, Dim{:changetime} => Trial)
                        set(k, DimensionalData.format(d, lookup(k, Trial)))
                    end
                    mints = lookup.(k, 𝑡)
                    _, minti = findmin(length.(mints))
                    mints = mints[minti]
                    k = map(k) do p # Match to smallest time interval. Should only differ by a sample or so
                        p[𝑡(At(mints))]
                    end
                    if stimulus == r"Natural_Images"
                        d = Trial(vcat(lookup.(k, Trial)...))
                        k = cat(k...; dims = d)
                    else
                        k = stack(SessionID(sessionids), k)
                    end
                    push!(ovars, string(v) => k)
                end
                @info "Saving data for $(structures[si])"
                for (k, v) in pairs(ovars)
                    if eltype(v) isa Number
                        outfile[string(structures[si]) * "/" * k] = v .|> Float32
                    else
                        outfile[string(structures[si]) * "/" * k] = v
                    end
                end
            end
        end
    end
    return unifilepath
end

# function produce_out(Q::AbstractToolsArray, config; path = datadir("calculations"))
#     out = load_calculations(Q; path, stimulus = config["stimulus"], vars = config["vars"])
#     @strdict out
# end
# produce_out(Q::AbstractToolsArray; kwargs...) = config -> produce_out(Q, config; kwargs...)
# produce_out(; kwargs...) = config -> produce_out(config; kwargs...)
# function produce_out(config; path = datadir("calculations"),
#                      structures = SpatiotemporalMotifs.structures)
#     Q = calcquality(path)[Structure = At(structures)]
#     return produce_out(Q, config; path)
# end
function load_uni(; stimulus, vars = sort([:V, :csd, :x, :ϕ, :r, :k, :ω]),
                  path = datadir("calculations"),
                  structures = SpatiotemporalMotifs.structures,
                  kwargs...)
    Q = calcquality(path)[Structure = At(structures)]
    unifilepath = unify_calculations(Q; stimulus)
    @info "Loading unified data $vars for $stimulus"
    commonvars = [:oursessions, :unidepths, :layerints, :layernames, :layernums]
    jldopen(unifilepath, "r") do unifile
        @assert sort(keys(unifile)) == sort(structures)
        out = progressmap(structures) do structure
            D = Dict()
            for v in vars
                D[v] = unifile[structure][string(v)]
            end
            for v in commonvars
                D[v] = unifile[structure][string(v)]
            end
            return D
        end
        return out
    end
end

function produce_unitdepths(session::AN.AbstractSession)
    sessionid = AN.getid(session)
    @info "Computing depths for session $sessionid"
    # * Want to get a dataframe mapping all units to a probe depth and a streamline
    #   depth.
    probes = AN.getprobes(session)
    unitmetrics = AN.getunitmetrics(session)
    cdf = AN.getchannels(session)
    probedepths = Dict{Any, Any}()
    streamlinedepths = Dict{Any, Any}()
    for probe in probes.id
        _cdf = cdf[cdf.probe_id .== probe, :]
        channels = _cdf.id
        _probedepths = Dict(channels .=>
                                AN.AllenNeuropixelsBase._getchanneldepths(_cdf,
                                                                          channels;
                                                                          method = :probe))
        _streamlinedepths = Dict(channels .=>
                                     AN.AllenNeuropixelsBase._getchanneldepths(_cdf,
                                                                               channels;
                                                                               method = :streamlines))

        p = Dict(channels .=> getindex.([_probedepths], channels))
        s = Dict(channels .=> getindex.([_streamlinedepths], channels))
        merge!(probedepths, p)
        merge!(streamlinedepths, s)
    end
    D = DataFrame([keys(probedepths) |> collect, values(probedepths) |> collect,
                      getindex.([streamlinedepths], keys(probedepths))
                  ],
                  [
                      :peak_channel_id,
                      :probedepth,
                      :streamlinedepth
                  ])
    D = innerjoin(D, unitmetrics, on = :peak_channel_id)
    D = innerjoin(D, cdf, on = :peak_channel_id => :id)
    D.ecephys_session_id .= sessionid
    return D
end

function produce_unitdepths(sessionids)
    Ds = map(sessionids) do sessionid
        session = AN.Session(sessionid)
        produce_unitdepths(session)
    end
    return vcat(Ds...)
end

function load_unitdepths(Q::AbstractDimArray; path = datadir("power_spectra"))
    if !all(Q)
        @warn "The quality matrix indicates missing values. Please be cautious."
    end
    D = map(lookup(Q, Dim{:stimulus})) do stimulus
        d = map(lookup(Q, SessionID)) do sessionid
            map(lookup(Q, Structure)) do structure
                filename = savepath((@strdict sessionid structure stimulus), "jld2",
                                    path)
                f = jldopen(filename, "r")
                if haskey(f, "error")
                    return nothing
                else
                    return f["unitdepths"]
                end
            end
        end
        d = filter(!isnothing, d)
        d = vcat(d...)
    end
    D = filter(!isempty, D)
    D = vcat(D...)
    D = filter(!isnothing, D)
    D = vcat(D...)
    return D
end
function load_unitdepths(; path = datadir("power_spectra"), kwargs...)
    Q = calcquality(path)
    load_unitdepths(Q; path, kwargs...)
end
