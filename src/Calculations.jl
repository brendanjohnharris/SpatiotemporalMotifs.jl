using DataFrames
using JLD2
using TimeseriesTools
using UnPack

function send_powerspectra(sessionid, stimulus, structure;
                           outpath = calcdir("power_spectra"),
                           plotpath = calcdir("plots", "power_spectra"))
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

    fstimulus = _params[:stimulus] isa Regex ? _params[:stimulus].pattern :
                _params[:stimulus]

    plotfile = joinpath(plotpath, "$(_params[:sessionid])",
                        "$(fstimulus)_$(_params[:structure]).pdf")
    pac_plotfile = joinpath(plotpath, "$(_params[:sessionid])",
                            "$(fstimulus)_$(_params[:structure])_pac.pdf")

    if isempty(ssession) # Only initialize session if we have to
        session = AN.Session(params[:sessionid])
        push!(ssession, session)
    else
        session = ssession[1]
    end
    # @info "Calculating $(stimulus) LFP in $(structure) for session $(params[:sessionid])"
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
            @info "Saved plot to `$plotfile`"
        end

        # * Calculate a comodulogram for these time series
        @info "Calculating comodulogram"
        c = x -> comodulogram(ustripall(x); fs = ustripall(samplingrate(LFP)),
                              fₚ = 3:0.25:12,
                              fₐ = 25:1:125, dp = 2, da = 20)
        N = min(360000, size(LFP, 1))
        C = map(c, eachslice(LFP[1:N, :], dims = 2))
        C = stack(dims(LFP[1:N, :], Chan), ToolsArray.(C), dims = 3)
        sLFP = mapslices(x -> surrogate(ustripall(parent(x)), IAAFT()), LFP[1:N, :],
                         dims = 1)
        sC = map(c, eachslice(sLFP, dims = 2))
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
            wsave(pac_plotfile, f)
            @info "Saved plot to `$pac_plotfile`"
        end

        # * Format data for saving
        streamlinedepths = AN.getchanneldepths(session, LFP; method = :streamlines)
        layerinfo = AN.Plots._layerplot(session, channels)

        @info "Calculating spontaneous order parameters"
        begin # * Spontaneous order parameter. Need to rectify for this
            LFP = set(LFP, Chan => Depth(depths))
            origts = deepcopy(times(LFP))
            LFP = rectify(LFP; dims = Depth)
            θ = bandpass(LFP, SpatiotemporalMotifs.THETA() .* u"Hz")
            ϕ = analyticphase(θ)
            ω = centralderiv(ϕ, dims = 𝑡, grad = phasegrad) .< 0u"Hz" # Mask negative freqs

            k = -centralderiv(ϕ, dims = Depth, grad = phasegrad)

            k[ω] .= NaN * unit(eltype(k))
            R = dropdims(nansafe(mean, dims = Depth)(sign.(k)); dims = Depth)

            k = [] # Free a bit of mem
            GC.gc()

            # * Surrogates
            idxs = randperm(size(ϕ, Depth))
            ϕs = set(ϕ, parent(ϕ[:, idxs])) # Spatially shuffle channels
            ωs = centralderiv(ϕs, dims = 𝑡, grad = phasegrad) .< 0u"Hz"
            ks = -centralderiv(ϕs, dims = Depth, grad = phasegrad)
            ks[ωs] .= NaN * unit(eltype(ks))
            sR = dropdims(nansafe(mean, dims = Depth)(sign.(ks)); dims = Depth)

            ϕs = []
            ωs = []
            ks = []
            GC.gc()
        end

        begin # * Spontaneous spike-LFP coupling. We set the LFP time indices back to their original, non-rectified values
            unitdepths = produce_unitdepths(session)
            unitdepths[!, :stimulus] .= stimulus

            @info "Calculating spike-LFP coupling"
            spiketimes = AN.getspiketimes(session, structure)
            units = AN.getunitmetrics(session)
            units = units[units.ecephys_unit_id .∈ [keys(spiketimes)], :]
            unitdepths = unitdepths[unitdepths.ecephys_unit_id .∈ [units.id], :]

            γ = bandpass(LFP, SpatiotemporalMotifs.GAMMA() .* u"Hz")
            r = abs.(hilbert(γ))
            r[ω] .= NaN * unit(eltype(r))
            ϕ[ω] .= NaN * unit(eltype(ϕ))

            ω = [] # Free mem
            GC.gc()

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
                    "unitdepths" => unitdepths,
                    "plotfiles" => relpath.([plotfile, pac_plotfile], [projectdir()]))
        tagsave(outfile, outD)

        GC.safepoint()
        GC.gc()
    catch e
        GC.safepoint()
        GC.gc()
        @warn e
        tagsave(outfile, Dict("error" => sprint(showerror, e)))
    end
    # @info "Finished calculations"
    GC.safepoint()
    GC.gc()
    return outfile
end

function joint_rectify(LFP, spiketimes, stimulustimes)
    orig_ts = times(LFP) # Unaligned times
    LFP = rectify(LFP, dims = 𝑡, tol = 3)

    _tmap = Timeseries(orig_ts, zeros(length(orig_ts)))
    intrvl = Interval(_tmap)

    tmap = deepcopy(_tmap)
    tmap[𝑡(Near(stimulustimes))] .= stimulustimes
    tmap = rectify(tmap, dims = 𝑡, tol = 3)
    @assert times(LFP) == times(tmap)
    # dev = tmap[findlast(tmap .> 1)] - times(tmap)[findlast(tmap .> 1)] # Alignment error
    stimulustimes = times(tmap[tmap .> 1]) |> collect # * The change times transformed to rectified time coordinates

    spiketimes = map(collect(spiketimes)) do (u, ts) # * Rectify the spike times
        ts = ts[ts .∈ [intrvl]] # Remove spikes outside the LFP interval
        tmap = deepcopy(_tmap)
        tmap[𝑡(Near(ts))] .= ts
        tmap = rectify(tmap, dims = 𝑡, tol = 3)
        ts = times(tmap[tmap .> 1]) |> collect
        return u => Float32.(ts)
    end |> Dict

    stimulustimes = stimulustimes .* u"s"
    tmap = set(tmap .* u"s", 𝑡 => times(tmap) .* u"s")

    return LFP, stimulustimes, spiketimes # Aligned to common rectified times
end

function what_files(stimulus)
    if stimulus == r"Natural_Images"
        files = [stimulus, "Natural_Images_nochange"] # * and "Natural_Images_omission"
    elseif stimulus == "Natural_Images_passive"
        files = [stimulus, "Natural_Images_passive_nochange"]
    else
        files = [stimulus]
    end
    return files
end
function what_lfp(stimulus, session, structure)
    probeid = AN.getprobe(session, structure)
    if stimulus == r"Natural_Images"
        _stimulus = r"Natural_Images"
        epoch = (:longest, :active)
    elseif stimulus == "Natural_Images_passive"
        _stimulus = r"Natural_Images"
        epoch = (:longest, :active => ByRow(==(false)))
    else # Flashes or others
        _stimulus = stimulus
        epoch = :longest
    end
    LFP = AN.formatlfp(session;
                       probeid, structure, stimulus = _stimulus, rectify = false, epoch)
    return LFP
end

function get_stimulustimes(stimulus, session)
    if stimulus == r"Natural_Images"
        trials = AN.getchangetrials(session)
        stimulustimes = trials.change_time_with_display_delay
    elseif stimulus == "Natural_Images_omission"
        stimuli = AN.getstimuli(session)
        stimuli = stimuli[.!ismissing.(stimuli.omitted), :]
        stimuli = stimuli[identity.(stimuli.omitted), :] # Just the omission trials
        stimuli = stimuli[stimuli.active, :] # * just the active-epoch omissions
        stimulustimes = stimuli.start_time
        trials = stimuli
    elseif stimulus == "Natural_Images_nochange"
        stimuli = AN.getstimuli(session)
        stimuli = stimuli[.!ismissing.(stimuli.omitted), :]
        stimuli = stimuli[.!stimuli.omitted, :]
        stimuli = stimuli[stimuli.active, :] # * just the active-epoch omissions
        stimuli = stimuli[.!ismissing.(stimuli.is_change), :]
        stimuli = stimuli[.!stimuli.is_change, :]
        stimulustimes = stimuli.start_time
        trials = stimuli
    elseif stimulus == "Natural_Images_passive"
        stimuli = AN.getstimuli(session)
        # stimuli = stimuli[stimuli.start_time .∈ [Interval(LFP)], :] # !!!
        stimuli = stimuli[.!ismissing.(stimuli.is_change), :]
        stimuli = stimuli[identity.(stimuli.is_change), :]
        stimuli = stimuli[.!stimuli.active, :] # * just the passive changes
        stimulustimes = stimuli.start_time
        trials = stimuli
    elseif stimulus == "Natural_Images_passive_nochange"
        stimuli = AN.getstimuli(session)
        stimuli = stimuli[.!ismissing.(stimuli.omitted), :]
        stimuli = stimuli[.!stimuli.omitted, :]
        stimuli = stimuli[.!stimuli.active, :] # * just the passive-epoch omissions
        stimuli = stimuli[.!ismissing.(stimuli.is_change), :]
        stimuli = stimuli[.!stimuli.is_change, :]
        stimulustimes = stimuli.start_time
        trials = stimuli
    else # Flashes and others
        trials = AN.stimulusintervals(session, stimulus)
        stimulustimes = trials.start_time # For flashes
    end
    return stimulustimes, trials
end

function compute_csd(LFP, stimulustimes)
    # CSD computed as
    prestim = LFP[𝑡 = -0.1u"s" .. 0u"s"]
    prestim_mean = mean(prestim, dims = (𝑡, :changetime))
    baseline_LFP = LFP .- prestim_mean # Baseline correct
    csd = .-centralderiv(centralderiv(baseline_LFP, dims = Depth); dims = Depth)
    return csd
end

function aligned_calculations(LFP; pass_θ, pass_γ, ΔT, doupsample, stimulustimes)
    θ = bandpass(LFP, pass_θ)
    doupsample > 0 && (θ = upsample(θ, doupsample, Depth))

    γ = bandpass(LFP, pass_γ)
    doupsample > 0 && (γ = upsample(γ, doupsample, Depth))

    a = hilbert(θ)
    aᵧ = hilbert(γ)

    ϕ = angle.(a) # _generalized_phase(ustrip(θ)) #
    # ϕᵧ = angle.(aᵧ) # _generalized_phase(ustrip(γ)) #

    ω = centralderiv(ϕ, dims = 𝑡, grad = phasegrad) # Angular Frequency

    k = -centralderiv(ϕ, dims = Depth, grad = phasegrad) # Wavenumber
    # We have a minus sign because the hilbert transform uses the convention of ϕ = ωt (for
    # univariate; phase increases with time for positive frequencies), which implies ϕ = ωt - kx (for multivariate).
    v = ω ./ k # Phase velocity of theta

    # ωᵧ = centralderiv(ϕᵧ, dims = 𝑡, grad = phasegrad) # Frequency
    # kᵧ = -centralderiv(ϕᵧ, dims = Depth, grad = phasegrad) # Wavenumber
    # ∂ωᵧ = centralderiv(ωᵧ; dims = 𝑡)
    # ∂kᵧ = centralderiv(kᵧ; dims = 𝑡)
    # vₚ = ωᵧ ./ kᵧ # Phase velocity
    # vᵧ = ∂ωᵧ ./ ∂kᵧ # Group velocity

    function alignmatchcat(x)
        x = align(x, stimulustimes, ΔT)[2:(end - 2)]
        x = matchdim(x; dims = 1)
        x = stack(Dim{:changetime}(times(x)), x)
    end

    V = LFP
    output = Dict{String, Any}()
    @pack! output = V, θ, γ, a, ϕ, ω, k, v, aᵧ

    # * Clean up and shrink size
    for k in keys(output)
        output[k] = alignmatchcat(output[k])
    end
    output["r"] = abs.(output["aᵧ"]) .* unit(eltype(output["V"]))

    csd = compute_csd(output["V"], stimulustimes)
    output["csd"] = csd

    for k in keys(output)
        if eltype(ustrip(output[k][1])) <: Complex
            output[k] = ComplexF32.(output[k])
        elseif eltype(ustrip(output[k][1])) <: Real
            output[k] = Float32.(output[k])
        else
            throw(error("Unexpected eltype for $k: $(eltype(output[k]))"))
        end
    end

    return output
end

function render_complete(outfiles, complete)
    map(complete, outfiles) do c, outfile
        _, params, _ = parse_savename(outfile; connector)
        desc = " | $(params["sessionid"]) | $(params["stimulus"]) | $(params["structure"])"
        if c
            printstyled("\n✓" * desc; color = :green)
        else
            printstyled("\n✗" * desc; color = :yellow)
        end
    end
    println("\n")
    return all(complete)
end

function send_calculations(D::Dict, session = AN.Session(D["sessionid"]);
                           pass_θ = THETA() .* u"Hz",
                           pass_γ = GAMMA() .* u"Hz",
                           ΔT = -0.5u"s" .. 0.75u"s",
                           doupsample = 0,
                           rewrite = false,
                           outpath)
    @unpack sessionid, structure, stimulus = D

    # * Infer outfiles
    files = what_files(stimulus)
    outfiles = map(files) do file
        _D = @strdict sessionid structure
        _D["stimulus"] = file
        return savepath(_D, "jld2", outpath)
    end

    complete = map(zip(files, outfiles)) do (file, outfile)
        rewrite && return false # * If we are rewriting, we always recalculate
        if !isfile(outfile)
            @info "Creating `$(outfile)`"
            return false
        end
        try
            jldopen(outfile, "r"; iotype = IOStream) do fl
                # * Check error
                if haskey(fl, "error")
                    @warn "Error in `$(outfile)`"
                    @error fl["error"]
                    return false
                end

                # * Check required keys
                if !has_calc_keys(fl)
                    @error "Missing keys in $(outfile)"
                    return false
                end

                # * Check loading
                try
                    fl["performance_metrics"]
                catch e
                    @warn "Error loading from $(outfile) "
                    @error e
                    return false
                end

                return true
            end
        catch e
            @warn "Error opening $(outfile)"
            @error e
            return false
        end
    end
    if !render_complete(outfiles, complete)
        files = files[.!complete]
        outfiles = outfiles[.!complete]

        performance_metrics = AN.getperformancemetrics(session)

        # * Format the LFP for this stimulus and any files we want to save
        LFP = what_lfp(stimulus, session, structure)
        depths = AN.getchanneldepths(session, LFP; method = :probe)
        streamlinedepths = AN.getchanneldepths(session, LFP; method = :streamlines)
        channels = lookup(LFP, Chan)
        LFP = set(LFP, Chan => Depth(depths))
        LFP = set(LFP, Depth => lookup(LFP, Depth) .* u"μm")
        LFP = rectify(LFP, dims = Depth)

        layerinfo = AN.Plots._layerplot(session, channels)
        spiketimes = AN.getspiketimes(session, structure)
        units = AN.getunitmetrics(session)
        units = units[units.ecephys_unit_id .∈ [keys(spiketimes)], :]

        for (file, outfile) in zip(files, outfiles)
            try
                stimulustimes, trials = get_stimulustimes(file, session)

                LFP_J, stimulustimes_J, spiketimes_J = joint_rectify(LFP, spiketimes,
                                                                     stimulustimes)

                trials.rectified_change_times = stimulustimes_J

                # * Format LFP dims
                LFP_J = set(LFP_J, 𝑡 => times(LFP_J) .* u"s")

                aligned_outputs = aligned_calculations(LFP_J; pass_θ, pass_γ, ΔT,
                                                       doupsample,
                                                       stimulustimes = stimulustimes_J)

                out = aligned_outputs
                out["trials"] = trials[2:(end - 2), :]

                out["channels"] = channels
                out["streamlinedepths"] = streamlinedepths
                out["units"] = units
                out["spiketimes"] = spiketimes_J # ! The RECTIFIED spiketimes
                out["layerinfo"] = layerinfo
                out["pass_θ"] = pass_θ
                out["pass_γ"] = pass_γ
                out["performance_metrics"] = performance_metrics # Can't save dicts with compression reliably, beware!

                # * Save
                @info "Saving `$outfile`"
                @tagsave outfile out
            catch e
                @error "Error saving `$outfile`"
                @error e
                tagsave(outfile, Dict("error" => sprint(showerror, e)))
                continue
            end
        end

        out = LFP = LFP_J = aligned_outputs = []
        GC.gc()
    end
    return outfiles
end

function send_calculations(sessionid;
                           structures = ["VISp",
                               "VISl",
                               "VISrl",
                               "VISal",
                               "VISpm",
                               "VISam"],
                           stimuli = [r"Natural_Images",
                               "flash_250ms",
                               "Natural_Images_passive"],
                           outpath = calcdir("calculations"),
                           kwargs...)
    session = AN.Session(sessionid)
    probestructures = AN.getprobestructures(session, structures)
    N = length(probestructures) * length(stimuli)
    for (i, (probeid, structure)) in enumerate(probestructures)
        for (j, stimulus) in enumerate(stimuli)
            n = j + (i - 1) * length(stimuli)
            printstyled("\n\nCalculating $(sessionid), $(structure), $(stimulus) ($n/$N)\n",
                        color = :green)
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

function test_calculations(args...; N = 10, kwargs...)
    i = 0
    while i < N
        X = randn(10000, 10000) # Requires at least some amount of memory (1GB?)
        sleep(12)
        i += 1
        @info "Poll: $(ENV["PBS_JOBID"])"
    end
end

function load_performance(; path = calcdir("calculations"), stimulus = r"Natural_Images")
    data, pfile = produce_or_load(Dict("stimulus" => stimulus),
                                  calcdir();
                                  filename = savepath("performance")) do conf
        stimulus = conf["stimulus"]
        Q = calcquality(path)[Structure = At(structures)]
        out = map([lookup(Q, Structure) |> first]) do structure
            out = map(lookup(Q, SessionID)) do sessionid
                if Q[SessionID = At(sessionid), Structure = At(structure)] == 0
                    return nothing
                end
                filename = savepath((@strdict sessionid structure stimulus), "jld2", path)
                f = jldopen(filename, "r"; iotype = IOStream)
                performance_metrics = f["performance_metrics"]
                close(f)
                return performance_metrics
            end
            out = filter(!isnothing, out)
        end
        performance = only(out)
        newsessions = performance |> DataFrame
        newsessions.sessionid = lookup(Q, SessionID) .|> Int
        return Dict("performance" => newsessions)
    end
    display(pfile)
    return data["performance"]
end

function _collect_calculations(; sessionid, structure, stimulus, path, subvars)
    outdict = Dict{String, Any}()
    filename = savepath((@strdict sessionid structure stimulus), "jld2", path)
    try
        jldopen(filename, "r"; iotype = IOStream) do f
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
                outdict[string(k)] = v
            end
            for (k, v) in ovars
                outdict[string(k)] = v .|> Float32
            end
        end
    catch e
        @error "Error in collecting $(filename)"
        @error e
        outdict["error"] = sprint(showerror, e)
    end
    return outdict
end

recursive_merge(x::AbstractDict...) = merge(recursive_merge, x...)

function collect_calculations(Q; path = calcdir("calculations"), stimulus,
                              rewrite = false,
                              outpath = calcdir())
    sessionids = lookup(Q, SessionID) .|> string
    outfilepath = savepath("out", Dict("stimulus" => stimulus), "jld2", outpath)
    if contains(stimulus |> string, "nochange")
        subvars = sort([:csd]) # Only use the LFP for inferring the CSD
    else
        subvars = sort([:V, :csd, :θ, :ϕ, :r, :k, :ω])
    end
    if isfile(outfilepath)
        good = jldopen(outfilepath, "r"; iotype = IOStream) do outfile
            good = all(structures .∈ [keys(outfile)])
            if !good
                @warn "Not all structures found in $(outfilepath). Recalculating."
                return good
            end
            # * Check if all structures have sessions
            for structure in structures
                grp = outfile.root_group[structure]
                _good = isempty(setdiff(sessionids, keys(grp)))
                if !_good
                    @warn "Missing session IDs in $(structure) for `$(stimulus)`"
                    good = good && _good
                    break
                end
                for sessionid in sessionids
                    subgrp = grp[sessionid]
                    _good = all(string.(subvars) .∈ [keys(subgrp)])
                    if !_good
                        @warn "Missing variables $(subvars) in $(sessionid) for $(structure) during `$(stimulus)`"
                        good = good && _good
                        break
                    end
                end
            end
            return good
        end
    else
        good = false
    end
    if rewrite && !good
        @info "Rewriting collected calculations for `$(stimulus)` from `$path`"
        rm(outfilepath, force = true)
    end

    if !isfile(outfilepath) # * This uses a lot of memory; use ~200 to be safe.
        @info "Collecting `$stimulus` from `$path`"
        begin # * Construct groups
            jldopen(outfilepath, "a+") do outfile
                map(lookup(Q, Structure)) do structure
                    gs = JLD2.Group(outfile, string(structure);
                                    est_num_entries = size(Q, SessionID),
                                    est_link_name_len = 10)
                    map(lookup(Q, SessionID)) do sessionid
                        g = JLD2.Group(gs, string(sessionid);
                                       est_num_entries = size(Q, SessionID))
                    end
                end
                return nothing
            end
        end
        @progress for structure in lookup(Q, Structure)
            @info "Collecting data for structure $(structure)"
            out = Vector{Dict{String, Any}}(undef, size(Q, SessionID))
            Threads.@threads for (i, sessionid) in collect(enumerate(lookup(Q, SessionID)))
                if hasdim(Q, :stimulus)
                    q = Q[SessionID = At(sessionid), Structure = At(structure),
                          stimulus = At(stimulus)]
                else
                    q = Q[SessionID = At(sessionid), Structure = At(structure)]
                end
                if q == 0
                    throw(error("Bad-quality data for $(sessionid), $(structure), $(stimulus)"))
                end
                out[i] = _collect_calculations(; sessionid, structure, stimulus, path,
                                               subvars)
            end
            # * Write groups recursively
            jldopen(outfilepath, "a+") do outfile
                for (i, sessionid) in enumerate(lookup(Q, SessionID))
                    _out = out[i]
                    for k in keys(_out)
                        key = join([string(structure), string(sessionid), string(k)], '/')
                        outfile[key] = _out[k]
                    end
                end
            end
            GC.gc()
        end
        @info "Finished collecting data, saved to `$outfilepath`"
    end
    return outfilepath
end

function load_calculations(Q; vars = sort([:V, :csd, :θ, :ϕ, :r, :k, :ω]), stimulus,
                           kwargs...)
    commonkeys = [
        :streamlinedepths, :layernames, :pass_γ, :pass_θ, :trials, :sessionid,
        :performance_metrics, :spiketimes
    ]
    vars = vcat(vars, commonkeys)
    outfilepath = collect_calculations(Q; stimulus, kwargs...)
    sessionids = jldopen(outfilepath, "r"; iotype = IOStream) do outfile
        grp = outfile.root_group
        @assert sort(keys(outfile)) == sort(structures)
        sessionids = map(structures) do structure
                         keys(grp[string(structure)]) |> sort
                     end |> unique |> only
        sessionids = sessionids[indexin(string.(lookup(Q, SessionID)), sessionids)]
        return sessionids
    end

    out = Vector{Vector{Dict{Symbol, Any}}}(undef, length(structures))
    @withprogress name="Loading calculations" begin
        threadlog = Threads.Atomic{Int}(0)
        threadmax = length(structures) * length(sessionids)
        lk = Threads.ReentrantLock()
        for (s, structure) in collect(enumerate(structures))
            out[s] = Vector{Dict{Symbol, Any}}(undef, length(sessionids))
            jldopen(outfilepath, "r"; iotype = IOStream) do outfile
                for (i, sessionid) in enumerate(sessionids)
                    out[s][i] = Dict{Symbol, Any}()

                    for v in vars
                        k = join([string(structure), string(sessionid), string(v)], '/')
                        out[s][i][v] = outfile[k]
                    end
                    lock(lk) do
                        Threads.atomic_add!(threadlog, 1)
                        @logprogress threadlog[] / threadmax
                    end
                end
            end
        end
    end
    return out
end

function unify_calculations(Q; stimulus, vars = sort([:V, :csd, :θ, :ϕ, :r, :k, :ω]),
                            rewrite = false, outpath = calcdir(), kwargs...)
    unifilepath = savepath("uni", Dict("stimulus" => stimulus), "jld2", outpath)
    # * Filter to posthoc sessions
    if !isfile(unifilepath) || rewrite
        @info "Loading calculations"
        out = load_calculations(Q; vars, stimulus, outpath, kwargs...)
        @info "Loaded calculations, unifying"
        session_table = load(calcdir("posthoc_session_table.jld2"), "session_table")
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
                @info "Unifying data for structure $(structures[si])"
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

function load_uni(; stimulus, vars = sort([:V, :csd, :θ, :ϕ, :r, :k, :ω]),
                  path = calcdir("calculations"),
                  structures = SpatiotemporalMotifs.structures,
                  kwargs...)
    Q = calcquality(path)[Structure = At(structures)]
    unifilepath = unify_calculations(Q; stimulus, kwargs...)
    @info "Loading unified data $vars for $stimulus"
    commonvars = [:oursessions, :unidepths, :layerints, :layernames, :layernums]
    jldopen(unifilepath, "r"; iotype = IOStream) do unifile
        @assert sort(keys(unifile)) == sort(structures)
        out = progressmap(structures) do structure
            D = Dict()
            for v in vars
                k = join([string(structure), string(v)], '/')
                D[v] = unifile[k]
            end
            for v in commonvars
                k = join([string(structure), string(v)], '/')
                D[v] = unifile[k]
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

function load_unitdepths(Q::AbstractDimArray; path = calcdir("power_spectra"))
    if !all(Q)
        @warn "The quality matrix indicates missing values. Please be cautious."
    end
    D = map(lookup(Q, Dim{:stimulus})) do stimulus
        d = map(lookup(Q, SessionID)) do sessionid
            map(lookup(Q, Structure)) do structure
                filename = savepath((@strdict sessionid structure stimulus), "jld2",
                                    path)
                jldopen(filename, "r") do f
                    if haskey(f, "error")
                        return nothing
                    else
                        return f["unitdepths"]
                    end
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
function load_unitdepths(; path = calcdir("power_spectra"), kwargs...)
    Q = calcquality(path)
    load_unitdepths(Q; path, kwargs...)
end
