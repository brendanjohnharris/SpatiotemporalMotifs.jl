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

begin # * Load the trial LFP for natural images
    session_table = load(datadir("posthoc_session_table.jld2"), "session_table")
    oursessions = session_table.ecephys_session_id
    path = datadir("calculations")
    Q = calcquality(path)[Structure = At(structures)]
    Q = Q[SessionID(At(oursessions))]
    @assert mean(Q[stimulus = At(r"Natural_Images")]) == 1
    out = load_calculations(Q; stimulus = r"Natural_Images", vars = [:V])
end

begin # * Calculate 1/f exponent for each trial
    addprocs(15)
    @everywhere using DrWatson
    @everywhere DrWatson.@quickactivate "SpatiotemporalMotifs"
    @everywhere using Unitful
    @everywhere using TimeseriesTools
    @everywhere using SpatiotemporalMotifs
end

begin # * 1/f across all sessions. Should take about 15 minutes over 50 workers
    oneoneff = map(out) do out_structure
        @info "Calculating 1/f exponent for $(metadata(out_structure[1][:V])[:structure])"
        pmap(out_structure) do out_session
            V = out_session[:V]
            V = V[𝑡 = SpatiotemporalMotifs.INTERVAL]
            map(SpatiotemporalMotifs.trial_1onf, eachslice(V, dims = (:Depth, :changetime)))
        end
    end
end
begin
    χ = map(oneoneff) do oneoneff_structure
        map(oneoneff_structure) do oneoneff_session
            getindex.(oneoneff_session, :χ)
        end
    end
    b = map(oneoneff) do oneoneff_structure
        map(oneoneff_structure) do oneoneff_session
            getindex.(oneoneff_session, :b)
        end
    end
    D = @strdict χ b
    tagsave(datadir("trial_1onf.jld2"), D)
end

begin # * Load data
    D = load(datadir("trial_1onf.jld2"))
    χ = D["χ"]
    b = D["b"]
    @assert length(χ) == length(b)
end

begin
    hitmiss = map(out) do structure
        map(structure) do session
            hitmiss = session[:trials].hit
        end
    end
end

begin # * Bin to common depths
    allsessions = getindex.(out[1], :sessionid)
    bins = intervals(0.0:0.1:0.9)

    Δχ = map(χ, hitmiss) do χ_structure, hm_structure
        x = map(χ_structure, hm_structure) do χ_session, hm_session
            a = χ_session[:, hm_session .== 1]
            b = χ_session[:, hm_session .== 0]
            d = median(a, dims = 2) .- median(b, dims = 2)
            y = dropdims(d, dims = 2)
            y = groupby(y, Depth => Bins(bins))
            y = mean.(y)
            set(y, Depth => mean.(lookup(y, Depth)))
        end
        ToolsArray(x, (SessionID(allsessions),)) |> stack
    end

    Δχs = map(χ, hitmiss, structures) do χ_structure, hm_structure, structure
        @info "Calculating shuffled Δχ for $(structure)"
        x = map(χ_structure, hm_structure) do χ_session, hm_session
            repeats = 1000
            z = map(1:repeats) do _
                # * Shuffle the hit/miss labels
                idxs = randperm(length(hm_session))
                hm_session = hm_session[idxs]
                a = χ_session[:, hm_session .== 1]
                b = χ_session[:, hm_session .== 0]
                d = median(a, dims = 2) .- median(b, dims = 2)
                y = dropdims(d, dims = 2)
                y = groupby(y, Depth => Bins(bins))
                y = mean.(y)
                set(y, Depth => mean.(lookup(y, Depth)))
            end
            ToolsArray(z, (Dim{:repeats}(1:repeats),)) |> stack
        end
        ToolsArray(x, (SessionID(allsessions),)) |> stack
    end

    𝑝 = map(Δχ, Δχs) do Δχ, Δχs
        ps = map(eachslice(Δχ, dims = Depth), eachslice(Δχs, dims = Depth)) do x, y
            # * Perform the Mann-Whitney U test
            x = median(x)
            y = median(y, dims = SessionID)
            n_extreme = count(abs.(y) .>= abs(x))
            p = (n_extreme + 1) / (length(y) + 1)
            return p
            # p = pvalue(MannWhitneyUTest(x[:], y[:]); tail = :both)
        end
    end
    𝑝 = ToolsArray(𝑝, (Structure(structures),)) |> stack
    𝑝[:] .= MultipleTesting.adjust(𝑝[:], BenjaminiHochberg())
end

begin # * Plotting
    # * Filter for oursessions....
    f = Figure()
    ax = Axis(f[1, 1])

    file = map(structures, Δχ) do structure, Δχ_structure
        Δχ_session = Δχ_structure[SessionID(At(oursessions))]
        μ, (σl, σh) = bootstrapmedian(Δχ_session, dims = 2)

        _μ = upsample(μ, 5)
        _σl, _σh = upsample(σl, 5), upsample(σh, 5)

        band!(ax, lookup(_μ, 1), _σl, _σh, label = structure,
              color = (structurecolormap[structure], 0.2))
        lines!(ax, _μ, label = structure, color = structurecolormap[structure])

        sigs = 𝑝[Structure = At(structure)] .< SpatiotemporalMotifs.PTHR

        scatter!(ax, μ[sigs], color = structurecolormap[structure], label = structure)

        scatter!(ax, μ[.!sigs],
                 strokecolor = structurecolormap[structure], label = structure,
                 color = :white, strokewidth = 4)
        return Δχ_session
    end
    display(f)
end
