#! /bin/bash
# -*- mode: julia -*-
#=
exec julia +1.10.10 -t auto --color=yes "${BASH_SOURCE[0]}" "$@"
=#

using DrWatson
@quickactivate "SpatiotemporalMotifs"

using SpatiotemporalMotifs
@preamble

begin # * Load quality
    path = calcdir("calculations")
    Q = calcquality(path)[Structure = At(structures)]
end

begin # * Check quality
    @info "Checking quality of calculations"
    outpath = datadir("calculations")
    map(lookup(Q, :stimulus), eachslice(Q, dims = :stimulus)) do stimulus, q
        qual = mean(q)
        if qual < 1
            @error "Quality check failed for stimulus: `$stimulus`, quality: $qual. Attempting to repair."
            badsessions = lookup(q, SessionID)[last.(findall(.!q) .|> Tuple)] |> unique
            if length(badsessions) > 5
                throw(error("Too many bad sessions: $(length(badsessions)) for stimulus: `$stimulus`. Exiting."))
                return
            else
                @warn "Attempting to repair bad sessions: $badsessions"
            end
            SpatiotemporalMotifs.send_calculations.(badsessions; outpath)
        end
    end
end

begin # * Run collect_calculations for each session
    map(lookup(Q, :stimulus), eachslice(Q, dims = :stimulus)) do stimulus, q
        out = load_calculations(q; stimulus, vars = [:csd], rewrite = true)
        out = []
        GC.gc()
    end
end

# * Now do unified data
begin
    map(lookup(Q, :stimulus), eachslice(Q, dims = :stimulus)) do stimulus, q
        if !contains(stimulus |> string, "nochange")
            uni = load_uni(; stimulus)
        end
        uni = []
        GC.gc()
    end
end

begin # * Load and save performance metrics
    path = calcdir("calculations")
    performance = load_performance(; path)
end
