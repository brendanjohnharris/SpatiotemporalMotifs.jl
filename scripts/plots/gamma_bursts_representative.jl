#! /bin/bash
#=
exec julia -t auto "${BASH_SOURCE[0]}" "$@"
=#
using DrWatson
@quickactivate "SpatiotemporalMotifs"

import TimeseriesTools: freqs
using FileIO
using Unitful
import DimensionalData: metadata
using MultivariateStats
import SpatiotemporalMotifs: plotdir, calcquality, layernum2name, savepath,
                             structurecolors, amplitudecolormap, structures, commondepths,
                             parselayernum,
                             load_calculations, unify_calculations,
                             plotlayerints!, plotlayermap!,
                             @preamble
@preamble
set_theme!(foresight(:physics))

stimulus = r"Natural_Images"
vars = [:ϕ, :aᵧ]
alpha = 0.7

session_table = load(datadir("session_table.jld2"), "session_table")
oursessions = session_table.ecephys_session_id

path = datadir("calculations")
Q = calcquality(path)[structure = At(structures)]
quality = mean(Q[stimulus = At(stimulus)])

config = @strdict stimulus vars
data, file = produce_or_load(config, datadir(); filename = savepath) do config
    out = load_calculations(Q; path, stimulus, vars)
    uni = unify_calculations(out; vars)
    @strdict out uni
end
out = data["out"]
uni = data["uni"]
GC.gc()

begin # * Single-trial LFP and wavenumber
    f = FourPanel()
    mouse = 13
    trial = 20
    structs = [1, 4]
    for (_i, i) in enumerate(structs)
        k = out[i][mouse][:ϕ][:, :, trial]
        x = abs.(out[i][mouse][:aᵧ][:, :, trial]) .* 1000 # mV to μV
        layernames = out[i][mouse].layernames
        k = set(k, Depth(lookup(layernames, :depth)))
        x = set(x, Depth(lookup(layernames, :depth)))

        ax = Axis(f[1, _i], yreversed = true)
        structure = metadata(x)[:structure]
        ax.title = structure
        p1 = plotlayermap!(ax, x, layernames; colormap = amplitudecolormap,
                           colorrange = (0, 0.1)) |> first

        ax = Axis(f[2, _i], yreversed = true)
        structure = metadata(k)[:structure]
        ax.title = structure
        p2 = plotlayermap!(ax, k, layernames; colormap = cyclic) |> first

        if _i == length(structs)
            c = Colorbar(f[1, _i + 1], p1)
            c.label = "γ amplitude (mV)"
            c = Colorbar(f[2, _i + 1], p2)
            c.label = "θ phase"
        end
    end
    addlabels!(f)
    display(f)
    wsave(plotdir("gamma_bursts_representative", "single_trial_lfp.pdf"), f)
end
