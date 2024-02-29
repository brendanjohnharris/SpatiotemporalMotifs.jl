using SpatiotemporalMotifs
using Documenter

DocMeta.setdocmeta!(SpatiotemporalMotifs, :DocTestSetup, :(using SpatiotemporalMotifs); recursive=true)

makedocs(;
    modules=[SpatiotemporalMotifs],
    authors="brendanjohnharris <brendanjohnharris@gmail.com> and contributors",
    sitename="SpatiotemporalMotifs.jl",
    format=Documenter.HTML(;
        canonical="https://brendanjohnharris.github.io/SpatiotemporalMotifs.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/brendanjohnharris/SpatiotemporalMotifs.jl",
    devbranch="main",
)
