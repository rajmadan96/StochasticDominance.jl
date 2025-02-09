using StochasticDominance
using Documenter

DocMeta.setdocmeta!(StochasticDominance, :DocTestSetup, :(using StochasticDominance); recursive=true)

makedocs(;
    modules=[StochasticDominance],
    authors="Rajmadan Lakshmanan",
    sitename="StochasticDominance.jl",
    format=Documenter.HTML(;
        canonical="https://rajmadan96.github.io/StochasticDominance.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/rajmadan96/StochasticDominance.jl.git",
    devbranch="main"
)


