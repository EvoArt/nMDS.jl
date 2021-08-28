using nMDS
using Documenter

DocMeta.setdocmeta!(nMDS, :DocTestSetup, :(using nMDS); recursive=true)

makedocs(;
    modules=[nMDS],
    authors="Arthur Newbury",
    repo="https://github.com/EvoArt/nMDS.jl/blob/{commit}{path}#{line}",
    sitename="nMDS.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://EvoArt.github.io/nMDS.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/EvoArt/nMDS.jl",
)
