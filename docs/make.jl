using Documenter, PolyFit

makedocs(;
    modules=[PolyFit],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/jishnub/PolyFit.jl/blob/{commit}{path}#L{line}",
    sitename="PolyFit.jl",
    authors="Jishnu Bhattacharya <jishnuonline@gmail.com>",
    assets=String[],
)

deploydocs(;
    repo="github.com/jishnub/PolyFit.jl",
)
