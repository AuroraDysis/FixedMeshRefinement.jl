using FixedMeshRefinement
using Documenter

DocMeta.setdocmeta!(FixedMeshRefinement, :DocTestSetup, :(using FixedMeshRefinement); recursive=true)

makedocs(;
    modules=[FixedMeshRefinement],
    authors="Zhen Zhong <auroradysis@gmail.com> and contributors",
    sitename="FixedMeshRefinement.jl",
    format=Documenter.HTML(;
        canonical="https://AuroraDysis.github.io/FixedMeshRefinement.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/AuroraDysis/FixedMeshRefinement.jl",
    devbranch="main",
)
