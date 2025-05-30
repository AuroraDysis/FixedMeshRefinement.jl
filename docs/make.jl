using BergerOligerAMR
using Documenter

DocMeta.setdocmeta!(BergerOligerAMR, :DocTestSetup, :(using BergerOligerAMR); recursive=true)

makedocs(;
    modules=[BergerOligerAMR],
    authors="Zhen Zhong <auroradysis@gmail.com> and contributors",
    sitename="BergerOligerAMR.jl",
    format=Documenter.HTML(;
        canonical="https://AuroraDysis.github.io/BergerOligerAMR.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/AuroraDysis/BergerOligerAMR.jl",
    devbranch="main",
)
