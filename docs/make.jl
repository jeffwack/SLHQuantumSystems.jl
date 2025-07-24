using SLHQuantumSystems
using Documenter

DocMeta.setdocmeta!(SLHQuantumSystems, :DocTestSetup, :(using SLHQuantumSystems); recursive=true)

makedocs(;
    modules=[SLHQuantumSystems],
    authors="Jeffrey Wack <jeffwack111@gmail.com> and contributors",
    sitename="SLHQuantumSystems.jl",
    format=Documenter.HTML(;
        canonical="https://jeffwack.github.io/SLHQuantumSystems.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "API" => "api.md",
        "system building" => "build.md",
    ],
)

deploydocs(;
    repo="github.com/jeffwack/SLHQuantumSystems.jl",
    devbranch="main",
)
