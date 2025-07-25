using SLHQuantumSystems
using Documenter
using Literate

# Generate documentation from Literate examples BEFORE makedocs
include("generate.jl")

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
        "Examples" => GENERATED_EXAMPLE_PAGES,
        "API" => "api.md",
        "system building" => "build.md",
        "Literate Workflow" => "literate-workflow.md",
    ],
)

deploydocs(;
    repo="github.com/jeffwack/SLHQuantumSystems.jl",
    devbranch="main",
)
