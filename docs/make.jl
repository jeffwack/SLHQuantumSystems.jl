using SLHQuantumSystems
using Documenter
using DocumenterCitations
using Literate

# Generate documentation from Literate examples before calling makedocs
include("generate.jl")

DocMeta.setdocmeta!(SLHQuantumSystems, :DocTestSetup, :(using SLHQuantumSystems); recursive=true)

bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"))

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
        "Linear quantum networks" => "LinearQuantumNetworks.md",
        "API" => "api.md",
        "Examples" => GENERATED_EXAMPLE_PAGES,
        "Literate Workflow" => "literate-workflow.md",
    ],
    workdir=joinpath(@__DIR__, ".."),
    plugins=[bib],
    checkdocs=:exports,
    # Set ENV["DOCS_DRAFT"] = "true" to skip running @example blocks (fast local previews)
    draft=get(ENV, "DOCS_DRAFT", "false") == "true",
)

deploydocs(;
    repo="github.com/jeffwack/SLHQuantumSystems.jl",
    devbranch="main",
)
