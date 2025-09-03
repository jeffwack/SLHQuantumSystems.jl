using Documenter, Literate

# Define paths following LiveServer + Literate recommended pattern
LITERATE_INPUT = joinpath(@__DIR__, "..", "examples")
LITERATE_OUTPUT = joinpath(@__DIR__, "src", "generated")

# Process literate files
generated_pages = []
for (root, _, files) ∈ walkdir(LITERATE_INPUT), file ∈ files
    splitext(file)[2] == ".jl" || continue
    ipath = joinpath(root, file)
    opath = splitdir(replace(ipath, LITERATE_INPUT=>LITERATE_OUTPUT))[1]
    
    # Generate markdown for documentation
    Literate.markdown(ipath, opath, documenter=true, execute=false)
    
    # Add to pages list
    base_name = splitext(file)[1]
    title = titlecase(replace(base_name, r"([a-z])([A-Z])" => s"\1 \2"))
    relative_path = joinpath("generated", "$(base_name).md")
    push!(generated_pages, title => relative_path)
end

# Export the generated pages for use in make.jl
const GENERATED_EXAMPLE_PAGES = generated_pages
