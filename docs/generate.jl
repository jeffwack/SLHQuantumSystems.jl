using Literate

# Define paths
EXAMPLEDIR = joinpath(@__DIR__, "..", "examples")
BUILDDIR = joinpath(@__DIR__, "src", "generated")

# Ensure the build directory exists
mkpath(BUILDDIR)

println("Generating documentation from Literate examples...")

# Automatically discover all .jl files in examples directory
example_files = filter(f -> endswith(f, ".jl"), readdir(EXAMPLEDIR))
generated_pages = []

# Process each example file
for example in example_files
    example_path = joinpath(EXAMPLEDIR, example)
    
    if isfile(example_path)
        println("Processing: $example")
        
        # Generate markdown for documentation
        Literate.markdown(
            example_path, 
            BUILDDIR; 
            documenter=true,  # Enable Documenter.jl integration
            execute=false     # Don't execute during docs build for reliability
        )
        
        # Create page entry for this example
        # Convert filename to title (e.g., "cascadedcavities.jl" -> "Cascaded Cavities")
        base_name = splitext(example)[1]
        title = titlecase(replace(base_name, r"([a-z])([A-Z])" => s"\1 \2"))
        page_path = "generated/$(base_name).md"
        
        push!(generated_pages, title => page_path)
        
        println("  → Generated markdown: $page_path")
    else
        @warn "Example file not found: $example_path"
    end
end

println("Literate generation complete!")
println("Generated $(length(generated_pages)) example pages")

# Export the generated pages for use in make.jl
# This allows make.jl to dynamically build the pages array
const GENERATED_EXAMPLE_PAGES = generated_pages