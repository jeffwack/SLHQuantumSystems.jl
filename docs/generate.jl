using Literate

# Define paths
EXAMPLEDIR = joinpath(@__DIR__, "..", "examples")
BUILDDIR = joinpath(@__DIR__, "src", "generated")

# Ensure the build directory exists
mkpath(BUILDDIR)

# Check if we should skip generation (for LiveServer compatibility)
# Skip generation if output files exist and are newer than source files
function should_skip_generation()
    example_files = filter(f -> endswith(f, ".jl"), readdir(EXAMPLEDIR))
    
    for example in example_files
        example_path = joinpath(EXAMPLEDIR, example)
        base_name = splitext(example)[1]
        output_path = joinpath(BUILDDIR, "$(base_name).md")
        
        # If output doesn't exist or is older than source, we need to generate
        if !isfile(output_path) || stat(example_path).mtime > stat(output_path).mtime
            return false
        end
    end
    
    return true
end

# Skip generation if files are up-to-date (prevents LiveServer infinite loop)
if should_skip_generation()
    println("Literate files are up-to-date, skipping generation...")
else
    println("Generating documentation from Literate examples...")
end

# Automatically discover all .jl files in examples directory
example_files = filter(f -> endswith(f, ".jl"), readdir(EXAMPLEDIR))
generated_pages = []

# Process each example file (only if needed)
for example in example_files
    example_path = joinpath(EXAMPLEDIR, example)
    base_name = splitext(example)[1]
    output_path = joinpath(BUILDDIR, "$(base_name).md")
    
    if isfile(example_path)
        # Only regenerate if source is newer than output or output doesn't exist
        needs_generation = !isfile(output_path) || stat(example_path).mtime > stat(output_path).mtime
        
        if needs_generation && !should_skip_generation()
            println("Processing: $example")
            
            # Generate markdown for documentation
            Literate.markdown(
                example_path, 
                BUILDDIR; 
                documenter=true,  # Enable Documenter.jl integration
                execute=false     # Don't execute during docs build for reliability
            )
            
            println("  → Generated markdown: $(base_name).md")
        end
        
        # Always add to pages list (even if not regenerated)
        title = titlecase(replace(base_name, r"([a-z])([A-Z])" => s"\1 \2"))
        page_path = "generated/$(base_name).md"
        push!(generated_pages, title => page_path)
        
    else
        @warn "Example file not found: $example_path"
    end
end

println("Literate generation complete!")
println("Generated $(length(generated_pages)) example pages")

# Export the generated pages for use in make.jl
# This allows make.jl to dynamically build the pages array
const GENERATED_EXAMPLE_PAGES = generated_pages
