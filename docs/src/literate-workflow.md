# Literate.jl Workflow

This documentation uses [Literate.jl](https://fredrikekre.github.io/Literate.jl/) to automatically generate documentation pages from executable Julia scripts. This approach ensures that a single file can be executed directly while also serving as formatted documentation with explanations.

## How It Works - Script Discovery and Generation

When the documentation is built, the following process occurs:

1. **Auto-discovery**: The build system scans `examples/` for all `.jl` files using `walkdir()`
2. **Conversion**: Each script is processed by Literate.jl to generate markdown in `docs/src/generated/`
3. **Integration**: Generated markdown files are automatically added to the documentation via `GENERATED_EXAMPLE_PAGES`

## Development with LiveServer

For live development, use LiveServer with proper configuration to avoid double triggers:

```julia
using LiveServer
servedocs(literate_dir="examples", skip_dir="docs/src/generated")
```

This setup:
- Watches source files for changes and rebuilds docs automatically
- Uses `literate_dir="examples"` to tell LiveServer where the Literate scripts are
- Uses `skip_dir="docs/src/generated"` to prevent watching generated files (avoids infinite loops)

## Writing Literate Scripts

Literate.jl uses a simple convention for mixing code and documentation:

- **Documentation**: Lines starting with `#` (with a space) become markdown
- **Code**: Regular Julia code becomes executable code blocks
- **Headers**: Use `# # Header` for markdown headers (note the extra `#`)

See [the Literate.jl docs](https://fredrikekre.github.io/Literate.jl/v2/) for more information.

## File Organization

```
SLHQuantumSystems.jl/
├── examples/                    # Literate scripts (input)
│   ├── cascadedcavities.jl     
│   ├── cavityfeedback.jl       
│   └── fabry_perot.jl          # Automatically discovered
├── docs/
│   ├── src/
│   │   ├── generated/          # Auto-generated markdown (don't edit!)
│   │   │   ├── cascadedcavities.md
│   │   │   ├── cavityfeedback.md
│   │   │   └── fabry_perot.md
│   │   ├── index.md
│   │   └── api.md
│   ├── generate.jl             # Literate processing script
│   └── make.jl                 # Documenter build script
```

## Build Process

The workflow follows the official LiveServer + Literate.jl pattern:

1. **generate.jl**: Uses `walkdir()` to discover all `.jl` files in `examples/`
2. **Literate.markdown()**: Converts each script to markdown with `documenter=true, execute=false`
3. **make.jl**: Includes `generate.jl`, then runs `makedocs()` with dynamically generated pages
4. **LiveServer**: Watches files but skips the generated directory to prevent double builds

## Adding New Examples

To add a new example to the documentation:

1. **Create the script**: Add a new `.jl` file to the `examples/` directory
2. **Use Literate format**: Write using the comment conventions above

No manual configuration is needed! The system automatically:
- Discovers the new file
- Converts it to markdown
- Adds it to the documentation navigation
- Creates a properly formatted title from the filename

## Manual Testing

You can test individual examples by running them directly:

```julia
julia> include("examples/cascadedcavities.jl")
```

This ensures that your examples actually work before they become documentation.

