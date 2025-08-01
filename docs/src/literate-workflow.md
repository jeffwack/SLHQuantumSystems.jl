# Literate.jl Workflow

This documentation uses [Literate.jl](https://fredrikekre.github.io/Literate.jl/) to automatically generate documentation pages from executable Julia scripts. This approach ensures that a single file can be executed directly while also serving as formatted documentation with explanations.

## How It Works - Script Discovery and Generation

When the documentation is built (`julia docs/make.jl`), the following process occurs:

1. **Auto-discovery**: The build system scans `examples/` for all `.jl` files
2. **Conversion**: Each script is processed by Literate.jl to generate markdown in docs/src/generated/
3. **Integration**: Generated markdown files are automatically added to the documentation by adding pages in make.jl

## Writing Literate Scripts

Literate.jl uses a simple convention for mixing code and documentation:

- **Documentation**: Lines starting with `#` (with a space) become markdown
- **Code**: Regular Julia code becomes executable code blocks
- **Headers**: Use `# # Header` for markdown headers (note the extra `#`)

See [the Literate.jl docs](https://fredrikekre.github.io/Literate.jl/v2/) for more information.

## File Organization

```
SLHQuantumSystems.jl/
├── examples/
│   ├── cascadedcavities.jl     # Literate script
│   └── newexample.jl           # Automatically discovered
├── docs/
│   ├── src/
│   │   ├── generated/          # Auto-generated (don't edit!)
│   │   │   ├── cascadedcavities.md
│   │   │   └── newexample.md
│   │   └── index.md
│   ├── generate.jl             # Generation script
│   └── make.jl                 # Documentation build
```

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

