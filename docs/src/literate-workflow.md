# Literate.jl Workflow

This documentation uses [Literate.jl](https://fredrikekre.github.io/Literate.jl/) to automatically generate documentation pages from executable Julia scripts. This approach ensures that examples in the documentation are always runnable and stay in sync with the codebase.

## Overview

The workflow automatically converts Julia scripts in the `examples/` directory into documentation pages. Each script can be executed directly while also serving as formatted documentation with explanations.

## How It Works

### 1. Script Discovery and Generation

When the documentation is built (`julia docs/make.jl`), the following process occurs:

1. **Auto-discovery**: The build system scans `examples/` for all `.jl` files
2. **Conversion**: Each script is processed by Literate.jl to generate markdown
3. **Integration**: Generated markdown files are automatically added to the documentation navigation

### 2. Build Process Flow

```
examples/*.jl → Literate.jl → docs/src/generated/*.md → Documenter.jl → HTML docs
```

The key insight is that markdown generation happens **before** Documenter runs, ensuring all files exist when the documentation is built.

## Writing Literate Scripts

### Comment Format

Literate.jl uses a simple convention for mixing code and documentation:

- **Documentation**: Lines starting with `#` (with a space) become markdown
- **Code**: Regular Julia code becomes executable code blocks
- **Headers**: Use `# # Header` for markdown headers (note the extra `#`)

### Example Structure

```julia
# # Example Title
#
# This is a documentation paragraph explaining what the example does.
# You can use **markdown formatting** and [links](https://example.com).
#
# ## Section Header
#
# More explanation here.

# This is regular Julia code
using SomePackage

# More documentation
# explaining the next code block

result = some_computation()
println(result)
```

### Best Practices

1. **Start with a descriptive title**: Use `# # Title` as the first line
2. **Explain before showing**: Add explanatory text before code blocks
3. **Use sections**: Break long examples into logical sections
4. **Add context**: Explain why each step is necessary
5. **Show results**: Include `println()` or display statements to show output

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
3. **Build docs**: Run `julia docs/make.jl` - the new example will automatically appear

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

## Generated Content

The build process creates:

- **Markdown files**: In `docs/src/generated/` (auto-generated, don't edit)
- **Navigation entries**: Automatically added to the documentation menu
- **Cross-references**: Proper links between documentation pages

### Generated Features

Literate.jl with `documenter=true` provides:

- **Proper code blocks**: Converts code to `@example` blocks for execution
- **Edit links**: "Edit on GitHub" links point to the source `.jl` file
- **Clean formatting**: Removes Documenter-specific syntax for other outputs

## Troubleshooting

### Common Issues

1. **Build failures**: Check that your script runs without errors when included directly
2. **Missing pages**: Ensure your script is in the `examples/` directory and ends with `.jl`
3. **Formatting issues**: Verify you're using `# ` (hash + space) for documentation lines

### Debugging

- **Test scripts directly**: `julia> include("examples/yourscript.jl")`
- **Check generated files**: Look in `docs/src/generated/` after building
- **Examine build output**: The generation process prints status messages

## Technical Details

### Generation Script Location

The generation logic is in `docs/generate.jl`, which:
- Scans the examples directory
- Processes each file with `Literate.markdown()`
- Creates the `GENERATED_EXAMPLE_PAGES` array for navigation

### Integration with Documenter

The `docs/make.jl` file:
1. Imports Literate.jl
2. Runs the generation script
3. Uses the generated page list in `makedocs()`

This ensures the markdown files exist before Documenter tries to process them.

---

*This workflow provides a maintainable way to keep documentation examples current and executable, reducing the burden of maintaining separate documentation and example code.*