"""
    get_qnumbers(expr)

Extract all quantum operators from an expression.

# Arguments
- `expr`: A symbolic expression

# Returns
- `Set`: Set of quantum operators found in the expression

"""
function get_qnumbers(expr)
    if SymbolicUtils.iscall(expr)
        return union(get_qnumbers.(SymbolicUtils.arguments(expr))...)
    else
        return (expr isa SecondQuantizedAlgebra.QNumber) ? Set{SecondQuantizedAlgebra.QNumber}([expr]) : Set{SecondQuantizedAlgebra.QNumber}()
    end
end

"""
    get_numsymbols(expr)

Extract all symbolic parameters from an expression.

# Arguments
- `expr`: A symbolic expression

# Returns
- `Set`: Set of symbolic parameters found in the expression

"""
function get_cnumbers(expr)
    if SymbolicUtils.iscall(expr)
        return union(get_cnumbers.(SymbolicUtils.arguments(expr))...)
    else
        return (expr isa SymbolicUtils.Symbolic) ? Set{SymbolicUtils.Symbolic}([expr]) : Set{SymbolicUtils.Symbolic}() 
    end
end

"""
    get_additive_terms(expr)

Extract additive terms from a quantum operator expression.

Takes an expression containing quantum operators and returns a list of terms
that contain no addition, only multiplication. Summing all returned terms
results in the original expression.

# Arguments
- `expr`: A symbolic expression containing quantum operators

# Returns
- `Vector`: List of terms without addition operators

"""
function get_additive_terms(expr)
    if SymbolicUtils.iscall(expr)
        op = SymbolicUtils.operation(expr)
        args = SymbolicUtils.arguments(expr)
        
        if op === (+)
            # If this is an addition, recursively get terms from each argument
            terms = []
            for arg in args
                append!(terms, get_additive_terms(arg))
            end
            return terms
        else
            # If this is not an addition (multiplication, function call, etc.), 
            # return the whole expression as a single term
            return [expr]
        end
    else
        # If not a tree (atom), return as single term
        return [expr]
    end
end

function ordered_qsymbols(expr)
    if SymbolicUtils.iscall(expr)
        return cat(ordered_qsymbols.(SymbolicUtils.arguments(expr))...,dims=1)
    else
        return (expr isa SecondQuantizedAlgebra.QNumber) ? SecondQuantizedAlgebra.QNumber[expr] : SecondQuantizedAlgebra.QNumber[]
    end
end


function islinear(expr)
    terms = get_additive_terms(expr)
    for term in terms 
        if length(ordered_qsymbols(term)) > 2
            return false
        end
    end
    return true
end


