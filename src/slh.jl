#= This file defines the SLH type as well as fundamental functions for combining them =#

"""
SLH(name, inputs, outputs, S, L, H)

An SLH triple describes an open quantum system. See Combes, arXiv.1611.00375

The name of the system should be unique. When multiple systems are combined, the names of their inputs and outputs will 
have the system name appended to them. The inputs and outputs describe 'ports' where signals leave and enter the system.
Quantum systems must have the same number of inputs and outputs, which we denote by n.

size(S) = (n, n) <- S is an nxn matrix

size(L) = (n,)

size(H) = ()

The two ways of combining SLH systems are concatenate() and feedbackreduce()
"""
struct SLH
    name::String
    subspaces::Vector{Subspace}
    inputs::Vector{String} #must have unique elements
    outputs::Vector{String} #must have unique elements
    S #size nxn
    L #size n
    H #has operators which act on hilbert

    #Coerce all operators to type QAdd, which can be iterated as (term,coeff) pairs
    function SLH(name, subspaces, inputs, outputs, S, L, H)
        L = [Lj isa SecondQuantizedAlgebra.Op ? 1 * Lj : Lj for Lj in L]
        H = H isa SecondQuantizedAlgebra.Op ? 1 * H : H

        # every operator must act on a declared subspace
        m = length(subspaces)
        for q in [L; H]
            q isa SecondQuantizedAlgebra.QAdd || continue
            for k in acts_on(q)
                if !(1 <= k <= m)
                    return error("$name: an operator acts on hilbert space $k, but only $m subspaces are declared")
                end
            end
        end

        return new(name, subspaces, inputs, outputs, S, L, H)
    end
end

#This constructor automatically makes all the modes 'generic' and names the inputs and outputs in1, in2,... out1, out2...
function SLH(name, S, L, H)

    n = maximum(acts_on(sum(L) + H))

    subspaces = fill(GenericMode(""), n)

    m = length(L)
    if m == 1
        inputs = ["in"]
        outputs = ["out"]
    else
        inputs = ["in$j" for j in 1:m]
        outputs = ["out$j" for j in 1:m]
    end

    return SLH(name, subspaces, inputs, outputs, S, L, H)

end

#This function is for SecondQuantizedAlgebra operators
function promote_op(op, aon_offset, topname)

    new_op_name = Symbol(topname, "_", operator_name(op))

    subspaceindex = aon_offset + only(acts_on(op))

    return SecondQuantizedAlgebra.set_acts_on(SecondQuantizedAlgebra.rename(op, new_op_name), subspaceindex)
end

#SQA's substitute applies all rules simultaneously and once, so shifting modes 1→2 and 2→3
#in the same expression cannot chain, and conj(η) is promoted along with η
function promote_expr(q, topname, offset)
    scalar_rules = Dict{Any, Any}(v => promote_name(v, topname) for v in Symbolics.get_variables(q))
    if !(q isa SecondQuantizedAlgebra.QAdd)
        return substitute(q, scalar_rules)
    end
    ops = unique(op for (term, _) in q for op in term.ops)
    op_rules = Dict{Any, Any}(op => promote_op(op, offset, topname) for op in ops)
    return substitute(q, merge(op_rules, scalar_rules))
end

"""
concatenate(name, syslist::Vector{SLH})

creates a composite system with no interconnections. Combes eq. 59

When systems are concatenated, the names of their inputs, outputs, operators,
parameters, and Hilbert spaces are 'promoted' by prepending the name of the 
system to the existing name. This prevents name collisions as long as all 
SLHSystems are created with a unique name.
"""
function concatenate(syslist, name) #Should this be changed to concatenate(name, systems...) ?

    ####
    # Subspaces
    newsubspaces = vcat([[promote_name(subsys, sys.name) for subsys in sys.subspaces] for sys in syslist]...)
    ###

    #We 'stack' the inputs and outputs of the systems we are combining.
    #first, we promote the names of inputs and outputs, to prevent naming collisions
    newinputs = [[sys.name * "_" * input for input in sys.inputs] for sys in syslist]
    inputs = cat(newinputs..., dims = 1)

    newoutputs = [[sys.name * "_" * output for output in sys.outputs] for sys in syslist]
    outputs = cat(newoutputs..., dims = 1)

    offsets = cumsum([0; [length(sys.subspaces) for sys in syslist[1:(end - 1)]]])

    #we concatenate all the scattering matrices block diagonally
    S = cat([sys.S for sys in syslist]...; dims = (1, 2))

    L = vcat([[promote_expr(Lj, sys.name, off) for Lj in sys.L] for (sys, off) in zip(syslist, offsets)]...)

    H = sum([promote_expr(sys.H, sys.name, off) for (sys, off) in zip(syslist, offsets)])

    return SLH(name, newsubspaces, inputs, outputs, S, L, H)
end

"""
feedbackreduce(A::SLH,output,input)

Connects the output port to the input port, reducing the number of outputs and inputs by one each. Combes eq 61.
"""
function feedbackreduce(A, output, input)

    x = findfirst(isequal(output), A.outputs)
    y = findfirst(isequal(input), A.inputs)

    newoutputs = A.outputs[eachindex(A.outputs) .!= x]
    newinputs = A.inputs[eachindex(A.inputs) .!= y]

    Sxbarybar = A.S[1:end .!= x, 1:end .!= y]
    Sxbary = A.S[1:end .!= x, y]
    Sxybar = permutedims(A.S[x, 1:end .!= y])
    Sxy = A.S[x, y]

    Sy = A.S[:, y]

    Lxbar = A.L[1:end .!= x]
    Lx = A.L[x]

    S = Sxbarybar + Sxbary * (1 - Sxy)^(-1) * Sxybar
    #Have to use fill here because the usual broadcast syntax on operator does not work
    L = Lxbar + Sxbary .* fill(((1 - Sxy)^(-1) * Lx), size(Sxbary))

    term1 = adjoint(A.L) * Sy
    term2 = (1 - Sxy)^(-1)
    term3 = Lx
    termA = adjoint(Lx)
    termB = ((1 - adjoint(Sxy))^(-1))
    termC = (adjoint(Sy) * A.L)

    Hprime = 1 / (2im) * (term1 * term2 * term3 - termA * termB * termC)

    #return (A.H, Hprime)
    H = simplify(A.H + Hprime)

    return SLH(A.name, A.subspaces, newinputs, newoutputs, S, L, H)
end
