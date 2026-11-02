#= This file defines the SLH type as well as fundamental functions for combining them =#

"""
    SLH(name, subspaces, inputs, outputs, S, L, H)

An SLH triple describes a quantum input-output system. See Combes, arXiv.1611.00375

The system acts on one [`Subspace`](@ref) per factor of its Hilbert space, and has `n` input
and `n` output ports, where signals enter and leave: S is n×n and L has n entries. The
operators in L and H are those of the subspaces, from [`operators`](@ref):

    subspaces = [OpticalMode("cavity"), MechanicalMode("mirror", PhaseSpace)]
    a, = operators(subspaces, "cavity")
    x, p = operators(subspaces, "mirror")

A system with no subspaces, such as a beamsplitter, is pure scattering: `SLH("bs", Subspace[],
["in1", "in2"], ["out1", "out2"], S, [0, 0], 0)`.

The names of the system, its subspaces and its ports cannot contain `_`, which joins them
into the names of everything in a composite system. The two ways of combining SLH systems
are [`concatenate`](@ref) and [`feedbackreduce`](@ref).
"""
struct SLH
    name::String
    subspaces::Vector{Subspace}
    inputs::Vector{String} #unique, length n
    outputs::Vector{String} #unique, length n
    S #size nxn, or a length 1 vector when n = 1
    L #size n
    H #has operators which act on hilbert

    # promoted = true for the systems built by concatenate and feedbackreduce, whose ports are already promoted
    function SLH(name, subspaces, inputs, outputs, S, L, H; promoted = false)

        check_name("system", name)
        promoted || foreach(port -> check_name("port", port), [inputs; outputs])

        #Coerce all operators to type QAdd, which can be iterated as (term,coeff) pairs
        L = [Lj isa SecondQuantizedAlgebra.Op ? 1 * Lj : Lj for Lj in L]
        H = H isa SecondQuantizedAlgebra.Op ? 1 * H : H

        # every operator acts on a declared subspace, belongs to its kind of Hilbert space, and is named after it
        m = length(subspaces)
        for op in unique(op for q in [L; H] for op in qops(q))
            k = only(acts_on(op))
            1 <= k <= m || error("$name: an operator acts on subspace $k, but only $m subspaces are declared")
            sub = subspaces[k]
            space = spacetype(op)
            sub isa Subspace{space} || error("$name: $op acts on subspace $(repr(sub.name)), which is not a $space")
            expected = subname(sub, operator_suffix(op))
            if operator_name(op) != expected
                return error("$name: $op acts on subspace $(repr(sub.name)), so it must be named $expected; get it from operators(subspaces, $(repr(sub.name)))")
            end
        end

        subnames = [sub.name for sub in subspaces]
        allunique(subnames) || return error("$name: subspace names must be unique, got $subnames")

        # one input, one output and one row and column of S per channel of L; a single channel may give S as [s]
        n = length(L)
        if !(length(inputs) == length(outputs) == n)
            return error("$name: L has $n channels, but there are $(length(inputs)) inputs and $(length(outputs)) outputs")
        end
        if !(size(S, 1) == size(S, 2) == n)
            return error("$name: S must be $n×$n to match the channels of L, got size $(size(S))")
        end

        # ports are connected and indexed by name
        allunique(inputs) || return error("$name: input names must be unique, got $inputs")
        allunique(outputs) || return error("$name: output names must be unique, got $outputs")

        return new(name, subspaces, inputs, outputs, S, L, H)
    end
end

# Names are joined by "_" into the names of everything in a composite system, "ifo_arm_cavity_a",
# so a name given by the user cannot contain "_"; otherwise the levels could not be told apart
function check_name(kind, name)
    occursin('_', name) || return nothing
    return error("the $kind name $(repr(name)) contains '_', which separates the levels of promoted names; use e.g. $(repr(replace(name, '_' => "")))")
end

# The operators in an entry of L or H, which may be a bare operator, a sum, or a c-number
qops(q::SecondQuantizedAlgebra.Op) = [q]
qops(q::SecondQuantizedAlgebra.QAdd) = [op for (term, _) in q for op in term.ops]
qops(q) = []

# The type of Hilbert space an operator acts on, and the suffix of its name, as in subspace_operators
function spacetype(op)
    (is_destroy(op) || is_create(op)) && return FockSpace
    (is_position(op) || is_momentum(op)) && return PhaseSpace
    is_transition(op) && return NLevelSpace
    return error("$op: no Subspace supports its Hilbert space yet")
end
function operator_suffix(op)
    (is_destroy(op) || is_create(op)) && return "a"
    is_position(op) && return "x"
    is_momentum(op) && return "p"
    is_transition(op) && return "σ"
    return error("$op: no Subspace supports its Hilbert space yet")
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
    concatenate(name, systems::SLH...) → SLH

Creates a composite system named `name` with no interconnections. Combes eq. 59

When systems are concatenated, the names of their inputs, outputs, operators,
parameters, and Hilbert spaces are 'promoted' by prepending the name of the
system to the existing name. This prevents name collisions as long as all
SLH systems are created with a unique name. A list of systems is splatted:
`concatenate(name, systems...)`.
"""
function concatenate(name, syslist::SLH...)

    # two systems of the same name would share their promoted ports, parameters and operators
    names = [sys.name for sys in syslist]
    allunique(names) || return error("$name: the systems being concatenated must have unique names, got $names")

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

    #we concatenate all the scattering matrices block diagonally, promoting their parameters as in L and H
    S = cat([[promote_expr(s, sys.name, 0) for s in sys.S] for sys in syslist]...; dims = (1, 2))

    L = vcat([[promote_expr(Lj, sys.name, off) for Lj in sys.L] for (sys, off) in zip(syslist, offsets)]...)

    H = sum([promote_expr(sys.H, sys.name, off) for (sys, off) in zip(syslist, offsets)])

    return SLH(name, newsubspaces, inputs, outputs, S, L, H; promoted = true)
end

# The indices of the output and input ports that feedbackreduce connects, for an SLH or a QuantumStateSpace
function connected_ports(sys, (output, input))
    x = findfirst(isequal(output), sys.outputs)
    y = findfirst(isequal(input), sys.inputs)
    if isnothing(x)
        swapped = output in sys.inputs ? " ($output is an input; the order is output => input)" : ""
        return error("$(sys.name): no output port named $(repr(output))$swapped. Outputs: $(sys.outputs)")
    end
    if isnothing(y)
        swapped = input in sys.outputs ? " ($input is an output; the order is output => input)" : ""
        return error("$(sys.name): no input port named $(repr(input))$swapped. Inputs: $(sys.inputs)")
    end
    return x, y
end

"""
    feedbackreduce(sys, output => input, ...)

Connects each output port to an input port, in order, reducing the number of outputs and inputs
by one each. Combes eq 61. For example, `feedbackreduce(sys, "A_out" => "B_in", "B_out" => "A_in")`.
Works on `SLH` and `QuantumStateSpace` systems.
"""
feedbackreduce(sys, connections::Pair...) = foldl(feedbackreduce, connections; init = sys)

function feedbackreduce(A::SLH, connection::Pair)

    output, input = connection
    x, y = connected_ports(A, connection)

    newoutputs = A.outputs[eachindex(A.outputs) .!= x]
    newinputs = A.inputs[eachindex(A.inputs) .!= y]

    Sxbarybar = A.S[1:end .!= x, 1:end .!= y]
    Sxbary = A.S[1:end .!= x, y]
    Sxybar = permutedims(A.S[x, 1:end .!= y])
    Sxy = A.S[x, y]
    # Combes eq 61 needs 1 - Sxy to be invertible; a symbolic Sxy that only simplifies to 1 is not caught
    if isequal(Sxy, 1)
        return error("$(A.name): connecting $output to $input closes a loop with unit gain (S[$x, $y] = 1), so 1 - S[$x, $y] is not invertible")
    end

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
    H = simplify(A.H + Hprime)

    return SLH(A.name, A.subspaces, newinputs, newoutputs, S, L, H; promoted = true)
end
