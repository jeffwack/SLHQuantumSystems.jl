#= This file defines a type for the representation of linear systems
=#

using ControlSystemsBase: AbstractStateSpace, Continuous
import ControlSystemsBase: system_name, input_names, output_names, state_names,
    series, parallel, append, feedback, StateSpace

abstract type BasisType end

"""
Annihilation-creation operator basis, the "complex mode form."
State vector ordered as (a₁, a₁†, a₂, a₂†, ...)
ABCD matrices are complex-valued.
"""
struct AnnihilationCreationBasis <: BasisType end

"""
Quadrature (amplitude/phase) operator basis.
State vector ordered as (x₁, p₁, x₂, p₂, ...).
A mode on a `FockSpace` has the dimensionless quadratures
x = (a+a†)/√2, p = i(a†-a)/√2, whatever its type.
A mode on a `PhaseSpace` has states x and p, in the units of H and its parameters.
ABCD matrices are real-valued for physical systems.
"""
struct QuadratureBasis <: BasisType end

struct QuantumStateSpace{TE, TB <: BasisType} <: AbstractStateSpace{TE}
    # SLH metadata
    name::String
    subspaces::Vector
    inputs::Vector{String}
    outputs::Vector{String}
    # ABCD matrices (Any-typed to support both Num and numeric)
    A
    B
    C
    D
    # Required by AbstractStateSpace interface
    timeevol::TE
    basis::TB
end

system_name(sys::QuantumStateSpace) = sys.name

# ── Names ─────────────────────────────────────────────────────────────────────
# Ports are named by Strings and states and their associated operators are named by Symbols.
# Each mode or port has two entries in the state or I/O vector, named by appending the basis suffixes.

suffixes(::AnnihilationCreationBasis) = ("a", "adag")
suffixes(::QuadratureBasis) = ("x", "p")

"""
    subspace_names(sys) → Vector{String}

The name of each subspace, in order. Empty for an unnamed subspace, which takes the
system's name on `concatenate`. In a `QuantumStateSpace` every subspace is a mode.
"""
subspace_names(sys::Union{SLH, QuantumStateSpace}) = [sub.name for sub in sys.subspaces]

"""
    subspace(sys, name) → Subspace

The subspace of `sys` named `name`.
"""
subspace(sys::Union{SLH, QuantumStateSpace}, name) = sys.subspaces[subspaceindex(sys.subspaces, name)]

function subspaceindex(subs, name)
    k = findfirst(sub -> sub.name == name, subs)
    isnothing(k) && error("no subspace named $(repr(name)). Subspaces: $([sub.name for sub in subs])")
    return k
end

"""
    hilbert(sys) → HilbertSpace
    hilbert(subspaces) → HilbertSpace

The Hilbert space of a system or of a list of subspaces: the space of a single subspace, or
the product of their spaces in order.
"""
hilbert(sys::Union{SLH, QuantumStateSpace}) = hilbert(sys.subspaces)
hilbert(subs::AbstractVector{<:Subspace}) = length(subs) == 1 ? hilbertspace(only(subs)) : ⊗(hilbertspace.(subs)...)

"""
    operators(sys, name) → Vector
    operators(subspaces, name) → Vector
    operators(sys) → Vector

The operators that act on the subspace named `name`, by the type of its Hilbert space and named
after the subspace: `[cavity_a]` for a `FockSpace`, `[mirror_x, mirror_p]` for a `PhaseSpace`,
and the matrix of transitions `atom_σᵢⱼ` for an `NLevelSpace`. Write H and L with them:

    subs = [OpticalMode("cavity"), MechanicalMode("mirror", PhaseSpace)]
    a, = operators(subs, "cavity")
    x, p = operators(subs, "mirror")

Without a name, the operators of every subspace.
"""
operators(sys::Union{SLH, QuantumStateSpace}, name) = operators(sys.subspaces, name)
function operators(subs::AbstractVector{<:Subspace}, name)
    k = subspaceindex(subs, name)
    at = length(subs) == 1 ? () : (k,)
    return subspace_operators(subs[k], hilbert(subs), at...)
end
operators(sys::Union{SLH, QuantumStateSpace}) = vcat((vec(operators(sys, sub.name)) for sub in sys.subspaces)...)

"""
    parameters(sys, name) → Vector
    parameters(subspaces, name) → Vector
    parameters(sys) → Vector

The parameters of the subspace named `name`, by its type and named after the subspace:
`[cavity_L]` for an `OpticalMode`, `[mirror_m]` for a `MechanicalMode`.

Without a name, every parameter of `sys`: the symbols in its S, L and H, or its ABCD
matrices, and the parameters of its subspaces.
"""
parameters(sys::Union{SLH, QuantumStateSpace}, name) = parameters(sys.subspaces, name)
parameters(subs::AbstractVector{<:Subspace}, name) = subspace_parameters(subs[subspaceindex(subs, name)])
parameters(sys::SLH) = symbols([vec(collect(sys.S)); sys.L; sys.H; subspace_parameters.(sys.subspaces)...])
parameters(sys::QuantumStateSpace) = symbols([vec(sys.A); vec(sys.B); vec(sys.C); vec(sys.D); subspace_parameters.(sys.subspaces)...])

symbols(xs) = sort(unique(Symbolics.wrap(v) for x in xs for v in Symbolics.get_variables(flattencomplex(x))); by = string)

"""
    state_names(sys::QuantumStateSpace) → Vector{Symbol}

The name of each entry of the state vector: `[:arm_x, :arm_p, …]` in the quadrature basis
and `[:arm_a, :arm_adag, …]` in the annihilation-creation basis. An unnamed mode gives `[:x, :p]`.
"""
state_names(sys::QuantumStateSpace) = [subname(sub, s) for sub in sys.subspaces for s in suffixes(sys.basis)]
state_names(sys::QuantumStateSpace, i) = state_names(sys)[i]

"""
    input_ports(sys) → Vector{String}

The name of each input port, e.g. `["in"]`. Each port has two entries in the input
vector; see [`input_names`](@ref).
"""
input_ports(sys::Union{SLH, QuantumStateSpace}) = sys.inputs

"""
    output_ports(sys) → Vector{String}

The name of each output port, e.g. `["out"]`. Each port has two entries in the output
vector; see [`output_names`](@ref).
"""
output_ports(sys::Union{SLH, QuantumStateSpace}) = sys.outputs

"""
    input_names(sys::QuantumStateSpace) → Vector{String}
    input_names(sys::QuantumStateSpace, i)

The name of each entry of the input vector, labelling the columns of D: `["in_x", "in_p", …]`
in the quadrature basis and `["in_a", "in_adag", …]` in the annihilation-creation basis. The port names
alone are [`input_ports`](@ref).
"""
input_names(sys::QuantumStateSpace) = expand_ports(input_ports(sys), sys.basis)

"""
    output_names(sys::QuantumStateSpace) → Vector{String}
    output_names(sys::QuantumStateSpace, i)

The name of each entry of the output vector, labelling the rows of D: `["out_x", "out_p", …]`
in the quadrature basis and `["out_a", "out_adag", …]` in the annihilation-creation basis. The port names
alone are [`output_ports`](@ref).
"""
output_names(sys::QuantumStateSpace) = expand_ports(output_ports(sys), sys.basis)
input_names(sys::QuantumStateSpace, i) = input_names(sys)[i]
output_names(sys::QuantumStateSpace, i) = output_names(sys)[i]

expand_ports(ports, basis) = [port * "_" * s for port in ports for s in suffixes(basis)]

# A plain StateSpace, without names or basis, for the ControlSystemsBase functions that
# build a new system (series, append, step, …). ss(sys) calls this.
StateSpace(sys::QuantumStateSpace) = StateSpace(sys.A, sys.B, sys.C, sys.D, sys.timeevol)
Base.convert(::Type{StateSpace}, sys::QuantumStateSpace) = StateSpace(sys)

# ── Composition ───────────────────────────────────────────────────────────────
# concatenate and feedbackreduce follow their SLH counterparts, including name promotion, so that
# building the state space commutes with composing: QuantumStateSpace(concatenate(name, slhs...)) is
# concatenate(name, QuantumStateSpace.(slhs)...), and likewise for feedbackreduce.

"""
    concatenate(name, systems::QuantumStateSpace...) → QuantumStateSpace

Place state spaces side by side with no interconnections, with block-diagonal ABCD matrices.
As for SLH systems, the modes, ports and parameters of each system are promoted by prepending
its name, so building the state space before or after concatenating gives the same system.
The systems must share a basis.
"""
function concatenate(name, syslist::QuantumStateSpace{TE, TB}...) where {TE, TB}
    check_name("system", name)
    names = [sys.name for sys in syslist]
    allunique(names) || return error("$name: the systems being concatenated must have unique names, got $names")
    basis = first(syslist).basis

    subspaces = vcat([[promote_name(mode, sys.name) for mode in sys.subspaces] for sys in syslist]...)
    inputs = vcat([[sys.name * "_" * input for input in sys.inputs] for sys in syslist]...)
    outputs = vcat([[sys.name * "_" * output for output in sys.outputs] for sys in syslist]...)

    # the SLH constructor's checks; distinct system names can still promote to the same name
    subnames = [sub.name for sub in subspaces]
    allunique(subnames) || return error("$name: subspace names must be unique, got $subnames")
    allunique(inputs) || return error("$name: input names must be unique, got $inputs")
    allunique(outputs) || return error("$name: output names must be unique, got $outputs")

    promoted = [promote_matrices(sys) for sys in syslist]
    A, B, C, D = (blockdiag([mats[k] for mats in promoted]) for k in 1:4)

    return QuantumStateSpace(name, subspaces, inputs, outputs, A, B, C, D, first(syslist).timeevol, basis)
end

# The methods above take systems of one kind, and state spaces in one basis
function concatenate(name, systems...)
    if all(sys -> sys isa QuantumStateSpace, systems)
        return error("$name: the systems being concatenated must share a basis; convert with QuantumStateSpace(sys, QuadratureBasis())")
    end
    kinds = join(unique(nameof.(typeof.(systems))), " and ")
    return error("$name: cannot concatenate $kinds; build each SLH system's QuantumStateSpace first")
end
concatenate(name) = error("$name: concatenate needs at least one system")

# cat(Ms...; dims = (1, 2)) fills the off-diagonal blocks with zero(eltype), which fails when the
# entries of different systems promote to Any
function blockdiag(Ms)
    T = mapreduce(eltype, promote_type, Ms)
    rows, cols = sum(M -> size(M, 1), Ms), sum(M -> size(M, 2), Ms)
    out = T <: Union{Float64, ComplexF64} ? zeros(T, rows, cols) : Array{Any}(zeros(rows, cols))
    i, j = 0, 0
    for M in Ms
        out[i .+ axes(M, 1), j .+ axes(M, 2)] = M
        i, j = i + size(M, 1), j + size(M, 2)
    end
    return out
end

# The ABCD matrices with each parameter renamed by promote_name, as concatenate does for S, L and H
function promote_matrices(sys::QuantumStateSpace)
    mats = (sys.A, sys.B, sys.C, sys.D)
    params = unique(v for M in mats for x in M for v in Symbolics.get_variables(flattencomplex(x)))
    isempty(params) && return mats
    rules = Dict{Any, Any}(v => promote_name(v, sys.name) for v in params)
    return Tuple(map(x -> Symbolics.wrap(substitute(x, rules)), M) for M in mats)
end

# Connect the output port to the input port, as for SLH systems, which removes both ports. Their
# two entries are joined, `u_input = y_output`, and the loop is solved:
#
#     K  = (I - D_xy)⁻¹
#     A' = A + B_y K C_x                 B' = B_ȳ + B_y K D_xȳ
#     C' = C_x̄ + D_x̄y K C_x              D' = D_x̄ȳ + D_x̄y K D_xȳ
#
# where x and y are the entries of the connected output and input, and x̄ and ȳ the rest.
function feedbackreduce(sys::QuantumStateSpace, connection::Pair)
    output, input = connection
    x, y = connected_ports(sys, connection)
    rx, cy = (2x - 1):2x, (2y - 1):2y
    r̄, c̄ = setdiff(axes(sys.D, 1), rx), setdiff(axes(sys.D, 2), cy)

    # flattencomplex makes each entry a single expression, since dividing a Complex{Num} compares its parts
    A, B, C, D = (map(flattencomplex, M) for M in (sys.A, sys.B, sys.C, sys.D))

    # D_xy is S_xy in the basis of the port, and det(I - D_xy) = |1 - S_xy|²
    M = I - D[rx, cy]
    d = Symbolics.wrap(expandsimplify(M[1, 1] * M[2, 2] - M[1, 2] * M[2, 1]))
    if issymzero(d)
        return error("$(sys.name): connecting $output to $input closes a loop with unit gain, so I - D between them is not invertible")
    end
    # elementwise, since a matrix divided by a bare symbolic value becomes a single symbolic array expression
    K = [M[2, 2] -M[1, 2]; -M[2, 1] M[1, 1]] ./ d

    newA = expandsimplify.(A + B[:, cy] * K * C[rx, :])
    newB = expandsimplify.(B[:, c̄] + B[:, cy] * K * D[rx, c̄])
    newC = expandsimplify.(C[r̄, :] + D[r̄, cy] * K * C[rx, :])
    newD = expandsimplify.(D[r̄, c̄] + D[r̄, cy] * K * D[rx, c̄])

    inputs = sys.inputs[eachindex(sys.inputs) .!= y]
    outputs = sys.outputs[eachindex(sys.outputs) .!= x]
    return QuantumStateSpace(sys.name, sys.subspaces, inputs, outputs, newA, newB, newC, newD, sys.timeevol, sys.basis)
end

issymzero(x) = (v = Symbolics.value(x); v isa Number && iszero(v))

# ControlSystemsBase's compositions connect every port at once and drop the names, and parallel
# connection copies a quantum field, so they are not defined for a QuantumStateSpace
const _COMPOSE_ERROR = "compose a QuantumStateSpace with concatenate and feedbackreduce, which keep its names, or use ss(sys) for ControlSystemsBase's compositions"

Base.:*(::QuantumStateSpace, ::QuantumStateSpace) = error(_COMPOSE_ERROR)
Base.:+(::QuantumStateSpace, ::QuantumStateSpace) = error(_COMPOSE_ERROR)
Base.:-(::QuantumStateSpace, ::QuantumStateSpace) = error(_COMPOSE_ERROR)
series(::QuantumStateSpace, ::QuantumStateSpace) = error(_COMPOSE_ERROR)
parallel(::QuantumStateSpace, ::QuantumStateSpace) = error(_COMPOSE_ERROR)
append(::QuantumStateSpace, ::QuantumStateSpace...) = error(_COMPOSE_ERROR)
feedback(::QuantumStateSpace, ::QuantumStateSpace; kwargs...) = error(_COMPOSE_ERROR)
feedback(::QuantumStateSpace; kwargs...) = error(_COMPOSE_ERROR)

# Default: build in the quadrature basis, which accepts both ladder and position/momentum operators
function QuantumStateSpace(sys::SLH)
    return _build_quadrature_ss(sys)
end

function QuantumStateSpace(sys::SLH, ::AnnihilationCreationBasis)
    return _build_annihilation_creation_ss(sys)
end

function QuantumStateSpace(sys::SLH, ::QuadratureBasis)
    return QuantumStateSpace(sys)
end

"""
    _build_annihilation_creation_ss(sys::SLH) → QuantumStateSpace{Continuous, AnnihilationCreationBasis}

Build the annihilation-creation-basis ABCD matrices of a linear SLH system, following Combes eqs. 116-119.
H and L are parsed term by term into

    H = a†Ω₋a + ½(a†Ω₊a†ᵀ + aᵀΞa)        L = Φ₋a + Φ₊a†ᵀ

Combes writes the a† rows of A as the conjugates of the a rows, which holds only for a
Hermitian H. Here they are built from i[H, a†] directly, so H need not be Hermitian.
"""
function _build_annihilation_creation_ss(sys::SLH)

    S = sys.S
    L = sys.L
    H = sys.H

    m = length(sys.subspaces)

    # We separate our terms into energy conserving (omegaminus) and non-conserving (omegaplus for a†a†, Ξ for aa).
    omegaminus = Array{Any}(zeros(m, m))
    omegaplus = Array{Any}(zeros(m, m))
    Ξ = Array{Any}(zeros(m, m))

    # a c-number H (no operators at all, e.g. 0) is not a QAdd, and only shifts the energy
    for (term, coeff) in (H isa SecondQuantizedAlgebra.QAdd ? H : ())

        ops = term.ops

        if isempty(ops)
            #cnumbers do not contribute to Heisenberg eqs
            continue
        end

        if length(ops) != 2
            return error("$ops term of the Hamiltonian is not quadratic")
        end

        for op in ops
            ladder_only(op, "the Hamiltonian")
        end

        c = SecondQuantizedAlgebra.to_num(coeff)
        creators = filter(is_create, ops)

        if length(creators) == 1
            j = only(acts_on(only(creators)))
            k = only(acts_on(only(filter(is_destroy, ops))))
            omegaminus[j, k] += c
        elseif length(creators) == 2
            j = only(acts_on(creators[1]))
            k = only(acts_on(creators[2]))
            omegaplus[j, k] += c
            omegaplus[k, j] += c
        elseif length(creators) == 0
            j = only(acts_on(ops[1]))
            k = only(acts_on(ops[2]))
            Ξ[j, k] += c
            Ξ[k, j] += c
        else
            return error("don't know how to handle $term")
        end
    end

    #Now we parse L
    n = length(L)
    phiminus = Array{Any}(zeros(n, m))
    phiplus = Array{Any}(zeros(n, m))

    for (j, channel) in enumerate(L)
        # a c-number L[j] (no operators at all) is not a QAdd
        if !(channel isa SecondQuantizedAlgebra.QAdd)
            isequal(channel, 0) && continue
            return error("L[$j] = $channel is a constant (a coherent drive), which the linear ABCD model cannot represent")
        end
        for (term, coeff) in channel
            ops = term.ops
            if isempty(ops)
                return error("L[$j] contains a constant term $(SecondQuantizedAlgebra.to_num(coeff)) (a coherent drive), which the linear ABCD model cannot represent")
            end
            if length(ops) != 1
                return error("$ops term of L[$j] is not linear")
            end
            op = only(ops)
            ladder_only(op, "L[$j]")
            c = SecondQuantizedAlgebra.to_num(coeff)
            k = only(acts_on(op))
            if is_destroy(op)
                phiminus[j, k] += c
            else
                phiplus[j, k] += c
            end
        end
    end

    #Now we need to construct the matrices according to Combes eqs 116-119

    Phi = [
        phiminus phiplus;
        conj.(phiplus) conj.(phiminus)
    ]

    # i[H, a] = -i(Ω₋a + Ω₊a†) and i[H, a†] = -i(-Ξa - Ω₋ᵀa†). For a Hermitian H, Ξ = Ω₊* and
    # Ω₋ᵀ = Ω₋*, which is Combes' form
    Omega = [
        omegaminus omegaplus;
        -Ξ -transpose(omegaminus)
    ]


    A = -0.5 * J(m) * Phi' * J(n) * Phi - 1.0im * Omega
    D = [
        S zeros(size(S));
        zeros(size(S)) conj.(S)
    ]
    B = -J(m) * Phi' * J(n) * D
    C = Phi


    #Finally we permute the matrices for our convention of 'interlaced' rather than 'doubled up' creation and annihilation operators
    A = P(m) * A * P(m)'
    B = P(m) * B * P(n)'
    C = P(n) * C * P(m)'
    D = P(n) * D * P(n)'

    return QuantumStateSpace(sys.name, sys.subspaces, sys.inputs, sys.outputs, A, B, C, D, Continuous(), AnnihilationCreationBasis())

end

function ladder_only(op, where)
    if is_position(op) || is_momentum(op)
        return error("$op in $where is a position or momentum operator; the annihilation-creation basis needs a length scale to define a, so use the quadrature basis")
    end
    if !(is_create(op) || is_destroy(op))
        return error("$op (acts_on $(only(acts_on(op)))) in $where is not a bosonic operator")
    end
end

"""
    _build_quadrature_ss(sys::SLH) → QuantumStateSpace{Continuous, QuadratureBasis}

Build the quadrature-basis ABCD matrices of a linear SLH system, following Zhang et al.,
Linear quantum systems: a tutorial (arXiv:2205.04080), section 2. The state is
z = (x₁, p₁, x₂, p₂, …) with [z, zᵀ] = iΘ, and H and L are parsed into

    H = ½ zᵀGz        L = Λz

Ladder operators are read as a = (x + ip)/√2, so H and L may mix both kinds.
"""
function _build_quadrature_ss(sys::SLH)

    m = length(sys.subspaces)
    n = length(sys.L)
    Θ(k) = kron(I(k), [0 1; -1 0])

    G = Array{Any}(zeros(2m, 2m))
    # a c-number H (no operators at all, e.g. 0) is not a QAdd, and only shifts the energy
    for (term, coeff) in (sys.H isa SecondQuantizedAlgebra.QAdd ? sys.H : ())
        ops = term.ops
        #cnumbers do not contribute to Heisenberg eqs
        isempty(ops) && continue
        if length(ops) != 2
            return error("$ops term of the Hamiltonian is not quadratic")
        end
        v1 = quadrature_vector(ops[1], m, "the Hamiltonian")
        v2 = quadrature_vector(ops[2], m, "the Hamiltonian")
        M = SecondQuantizedAlgebra.to_num(coeff) * v1 * transpose(v2)
        G += M + transpose(M)
    end

    Λ = Array{Any}(zeros(n, 2m))
    for (j, channel) in enumerate(sys.L)
        # a c-number L[j] (no operators at all) is not a QAdd
        if !(channel isa SecondQuantizedAlgebra.QAdd)
            isequal(channel, 0) && continue
            return error("L[$j] = $channel is a constant (a coherent drive), which the linear ABCD model cannot represent")
        end
        for (term, coeff) in channel
            ops = term.ops
            if isempty(ops)
                return error("L[$j] contains a constant term $(SecondQuantizedAlgebra.to_num(coeff)) (a coherent drive), which the linear ABCD model cannot represent")
            end
            if length(ops) != 1
                return error("$ops term of L[$j] is not linear")
            end
            Λ[j, :] += SecondQuantizedAlgebra.to_num(coeff) * quadrature_vector(only(ops), m, "L[$j]")
        end
    end

    # conj rather than real/imag, which drop complex parameters (see conjpart)
    re(x) = (x + conj(x)) / 2
    im_(x) = (x - conj(x)) / 2im

    C = Array{Any}(zeros(2n, 2m))
    D = Array{Any}(zeros(2n, 2n))
    for k in 1:n
        C[2k - 1, :] = sqrt(Num(2)) * re.(Λ[k, :])
        C[2k, :] = sqrt(Num(2)) * im_.(Λ[k, :])
        for l in 1:n
            s = sys.S[k, l]
            D[(2k - 1):2k, (2l - 1):2l] = [re(s) -im_(s); im_(s) re(s)]
        end
    end

    # C and D are real, and so is G for a Hermitian H; simplifying them first keeps the products below real
    G, C, D = (expandsimplify.(X) for X in (G, C, D))

    # C♯ = -Θ Cᵀ Θ, with C real
    Csharp = -Θ(m) * transpose(C) * Θ(n)
    A = expandsimplify.(Θ(m) * G - Csharp * C / 2)
    B = expandsimplify.(-Csharp * D)

    return QuantumStateSpace(sys.name, sys.subspaces, sys.inputs, sys.outputs, A, B, C, D, Continuous(), QuadratureBasis())
end

# An operator as a vector over z: x and p are unit vectors, and a = (x + ip)/√2
function quadrature_vector(op, m, where)
    v = Array{Any}(zeros(2m))
    k = only(acts_on(op))
    c = 1 / sqrt(Num(2))
    if is_position(op)
        v[2k - 1] = 1
    elseif is_momentum(op)
        v[2k] = 1
    elseif is_destroy(op)
        v[(2k - 1):2k] = [c, im * c]
    elseif is_create(op)
        v[(2k - 1):2k] = [c, -im * c]
    else
        return error("$op (acts_on $k) in $where is not a bosonic operator")
    end
    return v
end

function J(n::Int)
    Id = Matrix{Int}(I, n, n)
    Ze = Matrix{Int}(zeros(n, n))
    return [
        Id Ze;
        Ze -Id
    ]
end

function P(n::Int)
    # Create permutation matrix
    PP = zeros(Int, 2n, 2n)

    # Map annihilation operators: position i → position 2i-1
    for i in 1:n
        PP[2i - 1, i] = 1
    end

    # Map creation operators: position n+i → position 2i
    for i in 1:n
        PP[2i, n + i] = 1
    end

    return PP
end

#Does not substitute operators; substitutes into the ABCD matrices
function Symbolics.substitute(sys::QuantumStateSpace, dict)
    newA = substitute_matrix(sys.A, dict)
    newB = substitute_matrix(sys.B, dict)
    newC = substitute_matrix(sys.C, dict)
    newD = substitute_matrix(sys.D, dict)
    return QuantumStateSpace(sys.name, sys.subspaces, sys.inputs, sys.outputs, newA, newB, newC, newD, sys.timeevol, sys.basis)
end

# fold = Val(true) evaluates constant subexpressions such as sqrt(0.3); Symbolics 7 leaves them
# symbolic by default. A fully numeric result is narrowed from Matrix{Number} to a concrete
# eltype, which numerics downstream (e.g. lu inside freqresp) require. Partially symbolic
# results are returned as they are.
function substitute_matrix(M, dict)
    # flattencomplex: a complex (::Number) parameter can sit inside a part of a Complex{Num}
    # entry, and substituting a complex value into a part breaks complex(re, im)
    newM = Symbolics.value.(Symbolics.substitute.(flattencomplex.(M), [dict]; fold = Val(true)))
    if all(x -> x isa Number && !(x isa Num) && !(x isa Complex{Num}), newM)
        # Array{T}(…) rather than T.(…), which keeps the Any eltype of an empty matrix (a system with no modes)
        return all(x -> iszero(imag(x)), newM) ? Array{Float64}(real.(newM)) : Array{ComplexF64}(newM)
    end
    return newM
end

"""
    QuantumStateSpace(sys::QuantumStateSpace, basis)

The same system in `basis`. A system already in `basis` is returned as it is, and one in the
annihilation-creation basis converts to the quadrature basis. The reverse is not implemented,
since a mode on a `PhaseSpace` has no annihilation operator without a length scale.
"""
QuantumStateSpace(sys::QuantumStateSpace{TE, TB}, ::TB) where {TE, TB <: BasisType} = sys
function QuantumStateSpace(sys::QuantumStateSpace, basis::BasisType)
    return error("$(sys.name): converting from $(nameof(typeof(sys.basis))) to $(nameof(typeof(basis))) is not implemented")
end

# Matrix products can nest a complex-valued Num inside a Complex{Num}, which
# Symbolics.expand chokes on (unwrap tries to build a Complex{Real}). Folding the
# imaginary unit back into a single symbolic expression makes it expandable.
flattencomplex(z::Complex) = Symbolics.wrap(Symbolics.unwrap(real(z)) + im * Symbolics.unwrap(imag(z)))
flattencomplex(z) = z

# Type piracy: Base's conj(::Complex) assumes both parts are real, and Symbolics defines
# conj(::Num) = x. SQA's to_num puts complex (::Number) parameters inside a part of a
# Complex{Num}, so conj silently dropped them. Remove once SQA/Symbolics handle this.
conjpart(x::Num) = Symbolics.symtype(Symbolics.unwrap(x)) <: Real ? x : Num(conj(Symbolics.unwrap(x)))
Base.conj(z::Complex{Num}) = Complex(conjpart(real(z)), -conjpart(imag(z)))

expandsimplify(z) = simplify(expand(flattencomplex(z)))

function QuantumStateSpace(sys::QuantumStateSpace{TE, AnnihilationCreationBasis}, ::QuadratureBasis) where {TE}

    # (x, p) = T (a, a†) for every mode and port, with x = (a + a†)/√2 and p = i(a† - a)/√2
    c = 1 / sqrt(Num(2))
    T, Tinv = c * [1 1; -im im], c * [1 im; 1 -im]
    blocks(k, M) = kron(Matrix{Int}(I, k, k), M)

    m, n = length(sys.subspaces), length(sys.inputs)
    left, right = blocks(m, T), blocks(m, Tinv)
    leftIO, rightIO = blocks(n, T), blocks(n, Tinv)

    newA = expandsimplify.(left * sys.A * right)
    newB = expandsimplify.(left * sys.B * rightIO)
    newC = expandsimplify.(leftIO * sys.C * right)
    newD = expandsimplify.(leftIO * sys.D * rightIO)

    return QuantumStateSpace(sys.name, sys.subspaces, sys.inputs, sys.outputs, newA, newB, newC, newD, sys.timeevol, QuadratureBasis())
end
