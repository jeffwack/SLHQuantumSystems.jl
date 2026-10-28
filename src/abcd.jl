#= This file defines a type for the representation of linear systems
=#

using ControlSystemsBase: AbstractStateSpace, Continuous
import ControlSystemsBase: system_name, input_names, output_names, state_names,
    ssdata, nstates, ninputs, noutputs, series, feedback, StateSpace

abstract type BasisType end

"""
Ladder (creation/annihilation) operator basis.
State vector ordered as (a₁, a₁†, a₂, a₂†, ...) after interlacing.
ABCD matrices are complex-valued.
"""
struct LadderBasis <: BasisType end

"""
Quadrature (amplitude/phase) operator basis.
State vector ordered as (X₁, P₁, X₂, P₂, ...).
All modes use the same 1/√2 unitary transform from the ladder basis.
For optical modes: X = (a+a†)/√2, P = i(a†-a)/√2.
For mechanical modes: q = (b+b†)/√2, r = i(b†-b)/√2 (dimensionless, normalized
consistently with optical quadratures; SI conversion via x_zpf, p_zpf).
A mode written with `Position` and `Momentum` has states x and p, in the units of H.
ABCD matrices are real-valued for physical systems; all entries in rad/s.
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
# Ports are Strings and states are Symbols. Each mode or port has two entries in the
# state or I/O vector, named by appending the basis suffixes.

suffixes(::LadderBasis) = ("a", "adag")
suffixes(::QuadratureBasis) = ("x", "p")

"""
    mode_names(sys) → Vector{String}

The name of each mode, in state order. Empty for an unnamed mode, which takes the
system's name on `concatenate`.
"""
mode_names(sys::Union{SLH, QuantumStateSpace}) = [mode.name for mode in sys.subspaces]

"""
    state_names(sys::QuantumStateSpace) → Vector{Symbol}

The name of each entry of the state vector: `[:arm_x, :arm_p, …]` in the quadrature basis
and `[:arm_a, :arm_adag, …]` in the ladder basis. An unnamed mode gives `[:x, :p]`.
"""
state_names(sys::QuantumStateSpace) =
    [Symbol(isempty(m) ? s : m * "_" * s) for m in mode_names(sys) for s in suffixes(sys.basis)]
state_names(sys::QuantumStateSpace, i) = state_names(sys)[i]

"""
    input_ports(sys) → Vector{String}
    output_ports(sys) → Vector{String}

The name of each input or output port, e.g. `["in"]`. Each port has two entries in the
input or output vector; see [`input_names`](@ref).
"""
input_ports(sys::Union{SLH, QuantumStateSpace}) = sys.inputs
output_ports(sys::Union{SLH, QuantumStateSpace}) = sys.outputs

"""
    input_names(sys::QuantumStateSpace) → Vector{String}
    output_names(sys::QuantumStateSpace) → Vector{String}

The name of each entry of the input or output vector, labelling the columns and rows of
D: `["in_x", "in_p", …]` in the quadrature basis and `["in_a", "in_adag", …]` in the
ladder basis. The port names alone are [`input_ports`](@ref) and [`output_ports`](@ref).
"""
input_names(sys::QuantumStateSpace) = expand_ports(input_ports(sys), sys.basis)
output_names(sys::QuantumStateSpace) = expand_ports(output_ports(sys), sys.basis)
input_names(sys::QuantumStateSpace, i) = input_names(sys)[i]
output_names(sys::QuantumStateSpace, i) = output_names(sys)[i]

expand_ports(ports, basis) = [port * "_" * s for port in ports for s in suffixes(basis)]

# A plain StateSpace, without names or basis, for the ControlSystemsBase functions that
# build a new system (series, append, step, …). ss(sys) calls this.
StateSpace(sys::QuantumStateSpace) = StateSpace(sys.A, sys.B, sys.C, sys.D, sys.timeevol)
Base.convert(::Type{StateSpace}, sys::QuantumStateSpace) = StateSpace(sys)

const _COMPOSE_ERROR = """not implemented"""

series(::QuantumStateSpace, ::QuantumStateSpace) = error(_COMPOSE_ERROR)
feedback(::QuantumStateSpace, ::QuantumStateSpace) = error(_COMPOSE_ERROR)
feedback(::QuantumStateSpace) = error(_COMPOSE_ERROR)

# Default: build in the quadrature basis, which accepts both ladder and position/momentum operators
function QuantumStateSpace(sys::SLH)
    return _build_quadrature_ss(sys)
end

function QuantumStateSpace(sys::SLH, ::LadderBasis)
    return _build_ladder_ss(sys)
end

function QuantumStateSpace(sys::SLH, ::QuadratureBasis)
    return QuantumStateSpace(sys)
end

"""
    _build_ladder_ss(sys::SLH) → QuantumStateSpace{Continuous, LadderBasis}

Build the ladder-basis ABCD matrices of a linear SLH system, following Combes eqs. 116-119.
H and L are parsed term by term into

    H = a†Ω₋a + ½(a†Ω₊a†ᵀ + h.c.)        L = Φ₋a + Φ₊a†ᵀ

"""
function _build_ladder_ss(sys::SLH)

    S = sys.S
    L = sys.L
    H = sys.H

    m = length(sys.subspaces)

    # We separate our terms into energy conserving (omegaminus) and non-conserving (omegaplus).
    omegaminus = Array{Any}(zeros(m, m))
    omegaplus = Array{Any}(zeros(m, m))

    #Each element of args is a vector of symbols/numbers corresponding to one additive term
    for (term, coeff) in H

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
            #If H is Hermitian then this is already handled
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

    Omega = [
        omegaminus omegaplus;
        -conj.(omegaplus) -conj.(omegaminus)
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

    #=#Symbolic simplification
    A = Symbolics.simplify.(A)
    B = Symbolics.simplify.(B)
    C = Symbolics.simplify.(C)
    D = Symbolics.simplify.(D)
    =#

    return QuantumStateSpace(sys.name, sys.subspaces, sys.inputs, sys.outputs, A, B, C, D, Continuous(), LadderBasis())

end

function ladder_only(op, where)
    if is_position(op) || is_momentum(op)
        return error("$op in $where is a position or momentum operator; the ladder basis needs a length scale to define a, so use the quadrature basis")
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
    for (term, coeff) in sys.H
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

    # G, C and D are real; simplifying them first keeps the products below real
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

function dampterms(L, a)
    return sum([0.5 * (Li' * commutator(a, Li) - commutator(a, Li') * Li) for Li in L])
end

function eqsofmotion(H, L, x)
    eqs = simplify.(1.0im * commutator.([H], x) + dampterms.([L], x))
    return eqs
end

#Depreciate?
# Calculates the A matrix by calculating the equations of motion and placing coefficients term by term.
#=
function makedriftA(H,L,x)
    eqs = eqsofmotion(H,L,x)
    terms = get_additive_terms.(eqs)
    args = [SymbolicUtils.arguments.(term) for term in terms]
    
    N = length(x)

    A = Array{Any}(zeros(N, N))

    #now we can build A one row at a time
    for (ii,row) in enumerate(args)
        for arglist in row
            qsym = arglist[findfirst(q->q isa SecondQuantizedAlgebra.QNumber, arglist)] #find the operator (assumes there is only one!)
            idx = findfirst(q->q==qsym, x)                                             #find the index in the state vector
            coeff = *(filter(c->!(c isa SecondQuantizedAlgebra.QNumber), arglist)...) # the coefficient is everything else
            A[ii,idx] += coeff
        end
    end

    return A
end
=#

#Does not substitute operators
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
        return all(x -> iszero(imag(x)), newM) ? Float64.(real.(newM)) : ComplexF64.(newM)
    end
    return newM
end

toquadrature(sys::QuantumStateSpace{TE, QuadratureBasis}) where {TE} = sys

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

function toquadrature(sys::QuantumStateSpace{TE, LadderBasis}) where {TE}

    blockpairs = [quadrature_transform(mode, Dict()) for mode in sys.subspaces]

    left = cat([blockpair[1] for blockpair in blockpairs]...; dims = (1, 2))
    right = cat([blockpair[2] for blockpair in blockpairs]...; dims = (1, 2))

    oldA = sys.A
    oldB = sys.B
    oldC = sys.C
    oldD = sys.D

    n_ports = length(sys.inputs)

    blockpairsIO = [quadrature_transform(GenericMode(""), Dict()) for ii in 1:n_ports]
    leftIO = cat([blockpair[1] for blockpair in blockpairsIO]...; dims = (1, 2))
    rightIO = cat([blockpair[2] for blockpair in blockpairsIO]...; dims = (1, 2))

    newA = expandsimplify.(left * oldA * right)
    newB = expandsimplify.(left * oldB * rightIO)
    newC = expandsimplify.(leftIO * oldC * right)
    newD = expandsimplify.(leftIO * oldD * rightIO)

    return QuantumStateSpace(sys.name, sys.subspaces, sys.inputs, sys.outputs, newA, newB, newC, newD, sys.timeevol, QuadratureBasis())
end
