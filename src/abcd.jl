#= This file defines a type for the representation of linear systems
=#

using ControlSystems: AbstractStateSpace, Continuous
import ControlSystems: system_name, input_names, output_names,
    ssdata, nstates, ninputs, noutputs, series, feedback

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
input_names(sys::QuantumStateSpace) = sys.inputs
output_names(sys::QuantumStateSpace) = sys.outputs

# ── ControlSystems.jl interface ───────────────────────────────────────────────

ssdata(sys::QuantumStateSpace) = (sys.A, sys.B, sys.C, sys.D)
nstates(sys::QuantumStateSpace) = size(sys.A, 1)
ninputs(sys::QuantumStateSpace) = size(sys.B, 2)
noutputs(sys::QuantumStateSpace) = size(sys.C, 1)

const _COMPOSE_ERROR = """not implemented"""

series(::QuantumStateSpace, ::QuantumStateSpace) = error(_COMPOSE_ERROR)
feedback(::QuantumStateSpace, ::QuantumStateSpace) = error(_COMPOSE_ERROR)
feedback(::QuantumStateSpace) = error(_COMPOSE_ERROR)

# Default: build in ladder basis then convert to quadrature
function QuantumStateSpace(sys::SLH)
    return toquadrature(_build_ladder_ss(sys))
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
            if !(is_create(op) || is_destroy(op))
                return error("$op (acts_on $(only(acts_on(op)))) in the Hamiltonian is not a creation or annihilation operator")
            end
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
            c = SecondQuantizedAlgebra.to_num(coeff)
            k = only(acts_on(op))
            if is_destroy(op)
                phiminus[j, k] += c
            elseif is_create(op)
                phiplus[j, k] += c
            else
                return error("$op in L[$j] is not a creation or annihilation operator; the ladder basis requires Fock operators")
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
