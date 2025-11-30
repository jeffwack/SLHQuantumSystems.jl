#= This file defines a type for the representation of linear systems
=#

struct StateSpace
    name
    subspaces
    parameters
    inputs
    outputs
    A
    B
    C
    D
end

#This uses the Combes method of calculating Phi and Omega (rather than directly calculating the equations of motion)
function StateSpace(sys::SLH)

    S = sys.S
    L = sys.L
    H = sys.H
    
    # First we need to ensure that this is a linear quantum system.
    # This means is consists of bosonic modes with quadratic couplings
    hilb = SecondQuantizedAlgebra.hilbert(H)
    if hilb isa SecondQuantizedAlgebra.ProductSpace    
        for subspace in hilb.spaces
            if !(subspace isa FockSpace)
                return error("Hilbert space contains non-bosonic modes")
            end
        end
    else
        if !(hilb isa FockSpace)
            return error("Hilbert space contains non-bosonic modes")
        end
    end
    
    # We create a vector which contains the creation and annihilation operators
    # of the system. For now we assume we are in that basis rather than
    # quadratures.
    x = filter(q->q isa Destroy,state_vector(H))
    m = length(x)
    
    # Now we go through every pair of annihilation operators and look for the
    # corresponding terms in the Hamiltonian
    terms = get_additive_terms(H)
    
    args = [SymbolicUtils.arguments(term) for term in terms]
    #println([typeof.(term) for term in args])

    # We want to separate our terms into energy conserving (omegaminus) and non-conserving (omegaplus).
    # The operator ordering of second quantized algebra is only applied to non-commuting operators.
    
    omegaminus = Array{Any}(zeros(m, m))
    omegaplus = Array{Any}(zeros(m, m))

    #organized by term
    for term in args
        #println(term)
        # first we will check that this term is quadratic in the operators
        if length(filter(q->q isa SecondQuantizedAlgebra.QNumber,term)) != 2
            return error("$term of the Hamiltonian is not quadratic")
        end
        #now we can determine what 'kind' of term this is by counting the number of creation operators
        creators = filter(q->q isa SecondQuantizedAlgebra.Create,term)
        destroyers = filter(q->q isa SecondQuantizedAlgebra.Destroy,term)
        if length(creators) == 1
            cr = first(creators)
            j = first(findall(q->q==cr',x))

            de = first(destroyers)
            k = first(findall(q->q==de,x))
            omegaminus[j,k] = first(filter(arg -> arg isa Number || arg isa SymbolicUtils.BasicSymbolic, term)) ##TODO: remove 'or' here and in analogous statements. Where is BasicSymbolic being introduced?
        elseif length(creators) == 2 #&& creators[1] == creators[2]
            cr1 = creators[1]
            j = first(findall(q->q==cr1',x))

            cr2 = creators[2]
            k = first(findall(q->q==cr2',x))
            
            coef = first(filter(arg -> arg isa Number || arg isa SymbolicUtils.BasicSymbolic, term))
            omegaplus[j,k] = coef
            omegaplus[k,j] = coef
        elseif length(creators) == 0 #&& destroyers[1] == destroyers[2]
            de = first(destroyers)
            j = first(findall(q->q==de,x))
            #If H is Hermitian then this is already handled
        else
            return error("don't know how to handle $term")
        end
    end

    #organized by operator
    #=
    for (j,aj) in enumerate(x)
        filteredterms = filter(list -> aj' in list, args)
        for (k, ak) in enumerate(x)
            term = first(filter(list -> ak in list,filteredterms))
            #types = typeof.(term)
            coef = first(filter(arg -> arg isa SymbolicUtils.BasicSymbolic, term))
            omegaminus[j,k] = coef
        end
    end
    =#
    #Now we parse L
    n = length(L)
    phiminus = Array{Any}(zeros(n, m))
    phiplus = Array{Any}(zeros(n,m))

    for (j,op) in enumerate(L)
        terms = get_additive_terms(op)
        args = [SymbolicUtils.arguments(term) for term in terms]
        for term in args
            #first make sure each term is linear in operators
            ops = filter(q->q isa SecondQuantizedAlgebra.QNumber,term)
            if length(ops) != 1
                return error("$term of L is not linear")
            end
            op = first(ops)
            #=
            println(SecondQuantizedAlgebra.hilbert(op))
            println(getfield(op,:name))
            println(getfield(op,:aon))
            println(SecondQuantizedAlgebra.hilbert.(x))
            println(getfield.(x,:name))
            println(getfield.(x,:aon))
            println(findall(q->q==op,x))
            =#
            if op isa SecondQuantizedAlgebra.Destroy
                k = first(findall(q->q==op,x))
                phiminus[j,k] = first(filter(arg -> arg isa Number || arg isa SymbolicUtils.BasicSymbolic, term))
            else
                k = first(findall(q->q==op',x))
                phiplus[j,k] = first(filter(arg -> arg isa Number || arg isa SymbolicUtils.BasicSymbolic,term))
            end
        end
    end

    #Now we need to construct the matrices according to Combes eqs 116-119
    
    Phi = [phiminus phiplus;
            conj.(phiplus) conj.(phiminus)]

    Omega = [omegaminus omegaplus;
            -conj.(omegaplus) -conj.(omegaminus)]

    
    A = -0.5*J(m)*Phi'*J(n)*Phi - 1.0im*Omega
    D = [S zeros(size(S));
        zeros(size(S)) conj.(S)]
    B = -J(m)*Phi'*J(n)*D
    C = Phi
    

    #Finally we permute the matrices for our convention of 'interlaced' rather than 'doubled up' creation and annihilation operators
    A = P(m)*A*P(m)'
    B = P(m)*B*P(n)'
    C = P(n)*C*P(m)'
    D = P(n)*D*P(n)'

    #=#Symbolic simplification
    A = Symbolics.simplify.(A)
    B = Symbolics.simplify.(B)
    C = Symbolics.simplify.(C)
    D = Symbolics.simplify.(D)
    =#

    return StateSpace(sys.name, sys.subspaces,sys.parameters, sys.inputs,sys.outputs, A,B,C,D)
    
end

function J(n::Int)
    Id = Matrix{Int}(I, n, n)
    Ze = Matrix{Int}(zeros(n, n))
    return [Id Ze;
            Ze -Id]
end

function P(n::Int)
    # Create permutation matrix
    PP = zeros(Int, 2n, 2n)
    
    # Map annihilation operators: position i → position 2i-1
    for i in 1:n
        PP[2i-1, i] = 1
    end
    
    # Map creation operators: position n+i → position 2i
    for i in 1:n
        PP[2i, n+i] = 1
    end
    
    return PP
end

function dampterms(L,a)
    return sum([0.5*(Li'*commutator(a,Li) - commutator(a,Li')*Li) for Li in L])
end

function eqsofmotion(H,L,x)
    eqs = simplify.(1.0im*commutator.([H],x)+ dampterms.([L],x))
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

function stateidx(op)
    if op isa Destroy || op isa Position
        return 2*op.aon - 1
    elseif op isa Create || op isa Momentum
        return 2*op.aon
    else
        return error("Can't assign index to non-bosonic operator")
    end
end

function state_vector(H)
    hilb = SecondQuantizedAlgebra.hilbert(H)

    ops = collect(get_qnumbers(H))

    # we want to order the operators based on the order of their hilbert spaces
    # in the product space, that is by aon
    x = sort(ops, by = stateidx)
    return x
end

#Does not substitute operators
function Symbolics.substitute(sys::StateSpace, dict)
    newA = Symbolics.value.(Symbolics.substitute.(sys.A, [dict]))
    newB = Symbolics.value.(Symbolics.substitute.(sys.B, [dict]))
    newC = Symbolics.value.(Symbolics.substitute.(sys.C, [dict]))
    newD = Symbolics.value.(Symbolics.substitute.(sys.D, [dict]))
    params = sys.parameters
    newparams = Dict([(key,dict[params[key]]) for key in keys(params)])
    return StateSpace(sys.name, sys.subspaces, newparams,sys.inputs, sys.outputs, newA, newB, newC, newD)
end

function toquadrature(sys::StateSpace)

    blockpairs = [quadratureblocks(sys,mode) for mode in sys.subspaces]

    left = cat([blockpair[1] for blockpair in blockpairs]...;dims=(1,2))
    right = cat([blockpair[2] for blockpair in blockpairs]...;dims=(1,2))
    
    oldA = sys.A
    oldB = sys.B 
    oldC = sys.C 
    oldD = sys.D

    n_ports = length(sys.inputs)

    blockpairsIO = [quadratureblocks(sys,GenericMode("")) for ii in 1:n_ports]
    leftIO = cat([blockpair[1] for blockpair in blockpairsIO]...;dims=(1,2))
    rightIO = cat([blockpair[2] for blockpair in blockpairsIO]...;dims=(1,2))

    newA = simplify.(expand.(left*oldA*right))
    newB = simplify.(expand.(left*oldB*rightIO))
    newC = simplify.(expand.(leftIO*oldC*right))
    newD = simplify.(expand.(leftIO*oldD*rightIO))

    return StateSpace(sys.name, sys.subspaces,sys.parameters, sys.inputs, sys.outputs, newA, newB, newC, newD)
end

