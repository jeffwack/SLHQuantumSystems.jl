struct StateSpace
    name
    inputs
    outputs
    A
    B
    C
    D
end

struct PassiveStateSpace
    name
    inputs
    outputs
    A
    B
    C
    D
end

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
    x = filter!(q->q isa Destroy,state_vector(H))
    m = length(x)

    # Now we go through every pair of annihilation operators and look for the
    # corresponding terms in the Hamiltonian
    terms = get_additive_terms(H)
    args = [SymbolicUtils.arguments(term) for term in terms]
    
    omegaminus = Array{Any}(zeros(m, m))
    #omegaminus = Array{Any}(zeros(m, m))

    for (j,aj) in enumerate(x)
        filteredterms = filter(list -> aj' in list, args)
        for (k, ak) in enumerate(x)
            term = first(filter(list -> ak in list,filteredterms))
            #types = typeof.(term)
            coef = first(filter(arg -> arg isa SymbolicUtils.BasicSymbolic, term))
            omegaminus[j,k] = coef
        end
    end

    #Now we parse L
    phiminus = Array{Any}(zeros(length(L), m))

    for (j,op) in enumerate(L)
        terms = get_additive_terms(op)
        args = [SymbolicUtils.arguments(term) for term in terms]
        for (k, ak) in enumerate(x)
            term = first(filter(list -> ak in list,args))
            coef = first(filter(arg -> arg isa SymbolicUtils.BasicSymbolic, term))
            phiminus[j,k] = coef
        end
    end
    
    newS = 

    A = -0.5*phi'*phi - 1.0im*omega
    B = -omega'*S
    C = omega
    D = S

    return PassiveStateSpace(sys.name, sys.inputs,sys.outputs, A,B,C,D)
    
end
"""
    slh2abcd(sys::SLH)
Convert a linear quantum system from SLH representation to ABCD representation.

This implements eq. 111 from Combes (2017). The name of the new system, as well as the names of its imputs and outputs will be directly inherited from the input system.
"""
function slh2abcd(sys::SLH)
    hilb = SecondQuantizedAlgebra.hilbert(sys.H)
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

    S = sys.S
    L = sys.L
    H = sys.H

    N = 2*length(L)

    if !islinear(H)
        return error("Hamiltonian contains non-quadratic terms")
    end
    
    x = state_vector(H)
    
    A = makedriftA(H,L,x)
    B = makeinputB(L,x)
    C = Matrix(B')
    D = Matrix{Int}(I, N, N)

    return StateSpace(sys.name,sys.inputs,sys.outputs,A,B,C,D)

end

## TODO: implement Combes constructions
## make separate types for passive and active linear systems?
function makephi(H,x)
    terms = get_additive_terms.(H)
    phi = Array{Any}(zeros(N, N))
    for term in terms
        args = SymbolicUtils.arguments(term)

    end
end

#Since some of the damping terms end up in A, I should create a single function
#which calculates the Heisenberg-Langevin equations of motion and then creates
#A, B, C, and D.
#
function dampterms(L,a)
    return sum([0.5*(Li'*commutator(a,Li) - commutator(a,Li')*Li) for Li in L])
end

function eqsofmotion(H,L,x)
    eqs = simplify.(1.0im*commutator.([H],x)+ dampterms.([L],x))
    return eqs
end

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

function makeinputB(L,x)
    N = 2*length(L)
    M = length(x)
    B = Array{Any}(zeros(M,N))
    
    for (ii, a) in enumerate(x)
        for (jj, c) in enumerate(L)
            B[ii,2*jj-1] = simplify(-1*commutator(a,c'))
            B[ii,2*jj] = simplify(commutator(a,c))
        end
    end

    return B
end

function stateidx(op)
    if op isa Destroy
        return 2*op.aon - 1
    elseif op isa Create
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

function Symbolics.substitute(sys::StateSpace, dict)
    newA = Symbolics.substitute.(sys.A, [dict])
    newB = Symbolics.substitute.(sys.B, [dict])
    newC = Symbolics.substitute.(sys.C, [dict])
    newD = Symbolics.substitute.(sys.D, [dict])

    return StateSpace(sys.name, sys.inputs, sys.outputs, newA, newB, newC, newD)
end

function fresponse(sys, omegalist)
    return fresponse(Matrix{Complex}(sys.A),
                     Matrix{Complex}(sys.B),
                     Matrix{Complex}(sys.C),
                     Matrix{Complex}(sys.D),
                     omegalist)
end

function fresponse(A::Matrix{Complex},B::Matrix{Complex},C::Matrix{Complex},D::Matrix{Complex}, omegalist::Vector{Float64})
    Qlist = [Matrix{Complex}(1.0im*omega*I-A) for omega in omegalist]
    matrices = [C*inv(Q)*B + D for Q in Qlist]
    P = matrices[1]
    matrixoflists = [[M[i,j] for M in matrices] for i in 1:size(P,1), j in 1:size(P,2)] 
    return matrixoflists
end

function symbfresponse(sys::StateSpace)

    @variables s
    iden = Matrix{Int}(I, size(sys.A)...)
    G =  inv(s*iden - sys.A)
    return simplify.(sys.C*G*sys.B + sys.D)
end

function toquadrature(M::Matrix)
    #the size of the matrix must be even in both dimensions
    (m,n) = size(M)
    if mod(m,2) == 0
        p = Int(m/2)
    else
        error("the size of the matrix must be even in both dimensions")
    end
    if mod(n,2) == 0
        q = Int(n/2)
    else
        error("the size of the matrix must be even in both dimensions")
    end
    
    A = 1/sqrt(2)*[1 1; -1 1]
    
    left = cat(fill(A,p)...;dims=(1,2))
    right = cat(fill(inv(A),q)...;dims=(1,2))
    
    return left*M*right
end

function toquadrature(sys::StateSpace)
    oldA = sys.A
    oldB = sys.B 
    oldC = sys.C 
    oldD = sys.D


    newA = toquadrature(oldA)
    newB = toquadrature(oldB)
    newC = toquadrature(oldC) 
    newD = toquadrature(oldD)

    return StateSpace(sys.name, sys.inputs, sys.outputs, newA, newB, newC, newD)
end
