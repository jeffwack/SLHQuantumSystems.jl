struct StateSpace
    name
    inputs
    outputs
    A
    B
    C
    D
end

# This ignores S for now
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

    H = sys.H
    L = sys.L
    N = 2*length(L)

    if !islinear(H)
        return error("Hamiltonian contains non-quadratic terms")
    end
    
    x = state_vector(H)
    
    A = makedriftA(H,L,x)
    B = makeinputB(L,x)
    C = B'
    D = Matrix{Int}(I, N, N)

    return StateSpace(sys.name,sys.inputs,sys.outputs,A,B,C,D)

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
