function slh2abcd(sys::SLH)
    hilb = SecondQuantizedAlgebra.hilbert(sys.H)
    for subspace in hilb.spaces
        if !(subspace isa FockSpace)
            return error("Hilbert space contains non-bosonic modes")
        end
    end

    H = sys.H

    if !islinear(H)
        return error("Hamiltonian contains non-quadratic terms")
    end

    terms = get_additive_terms(H)
    
    x = state_vector(H)

    for term in terms
        #we want to extract the numerical coefficient and populate the relevant
        #term of the matrix.
    end
end

function makeA(H,x)
    eqs = simplify.(1.0im*commutator.([H],x))
    terms = get_additive_terms.(eqs)
    args = [SymbolicUtils.arguments.(term) for term in terms]
    
    N = length(x)

    A = Array{Any}(zeros(4, 4))

    #now we can build A one row at a time
    for (ii,row) in enumerate(args)
        for arg in row
            qsym = first(filter(q->q isa SecondQuantizedAlgebra.QNumber,arg))
            idx = findfirst(q->q==qsym, x) #find the index in the state vector
            coeff = *(filter(c->c isa SymbolicUtils.Symbolic,arg)...)
            A[ii,idx] = coeff
        end
    end

    return A
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


