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
struct SLH{S,L,H}
    name::Symbol #should be a symbol
    inputs::Vector{Symbol} #must have unique elements
    outputs::Vector{Symbol} #must have unique elements
    S::S #size nxn
    L::L #size n
    H::H #has operators which act on hilbert
end


"""
operators(sys)

returns all the quantum operators contained in the system's Hamiltonian.
"""
function operators(sys)
    return get_qsymbols(sys.H)
end

"""
parameters(sys)

returns all the symbolic numbers contained in the system's Hamiltonian and coupling vector L.
"""
function parameters(sys::SLH)
    return union(get_numsymbols(sys.H),get_numsymbols(sum(sys.L)))
end


function promote(parameter::SymbolicUtils.BasicSymbolic, topname)
    old_sym = parameter.metadata[Symbolics.VariableSource][2]
    new_sym = Symbol(topname,:_,old_sym)
    return cnumber(new_sym)
end

function promote_name(hilb::SecondQuantizedAlgebra.ConcreteHilbertSpace, name)
    names = fieldnames(typeof(hilb))
    fields = [getfield(hilb,name) for name in names] 
    
    oldname = pop!(fields)

    newname = Symbol(name,:_,oldname)
    pushfirst!(fields,newname)

    return typeof(hilb).name.wrapper(fields...)
end


function promote_name(hilb::SecondQuantizedAlgebra.ProductSpace, name)
    #We want to create a new hilbert space where the names of all the subspaces 
    #has name prepended to it. 
    
    new_spaces = promote.(hilb.spaces, [name])

    return tensor(new_spaces...)
end

#This function is for SecondQuantizedAlgebra operators
function promote_op(operator,old_product_space,new_product_space, topname)
    #we identify the hilbert space which our operator acts on, which should be a subspace of the product space
    
    if hasproperty(operator.hilbert,:spaces) #then it is a product space
        subspace = operator.hilbert.spaces[operator.aon]
    else #it is not a product space
        subspace = operator.hilbert
    end

    #this identifies the operator's hilbert space as a subspace of the product space
    subspaceindex = findfirst(isequal(subspace),old_product_space.spaces)

    if isnothing(subspaceindex)
        error("$operator does not act on a subspace of $product_space")
    end

    #these next two lines grab the necessary data to construct a new version of the operator on the product space
    #in the case of a Fock space this is just the name of the operator
    #for an NLevelSpace this is the operator name and the names of the two levels it transitions between
    middlefieldnames = fieldnames(typeof(operator))[2:end-2]
    middlefields = [getfield(operator,name) for name in middlefieldnames]

    old_op_name = pop!(middlefields)
    new_op_name = Symbol(topname,:_,old_op_name)

    #this calls the operator constructor with the old 'middle data' but on the larger hilbert space
    return typeof(operator).name.wrapper(new_product_space,new_op_name, middlefields...,subspaceindex)
end

"""
concatenate(name, syslist::Vector{SLH})

creates a composite system with no interconnections. Combes eq. 59

When systems are concatenated, the names of their inputs, outputs, operators,
parameters, and Hilbert spaces are 'promoted' by prepending the name of the 
system to the existing name. This prevents name collisions as long as all 
SLHSystems are created with a unique name.
"""
function concatenate(syslist,name)
    old_hilberts = [SecondQuantizedAlgebra.hilbert(sys.H) for sys in syslist]
    sys_names = [sys.name for sys in syslist]
    hilb_product = tensor(promote_name.(old_hilberts,sys_names)...)
    old_hilb_product = tensor(old_hilberts...) # we will pass this to promote to aid in

    #hilb_product = SecondQuantizedAlgebra.tensor([SecondQuantizedAlgebra.hilbert(sys.H) for sys in syslist]...)
    
    #We 'stack' the inputs and outputs of the systems we are combining.
    #first, we promote the names of inputs and outputs, to prevent naming collisions
    newinputs = [[Symbol(sys.name,:_,input) for input in sys.inputs] for sys in syslist]
    inputs = cat(newinputs...,dims = 1)

    newoutputs = [[Symbol(sys.name,:_,output) for output in sys.outputs] for sys in syslist]
    outputs = cat(newoutputs...,dims = 1)

    #next, we concate all the scattering matrices block diagonally
    Slist = [sys.S for sys in syslist]
    S = cat(Slist...;dims=(1,2))

    
    #find and promote old operators to new larger Hilbert space
    oldopsplusname = [(collect(operators(sys)),sys.name) for sys in syslist]
    #println(oldopsplusname)
    newops = [[promote_op(op,old_hilb_product,hilb_product,name) for op in oplist] for (oplist,name) in oldopsplusname]
    
    oldops = [tup[1] for tup in oldopsplusname]

    #promote names of parameters
    oldparamsplusname = [(collect(parameters(sys)),sys.name) for sys in syslist]
    
    newparams = [[promote(param,name) for param in paramlist] for (paramlist,name) in oldparamsplusname]
    
    oldparams = [tup[1] for tup in oldparamsplusname]

    oldsyms = [cat(ops,params,dims=1) for (ops,params) in zip(oldops,oldparams)]
    newsyms = [cat(ops,params,dims=1) for (ops,params) in zip(newops,newparams)] 
    
    rulelist = [Dict([old => new for (old,new) in zip(oldoplist,newoplist)]) for (oldoplist,newoplist) in zip(oldsyms,newsyms)]
    newHs = [substitute(sys.H,rules) for (rules,sys) in zip(rulelist,syslist)]

    H = sum(newHs)

    newLs = [[substitute(collapse,rules) for collapse in sys.L] for (rules,sys) in zip(rulelist,syslist)]

    L = vcat(newLs...)

    return SLH(name,inputs, outputs,S,L,H)
end

"""
feedbackreduce(A::SLH,output,input)

Connects the output port to the input port, reducing the number of outputs and inputs by one each. Combes eq 61.
"""
function feedbackreduce(A,output, input)

    x = findfirst(isequal(output),A.outputs)
    y = findfirst(isequal(input), A.inputs)

    newoutputs = A.outputs[eachindex(A.outputs) .!= x]
    newinputs = A.inputs[eachindex(A.inputs) .!= y]

    Sxbarybar = A.S[1:end .!= x, 1:end .!= y]
    Sxbary = A.S[1:end .!= x, y]
    Sxybar = permutedims(A.S[x, 1:end .!= y])
    Sxy = A.S[x,y]

    Sy = A.S[:,y]

    Lxbar = A.L[1:end .!= x]
    Lx = A.L[x]

    S = Sxbarybar + Sxbary*(1-Sxy)^(-1)*Sxybar
    #Have to use fill here because the usual broadcast syntax on operator does not work
    L = Lxbar +  Sxbary .* fill(((1-Sxy)^(-1)*Lx),size(Sxbary))

    term1 = adjoint(A.L)*Sy
    term2 = (1-Sxy)^(-1)
    term3 = Lx
    termA = adjoint(Lx)
    termB = ((1-adjoint(Sxy))^(-1))
    termC = (adjoint(Sy)*A.L)

    Hprime = 1/(2im)*(term1*term2*term3-termA*termB*termC)

    H = A.H + Hprime

    return SLH(A.name,newinputs,newoutputs,S,L,H)
end
