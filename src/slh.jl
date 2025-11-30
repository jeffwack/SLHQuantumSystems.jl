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
    parameters #:: #Set{SymbolicUtils.BasicSymbolic}
    operators #:: #Set{SecondQuantizedAlgebra.QNumber}
    inputs::Vector{String} #must have unique elements
    outputs::Vector{String} #must have unique elements
    S #size nxn
    L #size n
    H #has operators which act on hilbert
end

#This constructor automatically makes all the modes 'generic' and names the inputs and outputs in1, in2,... out1, out2...
function SLH(name,S,L,H)
    # We do this because SecondQuantizedAlgebra throws an error if the Hilbert space of H and L is not the same.
    hilb = SecondQuantizedAlgebra.hilbert(sum([sum(L),H]))
    
    if hilb isa SecondQuantizedAlgebra.ConcreteHilbertSpace
        n = 1
    elseif hilb isa SecondQuantizedAlgebra.ProductSpace
        n = length(hilb.spaces)
    end

    subspaces = fill(GenericMode(""),n)

    paramlist = union(get_cnumbers(H),get_cnumbers(sum(L)))
    params = Dict(zip(nameof.(paramlist),paramlist))

    oplist = union(get_qnumbers(H),get_qnumbers(sum(L)))
    ops = Dict(zip(getfield.(oplist,:name),oplist))

    m = length(L)
    if m == 1
        inputs = ["in"]
        outputs = ["out"]
    else
        inputs = ["in$j" for j in 1:m]
        outputs = ["out$j" for j in 1:m]
    end

    return SLH(name,subspaces,params,ops,inputs,outputs,S,L,H)
     
end
    
#This function is for SecondQuantizedAlgebra operators. It is intended to be general enough to work on 'any' operator
function promote_op(operator,aon_offset,new_product_space, topname)

    subspaceindex = aon_offset + operator.aon

    #these next two lines grab the necessary data to construct a new version of the operator on the product space
    #in the case of a Fock space this is just the name of the operator
    #for an NLevelSpace this is the operator name and the names of the two levels it transitions between
    middlefieldnames = fieldnames(typeof(operator))[2:end-2]
    middlefields = [getfield(operator,name) for name in middlefieldnames]

    old_op_name = popfirst!(middlefields)
    new_op_name = Symbol(topname,"_",old_op_name)
    
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
    
    #The new hilbert space is the tensor product of the old hilbert spaces
    hilb_product = tensor(old_hilberts...)

    #Now for the subspaces, parameters, operators, inputs, and outputs we will
    # prepend the name of the system to the original names before combining
    sys_names = [sys.name for sys in syslist]

    ####
    # Subspaces
    newsubspaces = vcat([[promote_name(subsys,sys.name) for subsys in sys.subspaces] for sys in syslist]...)
    ###

    ###
    # Parameters
    oldparams = [collect(values(sys.parameters)) for sys in syslist]
    newparams = [[promote_name(param,name) for param in paramlist] for (paramlist,name) in zip(oldparams,sys_names)]
    paramreplacement = [Dict([old=>new for (old,new) in zip(subsys_oldparam,subsys_newparams)]) for (subsys_oldparam,subsys_newparams) in zip(oldparams,newparams)]
    ###

    ###
    # Operators
    oldops = [collect(union(get_qnumbers(sys.H),get_qnumbers(sys.L))) for sys in syslist]
    
    #aon_offsets consists of the number of
    #'atomic' or concrete Hilbert spaces accumulated so far.
    aon_offsets = [0]
    for hilb in old_hilberts
        if hilb isa SecondQuantizedAlgebra.ProductSpace
            push!(aon_offsets, aon_offsets[end] + length(hilb.spaces))
        elseif hilb isa SecondQuantizedAlgebra.ConcreteHilbertSpace
            push!(aon_offsets, aon_offsets[end] + 1)
        else
            error("don't recognize this Hilbert space")
        end
    end
    
    newops = [[promote_op(op,offset,hilb_product,name) for op in oplist] for (oplist,name,offset) in zip(oldops,sys_names,aon_offsets)]
    opreplacement = [Dict([old => new for (old,new) in zip(subsys_oldops,subsys_newops)]) for (subsys_oldops,subsys_newops) in zip(oldops,newops)]
    ###   
  
    rulelist = [merge(paramrules,oprules) for (paramrules,oprules) in zip(paramreplacement,opreplacement)]
    
    #we concatenate all the scattering matrices block diagonally
    S = cat([sys.S for sys in syslist]...;dims=(1,2))

    L = vcat([[substitute(collapse,rules) for collapse in sys.L] for (rules,sys) in zip(rulelist,syslist)]...)

    H = sum([substitute(sys.H,rules) for (rules,sys) in zip(rulelist,syslist)])

    #We 'stack' the inputs and outputs of the systems we are combining.
    #first, we promote the names of inputs and outputs, to prevent naming collisions
    newinputs = [[sys.name*"_"*input for input in sys.inputs] for sys in syslist]
    inputs = cat(newinputs...,dims = 1)

    newoutputs = [[sys.name*"_"*output for output in sys.outputs] for sys in syslist]
    outputs = cat(newoutputs...,dims = 1)
    
    paramlist = vcat(newparams...)
    paramsdict = Dict(zip(nameof.(paramlist),paramlist))
    #println(newops)
    oplist = vcat(newops...)
    #println(oplist)
    opsdict = Dict(zip(getfield.(oplist,:name),oplist))
    
    return SLH(name,newsubspaces,paramsdict,opsdict,inputs,outputs,S,L,H)
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

    #return (A.H, Hprime)
    H = simplify(A.H + Hprime)

    return SLH(A.name,A.subspaces,A.parameters,A.operators,newinputs,newoutputs,S,L,H)
end
