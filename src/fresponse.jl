#=
This file provides methods for evaluating the frequency domain response of a linear system by numerically calculating the Laplace transform
=#

function resolvent(A,omegalist)
    return [inv(Matrix{Complex}(1.0im*omega*I - A)) for omega in omegalist]
end

function fresponse_state2output(sys::StateSpace, omegalist::Vector{Float64}, from::Int, to::Int)
    A = Matrix{Complex}(sys.A)
    C = Matrix{Complex}(sys.C)
    
    Rlist = resolvent(A,omegalist)

    return [C[to,:]'*R[:,from] for R in Rlist]
end

#below needs fixed
function fresponse_state2output(sys::StateSpace, omegalist::Vector{Float64}, from::Symbol, to::Symbol)
    j = stateidx(from)
    k = first(findall(s->s==to,sys.outputs))
    return fresponse_state2output(sys, omegalist,j,k)
end

function fresponse_allIO(sys::StateSpace, omegalist::Vector{Float64})

    A = Matrix{Complex}(sys.A)
    B = Matrix{Complex}(sys.B)
    C = Matrix{Complex}(sys.C)
    D = Matrix{Complex}(sys.D)
    
    Rlist = resolvent(A,omegalist)

    matrices = [C*R*B + D for R in Rlist]
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

