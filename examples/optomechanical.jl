using SLHQuantumSystems
using SecondQuantizedAlgebra
using Symbolics
using GLMakie
#=
#following 1.5.4 in Linear Dynamical Quantum Systems we will build the Hamiltonian in the quadrature basis
hilb = PhaseSpace(:light) ⊗ PhaseSpace(:mirror)

q1 = Position(hilb,:q1,1)
q2 = Momentum(hilb,:q2,1)
x = Position(hilb,:x,2)
p = Momentum(hilb,:p,2)

@rnumbers m ω g κ

H = m*ω^2*x^2/2 + p^2/(2*m) + g*x*q1
L = κ*(q1+1.0im*q2)
S = [1]

optomech = SLH(:opto, [:in], [:out], S, L, H)
#gives error: Hilbert space has non-bosinic modes.
=#

# We construct the Hamiltonian using creation and annihilation operators

hilb = FockSpace(:cavity)⊗FockSpace(:mirror)

subspaces = [OpticalMode(""),MechanicalMode("")]

a = Destroy(hilb,operatornames(subspaces[1])[1],1)
b = Destroy(hilb,operatornames(subspaces[2])[1],2)

(ω,l,κ) = rnumbers(parameternames(subspaces[1])...)
(Ω,m,Γ) = rnumbers(parameternames(subspaces[2])...)

g = rnumbers(:g)[1]#coupling parameter is defined separately 

#Hamiltonian (Chen 2013 eq 2.4)
H = Ω*b'*b - g*(b'+b)*(a' + a)
L = [κ*a,Γ*b]
S = [1 0; 0 1]

params = [ω,l,κ,Ω,m,Γ,g]
pdict = Dict(zip(nameof.(params),params))
slh = SLH("opto",subspaces,pdict,["l_in","m_in"],["l_out","m_out"], S, L, H)
aass = StateSpace(slh)
qss = toquadrature(aass)

#Now we want to substitute numerical values.

paramdict = Dict([ω => 0, l=>1, κ => sqrt(2*2*pi*100),Ω=>0.1,m=> 1000,g=>1000,Γ => 0.001])

numeric = substitute(qss,paramdict)

freq = collect(logrange(0.01,10000,1000))

N = fresponse_allIO(numeric,freq)
S = fresponse_state2output(numeric, freq, 1,2)

fig = Figure()
ax = Axis(fig[1,1],xscale=log10, yscale=log10)
scatter!(ax,freq,abs.(N[2,2]);label="2,2")
scatter!(ax,freq,abs.(N[2,1]);label="2,1")
scatter!(ax,freq,abs.(N[2,1]-N[2,2]);label="sum")
scatter!(ax,freq,abs.(S);label="signal")
axislegend(ax)
fig

