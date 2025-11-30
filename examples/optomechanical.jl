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

@variables ω l κ #consider replacing calls to parameternames
@variables Ω m Γ 

@variables g #coupling parameter is defined separately 

#Hamiltonian (Chen 2013 eq 2.4)
H = Ω*b'*b - g*(b'+b)*(a' + a)
L = [κ*a,Γ*b]
S = [1 0; 0 1]

params = [ω,l,κ,Ω,m,Γ,g]
pdict = Dict(zip(nameof.(params),params))

opdict = Dict(zip(getfield.([a,b],:name),[a,b]))

slh = SLH("opto",subspaces,pdict,opdict,["l_in","m_in"],["l_out","m_out"], S, L, H)
aass = StateSpace(slh)
qss = toquadrature(aass)

#Now we want to substitute numerical values.

paramdict = Dict([ω => 0, l=>10, κ => 20,Ω=>0.05,m=> 10,g=>100,Γ => 0.001])

numeric = substitute(qss,paramdict)

freq = collect(logrange(0.01,10000,1000))

N = fresponse_allIO(numeric,freq)
S = fresponse_state2output(numeric, freq, 2,2)


fig = Figure()
ax = Axis(fig[1,1],xscale=log10, yscale=log10)
ylims!(ax,10e-3,10e7)
scatter!(ax,freq,abs.(N[2,2]);label="shot noise")
scatter!(ax,freq,abs.(N[2,1]);label="radiation pressure noise")
scatter!(ax,freq,abs.(N[2,2]-N[2,1])./abs.(S);label="noise to signal ratio")
scatter!(ax,freq,abs.(S);label="signal response")
axislegend(ax; position = :rt)
fig
# save("optomechanical.png",fig)

