using SLHQuantumSystems
using SecondQuantizedAlgebra

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

hilb = FockSpace(:light)⊗FockSpace(:mirror)

al = Destroy(hilb,:l,1)
am = Destroy(hilb,:m,2)

@rnumbers ω g κ

H = ω*am'*am + g*(al' + al)*(am' + am)
L = [κ*al]
S = [1]

slh = SLH(:opto, [:in], [:out], S, L, H)
aass = StateSpace(slh)
qss = toquadrature(aass)
