# Here we define a one-sided Fabry-Perot cavity as an SLH system, and perform some simple
# numerical experiments on it to demonstrate basic simulations with physical
# parameters

using SecondQuantizedAlgebra
using QuantumOptics
using SLHQuantumSystems

hilb = FockSpace(:cav)
a = Destroy(hilb, :a)

@cnumbers ω κ

cavity = SLH(:cav,[:in],[:out],[1],[√κ*a],ω*a'*a)

nH = to_numeric(cavity.H, FockBasis(10))
