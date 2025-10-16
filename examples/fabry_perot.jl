using SecondQuantizedAlgebra
using SLHQuantumSystems
using Symbolics
using GLMakie

hilb = FockSpace(:cav)
a = Destroy(hilb, :a)

@cnumbers ω κ_L κ_R

cav = SLH(:cav,[:in_L, :in_R],[:out_L, :out_R],[1 0; 0 1],[√κ_L*a, √κ_R*a],ω*a'*a)

cavSS = slh2abcd(cav)

paramdict = Dict([ω => 0, κ_L => 3, κ_R => 2])

numcav = substitute(cavSS,paramdict)

freq = collect(logrange(0.01,100,200))

tf = fresponse(numcav,freq)

