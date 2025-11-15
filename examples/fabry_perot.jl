using SecondQuantizedAlgebra
using SLHQuantumSystems
using Symbolics
using GLMakie

hilb = FockSpace(:cav)
a = Destroy(hilb, :a)

@cnumbers ω κ_L κ_R

cav = SLH("cav",[1 0; 0 1],[√κ_L*a, √κ_R*a],ω*a'*a)

cavSS = StateSpace(cav)

quadcavSS = toquadrature(cavSS)

#aasymbtfs = symbfresponse(cavSS)
#quadsymbtfs = symbfresponse(quadcavSS)

paramdict = Dict([ω => 10, κ_L => 3, κ_R => 2])

numcav = substitute(cavSS,paramdict)

numcavquad = toquadrature(numcav)

bode(numcavquad,("in1",1),("out2",1),collect(logrange(1,10000,100)))
