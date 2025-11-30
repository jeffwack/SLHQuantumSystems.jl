using SLHQuantumSystems
using SecondQuantizedAlgebra
using Symbolics

mode = MechanicalMode("")

@variables Ω m Γ
b = Destroy(FockSpace(:mass),operatornames(mode)[1])
paramdict = Dict(zip(nameof.([Ω,m,Γ]),[Ω,m,Γ]))

opdict = Dict(getfield(b,:name)=>b)

slh = SLH("mass",[mode],paramdict,opdict,["in"],["out"],[1],[Γ*b],Ω*b'*b)

ss = StateSpace(slh)
qss = toquadrature(ss)
