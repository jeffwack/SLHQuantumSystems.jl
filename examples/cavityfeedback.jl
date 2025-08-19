# Here we will show the workflow for defining a system symbolicly and simulating
# it numerically

using SLHQuantumSystems
using SecondQuantizedAlgebra

# TODO: use system definitions to reduce boilerplate here

hilb = FockSpace(:cav)

a = Destroy(hilb,:a)
@cnumbers ω κ

H = ω*a'*a
L = [√κ*a]
S = [1]

cav = SLH(:C,[:in],[:out],S,L,H)

hilb2 = FockSpace(:cav) ⊗ NLevelSpace(:spin,2)

a = Destroy(hilb2,:a,1)

σ(n,k) = Transition(hilb2, :σ, n, k)

H = ω*a'*a + ω*σ(2,2) 
L = [√κ*a,√κ*σ(1,2)]
S = [1 0;
     0 1]

spincav = SLH(:SC,[:in1,:in2],[:out1,:out2],S,L,H)

sys = concatenate([cav,spincav],:loop)


