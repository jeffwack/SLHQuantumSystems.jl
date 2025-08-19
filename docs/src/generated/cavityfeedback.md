```@meta
EditURL = "../../../examples/cavityfeedback.jl"
```

````@example cavityfeedback
using SLHQuantumSystems
using SecondQuantizedAlgebra

hilb = FockSpace(:cavity)

a = Destroy(hilb,:a)
@cnumbers ω κ

H = ω*a'*a
L = [√κ*a]
S = [1]

cavityA = SLH(:A,[:in],[:out],S,L,H)
cavityB = SLH(:B,[:in],[:out],S,L,H)

sys = concatenate([cavityA,cavityB],:loop)

sys = feedbackreduce(sys,:A_out,:B_in)
sys = feedbackreduce(sys,:B_out,:A_in)
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

