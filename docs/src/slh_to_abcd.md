# SLH to ABCD Conversion

This example demonstrates how to convert an SLH representation of a linear 
quantum system to ABCD state-space form. This transformation is useful for
calculating the frequency response of the system.

## Linear Optomechanical System

We'll use a linear optomechanical system as our example, a model of a laser
interferometer detecting gravitational waves. This model is discussed in [a
review article written by Yanbei Chen](https://arxiv.org/abs/1302.1924)

```@example slh_to_abcd
using SLHQuantumSystems
using QuantumCumulants

# Create Hilbert spaces
hilb = FockSpace(:optical) ⊗ FockSpace(:mechanical)

# Create operators
a = Destroy(hilb,:a,1)
b = Destroy(hilb,:b,2)

# Define parameters
@cnumbers ω Δ g κ 
```

The system has:
- Hamiltonian: `H = ω*b†*b + Δ*a†*a - g*(a† + a)*(b† + b)` (linear optomechanical coupling)
- Dissipation: `L = [√κ * a]` (optical cavity loss)
- Scattering: `S = [1]` (trivial scattering)

```@example slh_to_abcd
# Build the system
H = ω*b'*b + Δ*a'*a + g*(a' + a)*(b' + b)
L = [√κ * a]
S = [1]

# Create SLH system
sys = SLH(:optomech, [:in], [:out], S, L, H)

println("System created:")
println("Name: ", sys.name)
println("Inputs: ", sys.inputs)
println("Outputs: ", sys.outputs)
println("Hamiltonian: ", sys.H)
println("Dissipation: ", sys.L)
println("Scattering: ", sys.S)
```

## Converting to ABCD Form
First, we will create an ordering for the systems operators
```@example slh_to_abcd
ops = ordered_ops(sys)
```

Now, we want to calculate the equations of motion for this system. For now we
will ignore the scattering matrix, and use the Heisenberg-Lengevin equations of
motion,

$\dot{a} = i[H,a] + \sum_i (b_i^\dagger + \frac{1}{2}L_i^\dagger) [a,L_i] - [a,L_i^\dagger](b_i + \frac{1}{2}L_i)$

Which is a multi-input generalization of equation 2.12 in [Gardiner and
Collett](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.31.3761). TODO:
adress numerical factors.

We can see that the two terms in the parenthesis correspond to input and
damping. The damping terms should end up in A and the input terms in B. Let's
generate A one row at a time.

```@example slh_to_abcd
op = ops[1]

using Symbolics
opdot = 1.0im*simplify(commutator(sys.H,a))
```

```@example slh_to_abcd
# Convert to ABCD representation
try
    (A, B, C, D) = SLH2ABCD(sys)
    
    println("\nABCD matrices:")
    println("A (drift): ", A)
    println("B (input coupling): ", B)
    println("C (output coupling): ", C)
    println("D (feedthrough): ", D)
catch e
    println("Error during conversion: ", e)
end
```

The ABCD matrices represent the system in state-space form:
- **A**: Drift matrix (describes internal system dynamics)
- **B**: Input coupling matrix (how inputs affect the system)
- **C**: Output coupling matrix (how system state affects outputs)
- **D**: Direct feedthrough matrix (direct input-to-output coupling)
