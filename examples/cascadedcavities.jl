# # Cascaded Cavities Example
#
# This example demonstrates how to create and compose quantum systems using 
# the SLH framework by building a cascaded cavity system.
#
# ## Setup
#
# First, we import the required packages:

using SLHQuantumSystems
using SecondQuantizedAlgebra

# ## Creating the Hilbert Space and Operators
#
# We start by defining a Fock space for our cavity and the associated operators:

# Create a Hilbert space and operators
hilb = FockSpace(:cavity)

@qnumbers a::Destroy(hilb)
@cnumbers ω κ

# ## Defining System Components
#
# In the SLH framework, each system is characterized by three components:
# - **S**: Scattering matrix (direct input-output coupling)
# - **L**: Coupling vector (system-environment interaction)  
# - **H**: System Hamiltonian (internal dynamics)
#
# For a simple cavity, we define:

# Define system components 
H = ω * a' * a      # Harmonic oscillator Hamiltonian
L = [sqrt(κ) * a]   # Coupling to environment (decay)
S = [1]             # No direct scattering

# ## Building the Cascaded System
#
# Now we create two identical cavity systems and connect them in a cascade:

# Create SLH systems
cavityA = SLH(:A, [:in], [:out], S, L, H)
cavityB = SLH(:B, [:in], [:out], S, L, H)

## System Composition
#
# We concatenate the systems in a chain configuration, then apply feedback
# to connect the output of cavity A to the input of cavity B:

sys = concatenate([cavityA, cavityB], :chain)
sys = feedbackreduce(sys, :A_out, :B_in)

# ## Results
#
# Let's examine the resulting system:

println("Combined system Hamiltonian:")
println(sys.H)

println("\nSystem operators:")
println(operators(sys))

println("\nSystem parameters:")
println(parameters(sys))

