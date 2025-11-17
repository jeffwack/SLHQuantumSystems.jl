"""
    cavity(name)

Create a basic optical cavity SLH system.

Creates a single-mode cavity with detuning and decay. The cavity has one input 
and one output port with direct transmission (S=1).

# Arguments
- `name`: Symbol identifying the cavity (used for operator and parameter naming)

# Returns
- `SLH`: System with Hamiltonian H = Δ·a†a and coupling L = [κa]
# Parameters
- `κ`: Cavity decay rate
- `Δ`: Cavity detuning from driving field
"""
function cavity(name)
    #Define the filter cavity 
    hilb = FockSpace(:cavity)
    a = Destroy(hilb,:a)

    mode = OpticalMode("")

    (κ,Δ,L) = rnumbers(parameternames(mode)...)

    paramdict = Dict(zip(nameof.([κ,Δ,L]), [κ,Δ,L]))
      
    return SLH(name,
                [mode],
                paramdict,
                ["in"],
                ["out"],
                [1],
                [κ*a],
                Δ*adjoint(a)*a)
end

#=
"""
    qed_cavity(name)

Create a cavity QED system with a two-level atom.

Creates a cavity containing a two-level atom with Jaynes-Cummings coupling
and an external driving field.

# Arguments  
- `name`: Symbol identifying the system (used for operator and parameter naming)

# Returns
- `SLH`: System with cavity-atom interaction and driving

# Parameters
- `Δ`: Cavity detuning
- `g`: Atom-cavity coupling strength  
- `κ`: Cavity decay rate
- `h`: External driving amplitude
"""
function qed_cavity(name)
    hcav = FockSpace(:cav)
    hspin = NLevelSpace(:cav,2)
    hilb = SecondQuantizedAlgebra.tensor(hcav,hspin)

    # Define the operators
    a = Destroy(hilb,:a)

    σ(i,j) = Transition(hilb, :σ, i, j, 2)

    (Δ, g, κ, h) = rnumbers(:Δ,:g,:κ,:h)

    # Hamiltonian
    H = Δ*a'*a + g*(a'*σ(1,2) + a*σ(2,1)) + h*(a + a')
    L = [κ*a]

    return SLH(name,
                [:In],
                [:Out],
                [1],
                L,
                H)
end
=#

"""
    squeezing_cavity(name)

Create a squeezing cavity SLH system.

Creates a cavity that generates squeezed light through a parametric 
interaction (two-mode squeezing Hamiltonian).

# Arguments
- `name`: Symbol identifying the cavity (used for operator and parameter naming)

# Returns  
- `SLH`: System with squeezing Hamiltonian H = iϵ(a†² - a²) and coupling L = [√κ·a]

# Parameters
- `κ`: Cavity decay rate
- `ϵ`: Squeezing strength 
"""
function squeezing_cavity(name)
    hilb = FockSpace(:squeezer)
    
    mode = GenericMode("")
    
    a = Destroy(hilb, :a)

    (κ,ϵ) = rnumbers(:κ,:ϵ)

    return SLH(name,
                [mode],
                Dict(zip(nameof.([κ,ϵ]),[κ,ϵ])),
                ["in"],
                ["out"],
                [1],
                [κ*a],
                1im*ϵ*(adjoint(a)^2- a^2))
end
#=
"""
    radiation_pressure_cavity(name)

Create an optomechanical cavity with radiation pressure coupling.

Creates a two-mode system with an optical cavity mode coupled to a 
mechanical oscillator through radiation pressure.

# Arguments
- `name`: Symbol identifying the system (used for operator and parameter naming)

# Returns
- `SLH`: Optomechanical system with optical and mechanical modes

# Parameters  
- `κ`: Cavity decay rate
- `Δ`: Cavity detuning
- `Ω`: Mechanical frequency
- `g`: Optomechanical coupling strength
"""
function radiation_pressure_cavity(name)
    hilb = SecondQuantizedAlgebra.tensor(FockSpace(:optical), FockSpace(:mechanical))
    a = Destroy(hilb,:a,1) 
    b = Destroy(hilb,:b,2) 

    (κ, Δ, Ω, g) = rnumbers(:κ,:Δ,:Ω,:g)

    return SLH(name,
                [:In],
                [:Out],
                [1],
                [κ*a],
                Δ*a'*a+Ω*b'*b - g*a'*a*(b'+b))
end
=#
