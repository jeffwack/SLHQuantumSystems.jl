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

    @variables κ Δ L  

    paramdict = Dict(zip(nameof.([κ,Δ,L]), [κ,Δ,L]))
    opdict = Dict(zip(getfield.([a],:name),[a]))
      
    return SLH(name,
                [mode],
                paramdict,
                opdict,
                ["in"],
                ["out"],
                [1],
                [κ*a],
                Δ*adjoint(a)*a)
end

"""
    squeezing_cavity(name)

Create a squeezing cavity SLH system.

Creates a cavity that generates squeezed light through a parametric 
interaction (two-mode squeezing Hamiltonian).

# Arguments
- `name`: Symbol identifying the cavity (used for operator and parameter naming)

# Returns  
- `SLH`: System with squeezing Hamiltonian H = iϵ(a†² - a²) and coupling L = [κ·a]

# Parameters
- `κ`: Cavity decay rate
- `ϵ`: Squeezing strength 
"""
function squeezing_cavity(name)
    hilb = FockSpace(:squeezer)
    
    mode = GenericMode("")
    
    a = Destroy(hilb, :a)

    @variables κ ϵ
    
    opdict = Dict(zip(getfield.([a],:name),[a]))

    return SLH(name,
                [mode],
                Dict(zip(nameof.([κ,ϵ]),[κ,ϵ])),
                opdict,
                ["in"],
                ["out"],
                [1],
                [κ*a],
                1im*ϵ*(adjoint(a)^2- a^2))
end
