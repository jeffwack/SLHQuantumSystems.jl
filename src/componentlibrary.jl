"""
    cavity(name)

Create a basic optical cavity SLH system.

Creates a single-mode, single-sided cavity with detuning and decay.

# Arguments
- `name`: String used for operator and parameter naming

# Returns
- `SLH`: System with Hamiltonian H = Δ·a†a and coupling L = [√κ·a]
# Parameters
- `κ`: Energy decay rate [rad/s]
- `Δ`: Cavity detuning [rad/s]

"""
function cavity(name)
    hilb = FockSpace(:cavity)
    a = Destroy(hilb, :a)
    mode = OpticalMode("")
    pnames = parameternames(mode)
    pvars = [Symbolics.variable(n) for n in pnames]
    paramdict = Dict(zip(pnames, pvars))
    κ_sym = paramdict[param_key(mode, :κ)]
    Δ_sym = paramdict[param_key(mode, :Δ)]
    opdict = Dict(:a => a)
    return SLH(
        name,
        [mode],
        ["in"],
        ["out"],
        [1],
        [sqrt(κ_sym) * a],
        Δ_sym * adjoint(a) * a
    )
end

"""
    squeezing_cavity(name)

Create a squeezing cavity SLH system.

Creates a degenerate parametric amplifier: a single-mode cavity that generates squeezed
light.

# Arguments
- `name`:  String used for operator and parameter naming

# Returns
- `SLH`: System with Hamiltonian H = (i/2)·ϵ·(a†² - a²) and coupling L = [√κ·a]

# Parameters
- `κ`: Energy decay rate [rad/s]
- `ϵ`: Pump parameter [rad/s]. 
"""
function squeezing_cavity(name)
    hilb = FockSpace(:squeezer)

    mode = GenericMode("")

    a = Destroy(hilb, :a)

    @variables κ ϵ

    opdict = Dict(zip(operator_name.([a]), [a]))

    return SLH(
        name,
        [mode],
        ["in"],
        ["out"],
        [1],
        [sqrt(κ) * a],
        (1im / 2) * ϵ * (adjoint(a)^2 - a^2)
    )
end
