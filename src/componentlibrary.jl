"""
    cavity(name)

Create a basic optical cavity SLH system.

Creates a single-mode, single-sided cavity with detuning and decay.

# Arguments
- `name`: The system name. Its parameters are named `κ` and `Δ` until `concatenate`
  promotes them to `name_κ` and `name_Δ`.

# Returns
- `SLH`: System with Hamiltonian H = Δ·a†a and coupling L = [√κ·a]
# Parameters
- `κ`: Energy decay rate [rad/s]
- `Δ`: Cavity detuning [rad/s]

"""
function cavity(name)
    subspaces = [OpticalMode("")]
    a, = operators(subspaces, "")
    @variables κ Δ
    return SLH(
        name,
        subspaces,
        ["in"],
        ["out"],
        [1],
        [sqrt(κ) * a],
        Δ * adjoint(a) * a
    )
end

"""
    squeezing_cavity(name)

Create a squeezing cavity SLH system.

Creates a degenerate parametric amplifier: a single-mode cavity that generates squeezed
light.

# Arguments
- `name`: The system name. Its parameters are named `κ` and `ϵ` until `concatenate`
  promotes them to `name_κ` and `name_ϵ`.

# Returns
- `SLH`: System with Hamiltonian H = (i/2)·ϵ·(a†² - a²) and coupling L = [√κ·a]

# Parameters
- `κ`: Energy decay rate [rad/s]
- `ϵ`: Pump parameter [rad/s].
"""
function squeezing_cavity(name)
    subspaces = [OpticalMode("")]
    a, = operators(subspaces, "")
    @variables κ ϵ
    return SLH(
        name,
        subspaces,
        ["in"],
        ["out"],
        [1],
        [sqrt(κ) * a],
        (1im / 2) * ϵ * (adjoint(a)^2 - a^2)
    )
end
