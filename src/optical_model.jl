"""
    Optical model parser and SLH compiler.

Loads backend-agnostic `optical-model` YAML files (the same format used by
the Python SFLU compiler) and compiles them into SLH systems.

Main entry points:
- `load_optical_model(path)` — parse YAML → `OpticalModel` struct
- `compile_to_slh(model; P_arm)` — compile `OpticalModel` → `SLH`
"""

using YAML
using PhysicalConstants.CODATA2018: ReducedPlanckConstant as ℏ_SI, SpeedOfLightInVacuum as c_SI

# ============================================================================
# Data types for the parsed optical model
# ============================================================================

"""
Mirror component from the optical-model YAML.
"""
struct OpticalMirror
    name::String
    T::Float64              # power transmissivity
    L::Float64              # power loss
    role::String            # "end", "input", "signal_extraction", "filter_cavity", "generic"
    loss_ports::Bool
    lambda_m::Union{Float64,Nothing}
    suscept::Union{Dict,Nothing}   # mechanical susceptibility (e.g. Dict("type"=>"FreeMass","mass"=>40))
    loss_in_transmission::Bool
end

"""
Squeezer component from the optical-model YAML.
"""
struct OpticalSqueezer
    name::String
    sqz_dB::Float64
    sqz_angle_deg::Float64
end

"""
A beam path (space) connecting two component ports.
"""
struct OpticalSpace
    name::String
    from_comp::String       # component name
    from_port::String       # port name (e.g. "front")
    to_comp::String
    to_port::String
    L_m::Float64            # length in meters
    detune_rad::Float64     # one-way detuning in radians
    gouy_rad::Float64       # one-way Gouy phase in radians
end

"""
A point-to-point connection (e.g. squeezer injection).
"""
struct OpticalConnection
    from_comp::String
    from_port::String
    to_comp::String
    to_port::String
    label::String
end

"""
Strain excitation entry.
"""
struct StrainExcitation
    component::String
    scale::Float64
end

"""
Complete parsed optical model.
"""
struct OpticalModel
    description::String
    mirrors::Dict{String,OpticalMirror}
    squeezers::Dict{String,OpticalSqueezer}
    spaces::Dict{String,OpticalSpace}
    connections::Vector{OpticalConnection}
    strain_excitation::Vector{StrainExcitation}
    loss_groups::Dict{String,Vector{String}}
    lambda_m::Float64       # laser wavelength (from first mirror with lambda_m)
end

# ============================================================================
# YAML parser
# ============================================================================

"""
    load_optical_model(path::AbstractString) → OpticalModel

Parse an `optical-model` YAML file into a structured `OpticalModel`.
"""
function load_optical_model(path::AbstractString)
    data = YAML.load_file(path)
    return load_optical_model(data)
end

"""
    load_optical_model(data::Dict) → OpticalModel

Parse a pre-loaded optical-model dict into a structured `OpticalModel`.
"""
function load_optical_model(data::Dict)
    kind = get(data, "kind", "")
    kind == "optical-model" || error("Expected kind 'optical-model', got '$kind'")

    description = get(data, "description", "")
    components = get(data, "components", Dict())

    # --- Parse components ---
    mirrors = Dict{String,OpticalMirror}()
    squeezers = Dict{String,OpticalSqueezer}()
    lambda_m_global = nothing

    for (name, comp) in components
        ctype = comp["type"]
        if ctype == "Mirror"
            params = get(comp, "params", Dict())
            T = get(params, "T", 0.0)
            L = get(params, "L", 0.0)
            role = get(comp, "role", "generic")
            loss_ports = get(comp, "loss_ports", false)
            lam = get(params, "lambda_m", nothing)
            if lam !== nothing && lambda_m_global === nothing
                lambda_m_global = lam
            end

            suscept = get(params, "suscept", nothing)
            if suscept === nothing && haskey(params, "mass")
                suscept = Dict("type" => "FreeMass", "mass" => params["mass"])
            end

            lit = get(params, "loss_in_transmission", false)

            mirrors[name] = OpticalMirror(name, T, L, role, loss_ports, lam, suscept, lit)

        elseif ctype == "Squeezer"
            params = get(comp, "params", Dict())
            sqz_dB = get(params, "sqz_dB", 0.0)
            sqz_angle_deg = get(params, "sqz_angle_deg", 0.0)
            squeezers[name] = OpticalSqueezer(name, sqz_dB, sqz_angle_deg)
        end
        # Other component types (Detector, Laser, etc.) are not yet needed for SLH
    end

    if lambda_m_global === nothing
        lambda_m_global = 1.064e-6  # default Nd:YAG
    end

    # --- Parse spaces ---
    spaces_raw = get(data, "spaces", Dict())
    spaces = Dict{String,OpticalSpace}()
    for (name, sdef) in spaces_raw
        from_str = sdef["from"]   # e.g. "ITM.front"
        to_str = sdef["to"]
        from_comp, from_port = _split_port(from_str)
        to_comp, to_port = _split_port(to_str)
        params = get(sdef, "params", Dict())
        L_m = get(params, "L_m", 0.0)
        detune_rad = get(params, "detune_rad", 0.0)
        gouy_rad = get(params, "gouy_rad", 0.0)
        spaces[name] = OpticalSpace(name, from_comp, from_port, to_comp, to_port,
                                     L_m, detune_rad, gouy_rad)
    end

    # --- Parse connections ---
    conns_raw = get(data, "connections", [])
    connections = OpticalConnection[]
    for conn in conns_raw
        from_comp, from_port = _split_port(conn["from"])
        to_comp, to_port = _split_port(conn["to"])
        label = get(conn, "label", "")
        push!(connections, OpticalConnection(from_comp, from_port, to_comp, to_port, label))
    end

    # --- Parse strain excitation ---
    strain_raw = get(data, "strain_excitation", [])
    strain_exc = StrainExcitation[]
    for exc in strain_raw
        comp = exc["component"]
        scale = get(exc, "scale", 1/sqrt(2))
        push!(strain_exc, StrainExcitation(comp, scale))
    end

    # --- Parse loss groups ---
    loss_groups_raw = get(data, "loss_groups", Dict())
    loss_groups = Dict{String,Vector{String}}()
    for (gname, members) in loss_groups_raw
        loss_groups[gname] = Vector{String}(members)
    end

    return OpticalModel(description, mirrors, squeezers, spaces, connections,
                         strain_exc, loss_groups, lambda_m_global)
end

"""Split "Component.port" into ("Component", "port")."""
function _split_port(s::AbstractString)
    idx = findlast('.', s)
    idx === nothing && error("Port string '$s' must be in 'Component.port' format")
    return (s[1:idx-1], s[idx+1:end])
end

# ============================================================================
# Physical parameter derivation
# ============================================================================

"""
    derive_cavity_params(model::OpticalModel, space_name::String)

Compute the cavity decay rate κ² [rad/s] for a cavity defined by a space.
The ITM (input coupler) is inferred from the space endpoints based on mirror roles.

Returns `(κ_itm, κ_etm_loss, L_m, detune_rad)` where:
- `κ_itm = √(T_itm · c / (4·L))` — input coupler amplitude decay rate
- `κ_etm_loss = √(L_etm · c / (4·L))` — ETM loss amplitude decay rate
- `L_m` — cavity length
- `detune_rad` — cavity detuning
"""
function derive_cavity_params(model::OpticalModel, space_name::String)
    space = model.spaces[space_name]
    c = c_SI.val

    from_mirror = get(model.mirrors, space.from_comp, nothing)
    to_mirror = get(model.mirrors, space.to_comp, nothing)

    from_mirror !== nothing && to_mirror !== nothing ||
        error("Space '$space_name' endpoints must both be mirrors")

    # The input coupler is the mirror with role "input"; the end mirror has role "end"
    if from_mirror.role == "input"
        itm, etm = from_mirror, to_mirror
    elseif to_mirror.role == "input"
        itm, etm = to_mirror, from_mirror
    else
        # Fallback: higher T is the coupler
        if from_mirror.T >= to_mirror.T
            itm, etm = from_mirror, to_mirror
        else
            itm, etm = to_mirror, from_mirror
        end
    end

    L = space.L_m
    κ_itm = sqrt(itm.T * c / (4 * L))
    κ_etm_loss = sqrt(etm.L * c / (4 * L))

    return (κ_itm=κ_itm, κ_etm_loss=κ_etm_loss, L_m=L, detune_rad=space.detune_rad,
            itm=itm, etm=etm)
end

"""
    derive_optomechanical_params(model::OpticalModel, mirror::OpticalMirror, P_arm::Real)

Compute optomechanical coupling parameters for a mirror with a mechanical susceptibility.

Returns `(g_om, mass, Ω_mech)` where:
- `g_om` — linearized optomechanical coupling [rad/s]
- `mass` — mirror mass [kg]
- `Ω_mech` — mechanical frequency [rad/s]
"""
function derive_optomechanical_params(model::OpticalModel, mirror::OpticalMirror, P_arm::Real)
    mirror.suscept !== nothing || error("Mirror '$(mirror.name)' has no mechanical susceptibility")

    mass = mirror.suscept["mass"]
    c = c_SI.val
    ℏ = ℏ_SI.val
    λ = model.lambda_m
    ω_l = 2π * c / λ

    # Find the arm length from the space containing this mirror
    L_arm = nothing
    for (_, space) in model.spaces
        if space.from_comp == mirror.name || space.to_comp == mirror.name
            L_arm = space.L_m
            break
        end
    end
    L_arm !== nothing || error("Could not find arm length for mirror '$(mirror.name)'")

    # Mechanical frequency — for FreeMass use a small pendulum frequency
    stype = get(mirror.suscept, "type", "FreeMass")
    if stype == "FreeMass"
        Ω_mech = 2π * 1.0   # 1 Hz pendulum (standard convention)
    else
        Ω_mech = get(mirror.suscept, "frequency", 2π * 1.0)
    end

    # Zero-point fluctuation
    x_zpf = sqrt(ℏ / (2 * mass * Ω_mech))

    # Intracavity photon number
    N_bar = P_arm / (ℏ * ω_l)

    # Single-photon coupling
    g0 = (ω_l / L_arm) * x_zpf

    # Linearized coupling
    g_om = g0 * sqrt(N_bar)

    return (g_om=g_om, mass=mass, Ω_mech=Ω_mech, x_zpf=x_zpf, N_bar=N_bar,
            g0=g0, ω_l=ω_l, L_arm=L_arm)
end

# ============================================================================
# SLH compiler
# ============================================================================

"""
    compile_to_slh(model::OpticalModel; P_arm::Real)

Compile an `OpticalModel` into an SLH system.

# Arguments
- `model`: Parsed optical model from `load_optical_model`
- `P_arm`: Circulating arm power in watts (must be provided externally since
  it depends on laser power and recycling gains not in the optical model)

# Returns
- `(slh, params)` where `slh` is the symbolic SLH system and `params` is a
  Dict mapping symbolic variables to their numerical values.
"""
function compile_to_slh(model::OpticalModel; P_arm::Real)
    c = c_SI.val
    ℏ = ℏ_SI.val
    λ = model.lambda_m
    ω_l = 2π * c / λ

    # Identify the arm cavity space (contains the end mirror with suscept)
    arm_space_name = nothing
    arm_etm = nothing
    for (sname, space) in model.spaces
        from_m = get(model.mirrors, space.from_comp, nothing)
        to_m = get(model.mirrors, space.to_comp, nothing)
        if from_m !== nothing && from_m.suscept !== nothing
            arm_space_name = sname
            arm_etm = from_m
            break
        end
        if to_m !== nothing && to_m.suscept !== nothing
            arm_space_name = sname
            arm_etm = to_m
            break
        end
    end
    arm_space_name !== nothing || error("No arm cavity found (no mirror with susceptibility)")

    arm_params = derive_cavity_params(model, arm_space_name)
    om_params = derive_optomechanical_params(model, arm_etm, P_arm)

    # ── Build SLH subsystems ──

    # 1. Arm cavity: optomechanical system with optical + mechanical modes
    arm_slh = _build_arm_cavity_slh(arm_params, om_params, model)

    subsystems = SLH[arm_slh]
    subsystem_names = [arm_space_name]

    # 2. Signal extraction cavity (if present)
    sec_space_name = nothing
    for (sname, space) in model.spaces
        sname == arm_space_name && continue
        from_m = get(model.mirrors, space.from_comp, nothing)
        to_m = get(model.mirrors, space.to_comp, nothing)
        if from_m !== nothing && to_m !== nothing
            # SEC connects to the arm ITM side
            if from_m.role == "signal_extraction" || to_m.role == "signal_extraction"
                sec_space_name = sname
                break
            end
        end
    end

    sec_params_derived = nothing
    if sec_space_name !== nothing
        sec_space = model.spaces[sec_space_name]
        sec_slh, sec_params_derived = _build_sec_cavity_slh(sec_space, model)
        push!(subsystems, sec_slh)
        push!(subsystem_names, sec_space_name)
    end

    # 3. Squeezer (if present)
    sqz_slh_entry = nothing
    sqz_params_derived = nothing
    for (sqz_name, sqz) in model.squeezers
        sqz_slh_entry, sqz_params_derived = _build_squeezer_slh(sqz_name, sqz)
        push!(subsystems, sqz_slh_entry)
        push!(subsystem_names, sqz_name)
        break  # only one squeezer supported for now
    end

    # ── Concatenate all subsystems ──
    if length(subsystems) == 1
        sys = subsystems[1]
    else
        sys = concatenate(subsystems, "det")
    end

    # ── Feedbackreduce to close loops ──
    # Connect SEC output → Arm input and Arm output → SEC input
    if sec_space_name !== nothing
        # The arm cavity has ports "l_in" and "l_out"
        # The SEC has ports "in" and "out"
        # After concatenation these become prefixed with subsystem names
        arm_prefix = arm_space_name
        sec_prefix = sec_space_name

        # Close the loop: SEC out → Arm light in, Arm light out → SEC in
        sys = feedbackreduce(sys, "$(sec_prefix)_out", "$(arm_prefix)_l_in")
        sys = feedbackreduce(sys, "$(arm_prefix)_l_out", "$(sec_prefix)_in")
    end

    # Connect squeezer if present
    if sqz_slh_entry !== nothing && !isempty(model.connections)
        for conn in model.connections
            sqz = get(model.squeezers, conn.from_comp, nothing)
            sqz === nothing && continue
            # Find which subsystem the target mirror belongs to
            # The squeezer output connects to the SEC input side
            sqz_prefix = conn.from_comp
            # Determine target — typically injects into SEC
            if sec_space_name !== nothing
                sys = feedbackreduce(sys, "$(sqz_prefix)_out", "$(sec_prefix)_in")
            end
        end
    end

    # ── Build numerical parameter dict ──
    params = _build_param_dict(sys, arm_params, om_params, sec_params_derived, sqz_params_derived, model)

    return (slh=sys, params=params)
end

# ============================================================================
# SLH builders for individual subsystems
# ============================================================================

function _build_arm_cavity_slh(arm_params, om_params, model)
    hilb = FockSpace(:cavity) ⊗ FockSpace(:mirror)
    subspaces = [OpticalMode(""), MechanicalMode("")]

    a = Destroy(hilb, :a, 1)
    b = Destroy(hilb, :b, 2)

    @variables ω l κ Ω m_mirror Γ g

    # Hamiltonian: mechanical oscillator + optomechanical coupling
    H = Ω * b' * b - g * (b' + b) * (a' + a)

    # Collapse operators: cavity decay through ITM + mechanical damping
    L_ops = [κ * a, Γ * b]
    S_mat = [1 0; 0 1]

    pdict = Dict(zip(nameof.([ω, l, κ, Ω, m_mirror, Γ, g]),
                      [ω, l, κ, Ω, m_mirror, Γ, g]))
    opdict = Dict(zip(getfield.([a, b], :name), [a, b]))

    return SLH("arm", subspaces, pdict, opdict,
               ["l_in", "m_in"], ["l_out", "m_out"], S_mat, L_ops, H)
end

function _build_sec_cavity_slh(sec_space, model)
    sec_cav = cavity("sec")

    from_m = get(model.mirrors, sec_space.from_comp, nothing)
    to_m = get(model.mirrors, sec_space.to_comp, nothing)

    # Identify SRM (signal_extraction role)
    if from_m !== nothing && from_m.role == "signal_extraction"
        srm = from_m
    elseif to_m !== nothing && to_m.role == "signal_extraction"
        srm = to_m
    else
        srm = from_m  # fallback
    end

    c = c_SI.val
    κ_sec = sqrt(srm.T * c / (4 * sec_space.L_m))
    Δ_sec = sec_space.detune_rad  # detuning in rad (per round trip)

    return sec_cav, (κ_sec=κ_sec, Δ_sec=Δ_sec)
end

function _build_squeezer_slh(name, sqz::OpticalSqueezer)
    sqz_cav = squeezing_cavity(name)

    # Convert dB to squeezing parameter ε
    # For a degenerate OPA: squeeze parameter r = sqz_dB * ln(10)/20
    # The relationship between ε and the output squeezing depends on κ
    r = sqz.sqz_dB * log(10) / 20

    return sqz_cav, (r=r, sqz_angle_deg=sqz.sqz_angle_deg)
end

# ============================================================================
# Parameter dict builder
# ============================================================================

function _build_param_dict(sys, arm_params, om_params, sec_params, sqz_params, model)
    c = c_SI.val
    ℏ = ℏ_SI.val
    λ = model.lambda_m
    ω_l = 2π * c / λ

    params = Dict{Any,Any}()

    # Arm cavity parameters — find the symbolic variables in the system
    for (key, val) in sys.parameters
        skey = string(key)

        # Arm parameters (may be prefixed with "arm_" or "Arm_" after concatenation)
        if endswith(skey, "ω") || skey == "ω"
            params[val] = ω_l
        elseif endswith(skey, "l") || skey == "l"
            params[val] = om_params.L_arm
        elseif endswith(skey, "κ") || skey == "κ"
            # Could be arm κ or SEC κ — distinguish by prefix
            # κ_itm and κ_sec are already √(T·c/(4L)), ready for L = κ·a
            if startswith(skey, "sec") || startswith(skey, "SEC")
                params[val] = sec_params !== nothing ? sec_params.κ_sec : 0.0
            else
                params[val] = arm_params.κ_itm
            end
        elseif endswith(skey, "Ω") || skey == "Ω"
            params[val] = om_params.Ω_mech
        elseif skey == "m_mirror" || endswith(skey, "m_mirror")
            params[val] = om_params.mass
        elseif endswith(skey, "Γ") || skey == "Γ"
            params[val] = 0.0   # no mechanical damping
        elseif endswith(skey, "g") || skey == "g"
            params[val] = om_params.g_om
        elseif endswith(skey, "Δ") || skey == "Δ"
            if startswith(skey, "sec") || startswith(skey, "SEC")
                params[val] = sec_params !== nothing ? sec_params.Δ_sec : 0.0
            else
                params[val] = 0.0
            end
        elseif endswith(skey, "ϵ") || skey == "ϵ"
            if sqz_params !== nothing
                params[val] = sqz_params.r
            else
                params[val] = 0.0
            end
        end
    end

    return params
end
