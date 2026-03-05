using SLHQuantumSystems
using SecondQuantizedAlgebra
using Symbolics
using GLMakie

# ============================================================================
# LIGO Physical Parameters (from GWINC A+ configuration)
# Reference: pygwinc/gwinc/ifo/Aplus/ifo.yaml
# ============================================================================

# Physical constants
const c = 299792458.0              # m/s (speed of light)
const ℏ = 1.054571817e-34         # J⋅s (reduced Planck constant)

# Laser parameters
const λ = 1.064e-6                 # m (Nd:YAG wavelength)
const ω_l = 2π * c / λ             # rad/s (laser angular frequency)

# Cavity parameters
const L_arm = 3995.0               # m (arm cavity length)
const T_ITM = 0.014                # power transmittance of ITM (1.4%)

# Derived cavity parameters
const κ_cavity = T_ITM * c / (4*L_arm)   #amplitude coupling rate out of a linear cavity

# Circulating power (simplified, single cavity without power recycling)
# For full LIGO with power recycling: P_circ ≈ 750 kW
# For impedance-matched single cavity:
const P_circ = 750e3    # ≈ 750 kW

const g_optomech = 2*sqrt(ω_l*P_circ/(L_arm*c*ℏ)) 

# Test mass parameters
const m_mirror = 39.6              # kg (fused silica test mass)
const Ω_mech = 2π * 0.0         # rad/s (pendulum mode ~1 Hz)

println("=== LIGO Physical Parameters ===")
println("Arm length L = $(L_arm) m")
println("Cavity bandwidth κ/(2π) = $(round(κ_cavity, digits=1)) Hz")
println("Circulating power P_circ = $(round(P_circ/1e3, digits=1)) kW")
println("Effective coupling g = $(round(g_optomech/(2π), sigdigits=3)) Hz")
println()

# We construct the Hamiltonian using creation and annihilation operators

hilb = FockSpace(:cavity)⊗FockSpace(:mirror)

subspaces = [OpticalMode(""),MechanicalMode("")]

a = Destroy(hilb,operatornames(subspaces[1])[1],1)
b = Destroy(hilb,operatornames(subspaces[2])[1],2)

@variables ω l κ #consider replacing calls to parameternames
@variables Ω m Γ 

@variables g #coupling parameter is defined separately 

#Hamiltonian (Chen 2013 eq 2.4)
H = Ω*b'*b - g*(b'+b)*(a' + a)
L = [κ*a,Γ*b]
S = [1 0; 0 1]

params = [ω,l,κ,Ω,m,Γ,g]
pdict = Dict(zip(nameof.(params),params))

opdict = Dict(zip(getfield.([a,b],:name),[a,b]))

slh = SLH("opto",subspaces,pdict,opdict,["l_in","m_in"],["l_out","m_out"], S, L, H)
aass = StateSpace(slh)
qss = toquadrature(aass)

# ============================================================================
# Substitute physical parameter values into the model
# ============================================================================

paramdict = Dict([
    ω => 0.0,           # cavity detuning [rad/s] - on resonance
    l => L_arm,         # cavity length [m]
    κ =>  sqrt(κ_cavity),      # cavity decay rate [rad/s]
    Ω => Ω_mech,        # mechanical frequency [rad/s]
    m => m_mirror,      # mirror mass [kg]
    g => g_optomech,         # optomechanical coupling [rad/s]
    Γ => 0         # mechanical damping [rad/s]
])

println("=== Model Parameters ===")
println("Cavity detuning ω = 0 (on resonance)")
println("Cavity length l = $(L_arm) m")
println("κ/(2π) = $(round(κ_cavity/(2π), digits=1)) Hz")
println("Ω/(2π) = $(round(Ω_mech/(2π), digits=2)) Hz")
println("g/(2π) = $(round(g_optomech/(2π), sigdigits=3)) Hz")
println()

numeric = substitute(qss,paramdict)

# Frequency array in angular units (rad/s)
# SLHQuantumSystems transfer functions use angular frequency
freq_Hz = collect(logrange(0.5, 10000, 1000))  # Hz (for plotting)
freq = 2π .* freq_Hz  # rad/s (for calculations)

# Calculate transfer functions (uses angular frequency)
N = fresponse_allIO(numeric,freq)
S = fresponse_state2output(numeric, freq, 2,2)

# ============================================================================
# Calibrate to strain units
# ============================================================================
# The SLH formalism uses normalized units with ℏ=1 where:
# - Field quadratures are dimensionless
# - Mechanical position x is dimensionless
# - Vacuum noise PSD = 1/2 per quadrature
#
# Transfer functions:
# - N[i,j](ω): input field quadrature j → output field quadrature i
# - S(ω): mechanical position (normalized) → output field quadrature
#
# To convert to physical strain:
# 1. Field noise PSD = |N|² × (1/2) in ℏ=1 units
# 2. Equivalent displacement noise PSD = Field noise / |S|²
# 3. Physical displacement: multiply by x_zpf (zero-point fluctuation)
# 4. Strain: divide by L
#
# Result: Strain ASD = (x_zpf/L) × |N| / (√2 |S|)
#=
# Scaling factor from normalized to physical strain
const x_zpf = sqrt(ℏ/(2*m_mirror*Ω_mech))
strain_scale = x_zpf / L_arm
println("Strain scaling factor: $(strain_scale)")
println("  x_zpf = $(x_zpf) m")
println("  L_arm = $(L_arm) m")
println()
=#

#=
# Shot noise (from phase quadrature vacuum N[2,2])
N_shot_ASD =  abs.(N[2,2]) ./ (sqrt(2) .* abs.(S))

# Radiation pressure noise (from amplitude quadrature vacuum N[2,1])
N_rad_ASD =  abs.(N[2,1]) ./ (sqrt(2) .* abs.(S))

# Total quantum noise (add PSDs, then take sqrt for ASD)
# Note: N[2,1] and N[2,2] are independent noise sources
N_total_ASD = sqrt.(abs.(N[2,1]).^2 .+ abs.(N[2,2]).^2) ./ (sqrt(2) .* abs.(S))
=#

# Eq 6.23 part 1 from Linear Quantum Dynamical Systems
#= 
lam = 4*g_optomech^2/sqrt(κ_cavity)

N_rad_ASD = sqrt.(lam/(2*m_mirror^2*L_arm^2) .* abs2.(N[2,1])./freq.^4)

N_shot_ASD = sqrt.(1/(2*lam*L_arm^2).*abs2.(N[2,2]))

N_total_ASD = N_rad_ASD + N_shot_ASD
=#

N_total_ASD = 2 .* sqrt.(abs2.(N[2,1]).*abs2.(N[2,2])./(4*m_mirror^2*L_arm^4 .* freq.^4))


# Standard Quantum Limit for comparison
# h_SQL(ω) = √[8ℏ / (m * (ω * L)²)] where ω is angular frequency
h_SQL = sqrt.(8 * ℏ ./ (m_mirror .* (freq .* L_arm).^2))

println("=== Strain Sensitivity at 100 Hz ===")
idx_100Hz = argmin(abs.(freq .- 2π*100))
#println("Shot noise: $(round(N_shot_ASD[idx_100Hz], sigdigits=3)) 1/√Hz")
#println("Rad pressure: $(round(N_rad_ASD[idx_100Hz], sigdigits=3)) 1/√Hz")
println("Total quantum: $(round(N_total_ASD[idx_100Hz], sigdigits=3)) 1/√Hz")
println("SQL: $(round(h_SQL[idx_100Hz], sigdigits=3)) 1/√Hz")
println()

# ============================================================================
# Plotting
# ============================================================================

fig = Figure(size=(1200, 800))

# Top panel: Strain sensitivity (ASD)
ax1 = Axis(fig[1,1],
    xscale=log10, yscale=log10,
    xlabel="Frequency [Hz]",
    ylabel="Strain Sensitivity [1/√Hz]",
    title="Optomechanical Fabry Perot Quantum Noise (using aLIGO parameters)")

# Plot using freq_Hz (Hz) for x-axis, all y-values calculated with freq (rad/s)
#lines!(ax1, freq_Hz, N_shot_ASD, label="shot noise", linewidth=2)
#lines!(ax1, freq_Hz, N_rad_ASD, label="radiation pressure noise", linewidth=2)
lines!(ax1, freq_Hz, N_total_ASD, label="total quantum noise", linewidth=3)
lines!(ax1, freq_Hz, h_SQL, label="standard quantum limit",
       linewidth=2, linestyle=:dash, color=:gray)

axislegend(ax1; position=:rt)


fig
#save("optomechanical_ligo_strain.png", fig)


