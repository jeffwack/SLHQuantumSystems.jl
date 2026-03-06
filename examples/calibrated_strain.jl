using SLHQuantumSystems
using SecondQuantizedAlgebra
using Symbolics
using ControlSystems
using Plots
using LinearAlgebra

# ============================================================================
# LIGO Physical Parameters (from GWINC A+ configuration)
# ============================================================================

const c_phys   = 299792458.0           # m/s
const ℏ        = 1.054571817e-34       # J·s
const λ        = 1.064e-6              # m  (Nd:YAG)
const ω_l      = 2π * c_phys / λ       # rad/s

const L_arm    = 3995.0                # m  (arm cavity length)
const T_ITM    = 0.014                 # power transmittance of ITM
const κ_cavity = T_ITM * c_phys / (4*L_arm)   # amplitude decay rate [rad/s]

const P_circ   = 750e3                 # W  (circulating power)
const m_mirror = 39.6                  # kg
const Ω_mech   = 2π * 1.0             # rad/s — physical pendulum (~1 Hz);
                                       #   free-mass approx valid for f >> 1 Hz

# Linearized optomechanical coupling G = g_OM * √N̄
#   g_OM = (ω_cav/L) · x_zpf          single-photon coupling [rad/s]
#   x_zpf = √(ℏ/2mΩ)                  zero-point motion [m]
#   N̄ = P_circ / (ℏ ω_cav)            intracavity photon number
const x_zpf       = sqrt(ℏ / (2 * m_mirror * Ω_mech))
const g_OM        = (ω_l / L_arm) * x_zpf          # [rad/s]
const N_bar       = P_circ / (ℏ * ω_l)             # [dimensionless]
const g_optomech  = g_OM * sqrt(N_bar)              # linearized coupling [rad/s]

println("=== LIGO Parameters ===")
println("κ/(2π)    = $(round(κ_cavity/(2π),   digits=1)) Hz")
println("Ω/(2π)    = $(round(Ω_mech/(2π),     digits=2)) Hz  (pendulum, below signal band)")
println("g_OM/(2π) = $(round(g_OM/(2π),       sigdigits=3)) Hz  (single-photon coupling)")
println("N̄         = $(round(N_bar,           sigdigits=3)) photons")
println("G/(2π)    = $(round(g_optomech/(2π), sigdigits=3)) Hz  (linearized coupling)")
println("G/κ       = $(round(g_optomech/κ_cavity, sigdigits=3))  (strong-coupling regime)")

# ============================================================================
# SLH model: optomechanical cavity (Chen 2013 eq 2.4)
# ============================================================================

hilb      = FockSpace(:cavity) ⊗ FockSpace(:mirror)
subspaces = [OpticalMode(""), MechanicalMode("")]

a = Destroy(hilb, operatornames(subspaces[1])[1], 1)
b = Destroy(hilb, operatornames(subspaces[2])[1], 2)

@variables ω l κ Ω m Γ g

H     = Ω*b'*b - g*(b'+b)*(a' + a)
L_ops = [κ*a, Γ*b]
S_mat = [1 0; 0 1]

pdict  = Dict(zip(nameof.([ω,l,κ,Ω,m,Γ,g]), [ω,l,κ,Ω,m,Γ,g]))
opdict = Dict(zip(getfield.([a,b], :name), [a,b]))

slh = SLH("opto", subspaces, pdict, opdict,
          ["l_in","m_in"], ["l_out","m_out"], S_mat, L_ops, H)

qss = toquadrature(QuantumStateSpace(slh))

# ============================================================================
# Numerical substitution
# ============================================================================

paramdict = Dict(
    ω => 0.0,              # on resonance
    l => L_arm,
    κ => sqrt(κ_cavity),   # coupling amplitude [√(rad/s)]; κ² is decay rate
    Ω => Ω_mech,
    m => m_mirror,
    g => g_optomech,
    Γ => 0.0               # no mechanical damping → marginally stable (poles at ±iΩ)
)

numeric = substitute(qss, paramdict)

poles_numeric = eigvals(Matrix{ComplexF64}(numeric.A))
println("\n=== System Poles ===")
for p in sort(poles_numeric, by=real)
    println("  $(round(real(p), sigdigits=4)) + $(round(imag(p), sigdigits=4))i  rad/s")
end

# ============================================================================
# Frequency grid (angular, rad/s)
# ============================================================================

freq_Hz = collect(logrange(0.1, 20_000.0, 500))
freq    = 2π .* freq_Hz

# ============================================================================
# Transfer function matrix G(Ω)
#
# freqresp returns shape (nout, nin, nω):  G[:,:,k] is the 4×4 complex
# transfer matrix at angular frequency freq[k].
#
# Quadrature port ordering (interlaced, after toquadrature):
#   index 1: l_out/in amplitude (x)
#   index 2: l_out/in phase    (p)  ← homodyne readout
#   index 3: m_out/in position (x)
#   index 4: m_out/in momentum (p)
# ============================================================================

G = freqresp(numeric, freq)    # (4, 4, nω) ComplexF64

# ============================================================================
# Output noise spectral density
#
#   Sout(Ω) = G(Ω) Sin G†(Ω)
#
# Vacuum input in ℏ=1 units: each quadrature carries noise 1/2, so
#   Sin = (1/2) I
#
# For a lossless system G is symplectic (not necessarily unitary), so
# individual Sout diagonal entries may exceed 1/2 — noise is redistributed
# among output quadratures by the optomechanical interaction.
# ============================================================================

nin        = size(G, 2)
Sin_vacuum = Matrix(I, nin, nin) / 2

Sout = [G[:,:,k] * Sin_vacuum * G[:,:,k]' for k in axes(G, 3)]

# ============================================================================
# Sanity check: shot noise of vacuum = 1/2
#
# With g=0 the optical and mechanical ports fully decouple.  Each subsystem
# is lossless and the transfer matrix is unitary at every frequency, so every
# diagonal element of Sout must equal exactly 1/2.
# ============================================================================

let
    p0    = merge(paramdict, Dict(g => 0.0))
    G0    = freqresp(substitute(qss, p0), freq)
    Sout0 = [G0[:,:,k] * Sin_vacuum * G0[:,:,k]' for k in axes(G0, 3)]

    max_err = maximum(
        abs(real(Sout0[k][i,i]) - 0.5)
        for k in eachindex(Sout0), i in 1:size(Sout0[1], 1)
    )

    @assert max_err < 1e-8 "Vacuum shot noise test FAILED  (max_err = $max_err)"
    println("\n✓ Vacuum shot noise test passed  (g=0, max diag error = $(round(max_err, sigdigits=2)))")
end

# ============================================================================
# Optical output noise
#
# The optical phase quadrature (index 2) is the homodyne readout channel.
# With g≠0 the radiation pressure coupling can amplify noise in this
# quadrature above vacuum level, reflecting the optomechanical back-action.
# ============================================================================

S_opt_amp   = real.([Sout[k][1,1] for k in eachindex(Sout)])   # amplitude PSD
S_opt_phase = real.([Sout[k][2,2] for k in eachindex(Sout)])   # phase PSD (homodyne)

ASD_amp   = sqrt.(abs.(S_opt_amp))
ASD_phase = sqrt.(abs.(S_opt_phase))

idx_100 = argmin(abs.(freq_Hz .- 100.0))
println("\n=== Optical output noise at 100 Hz ===")
println("  amplitude quad Sout = $(round(S_opt_amp[idx_100],   sigdigits=4))")
println("  phase quad     Sout = $(round(S_opt_phase[idx_100], sigdigits=4))")
println("  vacuum reference    = 0.5")

# ============================================================================
# Strain sensitivity
#
# A gravitational wave of strain h(ω) causes a differential arm-length change
# δL = h·L_arm/2, which acts on the mirror as a tidal acceleration
#
#   ẍ_phys = (L_arm/2)·ḧ   →   δx_phys = (L_arm/2)·h   (free-mass limit)
#
# This is equivalent to a driving force on the mechanical momentum (state 4)
#
#   F_phys(ω) = m·(L_arm/2)·(−ω²)·h
#
# In SLH quadrature units (p_zpf = m·Ω·x_zpf):
#
#   F_SLH(ω) = F_phys / (2·p_zpf)
#            = −ω²·L_arm / (4·Ω_mech·x_zpf) · h
#
# The signal transfer function from strain h to optical phase output (index 2)
# is therefore
#
#   H_signal(ω) = T_p(ω) · ω²·L_arm / (4·Ω_mech·x_zpf)
#
# where T_p = fresponse_state2output(·, 4, 2) is the transfer from a
# unit drive on the momentum state to the homodyne output.
# ============================================================================

T_mech2opt = fresponse_state2output(numeric, freq, 4, 2)   # complex, length nω

# Signal transfer function (dimensionless: output response per unit strain h)
H_signal = T_mech2opt .* (freq.^2 .* L_arm ./ (4 .* Ω_mech .* x_zpf))

# Strain-referred noise ASD  [1/√Hz, SI via x_zpf]
h_ASD = sqrt.(abs.(S_opt_phase)) ./ abs.(H_signal)

# Standard Quantum Limit for a free mass  h_SQL = √(8ℏ / (m ω² L²))
h_SQL = sqrt.(8 .* ℏ ./ (m_mirror .* freq.^2 .* L_arm^2))

idx_100 = argmin(abs.(freq_Hz .- 100.0))
println("\n=== Strain sensitivity at 100 Hz ===")
println("  h_ASD = $(round(h_ASD[idx_100], sigdigits=3)) /√Hz")
println("  h_SQL = $(round(h_SQL[idx_100], sigdigits=3)) /√Hz")
println("  ratio h/h_SQL = $(round(h_ASD[idx_100]/h_SQL[idx_100], sigdigits=3))")

# ============================================================================
# Plotting
# ============================================================================

plot(freq_Hz, ASD_amp,
    xscale    = :log10,
    yscale    = :log10,
    xlabel    = "Frequency [Hz]",
    ylabel    = "Output noise ASD  [1/√Hz,  ℏ=1 units]",
    title     = "Optomechanical cavity — vacuum output noise",
    label     = "l_out amplitude quad",
    linewidth = 2)

plot!(freq_Hz, ASD_phase,
    label     = "l_out phase quad (homodyne)",
    linewidth = 2)

hline!([sqrt(0.5)],
    label     = "vacuum reference  √½ ≈ 0.707",
    linestyle = :dash,
    color     = :gray)

# ============================================================================
# LIGO sensitivity bucket
# ============================================================================

plot(freq_Hz, h_ASD,
    xscale    = :log10,
    yscale    = :log10,
    xlabel    = "Frequency [Hz]",
    ylabel    = "Strain sensitivity  [1/√Hz]",
    title     = "Optomechanical cavity — strain sensitivity",
    label     = "quantum noise (shot + radiation pressure)",
    linewidth = 2,
    color     = :blue)

plot!(freq_Hz, h_SQL,
    label     = "SQL (free mass)",
    linestyle = :dash,
    linewidth = 2,
    color     = :red)
