# # Strain sensitivity of a single-arm Fabry–Perot cavity

using SLHQuantumSystems
using SecondQuantizedAlgebra
using Symbolics
using ControlSystems
using CairoMakie
using PhysicalConstants.CODATA2018: ReducedPlanckConstant, SpeedOfLightInVacuum

ħ = ReducedPlanckConstant.val
c = SpeedOfLightInVacuum.val

λ = 1.064e-6           # laser wavelength [m]
L_arm = 3995.0         # arm length [m]
T_ITM = 0.014          # ITM power transmittance
P_circ = 750.0e3       # circulating power [W]
m = 39.6 / 2;          # reduced mass of the Michelson [kg]

# ## Model
#
# Cavity mode `a`, and a free mirror with position `x` [m] and momentum `p`,
# with `[x, p] = i`. The mass enters as `μ = m/ħ`, and `g` is the linearized
# radiation-pressure coupling:
#
# ```math
# H = \frac{p^2}{2μ} - g\, x (a + a^†), \qquad L = [\sqrt{κ}\, a]
# ```

hilb = FockSpace(:cavity) ⊗ PhaseSpace(:mirror)
a = Destroy(hilb, :a, 1)
x = Position(hilb, :x, 2)
p = Momentum(hilb, :p, 2)

@variables κ μ g

slh = SLH(
    "arm", [OpticalMode("cavity"), MechanicalMode("mirror")],
    ["in"], ["out"],
    [1], [sqrt(κ) * a],
    p^2 / (2μ) - g * x * (a' + a)
)

qss = QuantumStateSpace(slh)

# The drift matrix. Amplitude `cavity_x` pushes the mirror (`mirror_p`); mirror
# position `mirror_x` shifts the phase `cavity_p`.

state_names(qss), qss.A

# ## Parameters

ω_l = 2π * c / λ
N = P_circ * (2 * L_arm / c) / (ħ * ω_l)   # intracavity photon number

numeric = substitute(
    qss, Dict(
        κ => T_ITM * c / (2 * L_arm),  # FWHM linewidth
        μ => m / ħ,
        g => (ω_l / L_arm) * sqrt(N),  # g₀√N / x_zpf [Aspelmeyer_Kippenberg_Marquardt_2014](@cite)
    )
)

# ## Signal and noise
#
# Strain `h` acts on the mirror as a force `F = m ω² L h / 2`, which drives the
# mirror momentum: `d(mirror_p)/dt = F / ħ`.

f = logrange(1.0, 10_000.0, 400)
ω = 2π .* f

T = transfer_matrix(numeric, ω)

signal = T["out_p", :mirror_p] .* (m .* ω .^ 2 .* L_arm ./ 2) ./ ħ

S_shot = abs2.(T["out_p", "in_p"]) ./ 2   # vacuum: ½ in each input quadrature
S_rpn = abs2.(T["out_p", "in_x"]) ./ 2
h_shot = sqrt.(S_shot) ./ abs.(signal)
h_rpn = sqrt.(S_rpn) ./ abs.(signal)
h_total = sqrt.(S_shot .+ S_rpn) ./ abs.(signal);

# ## Closed form
#
# The conventional-interferometer result of
# [Kimble_Levin_Matsko_Thorne_Vyatchanin_2001](@cite), with Michelson mass `2m`:
#
# ```math
# S_h = \frac{h_\text{SQL}^2}{2}\left(\frac{1}{\mathcal{K}} + \mathcal{K}\right), \quad
# h_\text{SQL}^2 = \frac{4ħ}{m ω^2 L^2}, \quad
# \mathcal{K} = \frac{8 ω_l P γ}{m L c\, ω^2 (ω^2 + γ^2)}
# ```
#
# where `γ = κ/2` is the cavity half-linewidth.

γ = T_ITM * c / (4 * L_arm)
h_SQL = sqrt.(4ħ ./ (m .* ω .^ 2 .* L_arm^2))
𝒦 = 8 * ω_l * P_circ * γ ./ (m * L_arm * c .* ω .^ 2 .* (ω .^ 2 .+ γ^2))
h_kimble = h_SQL .* sqrt.((1 ./ 𝒦 .+ 𝒦) ./ 2)

fig = Figure()
ax = Axis(
    fig[1, 1]; xscale = log10, yscale = log10,
    xlabel = "Frequency [Hz]", ylabel = "Strain ASD [1/√Hz]"
)
lines!(ax, f, h_total; label = "total")
lines!(ax, f, h_shot; label = "shot", linestyle = :dash)
lines!(ax, f, h_rpn; label = "radiation pressure", linestyle = :dash)
lines!(ax, f, h_kimble; label = "closed form", color = :black, linestyle = :dot)
lines!(ax, f, h_SQL; label = "SQL", color = :gray)
axislegend(ax; position = :lb)
fig
