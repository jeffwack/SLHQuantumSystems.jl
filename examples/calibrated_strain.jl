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

# The mirror is a free mass. `Ω_m` and `Γ_m` only regulate the model;
# the result does not depend on them.

Ω_m = 2π * 5.0e-6      # mechanical resonance [rad/s]
Γ_m = 1.0e-23;         # mechanical port coupling [rad/s]

# ## Model
#
# Cavity mode `a`, mirror mode `b`, linearized radiation-pressure coupling `g`:
#
# ```math
# H = Ω b^†b - g (b + b^†)(a + a^†), \qquad L = [\sqrt{κ}\, a, \sqrt{Γ}\, b]
# ```

hilb = FockSpace(:cavity) ⊗ FockSpace(:mirror)
a = Destroy(hilb, :a, 1)
b = Destroy(hilb, :b, 2)

@variables κ Ω Γ g

slh = SLH(
    "arm", [OpticalMode(""), MechanicalMode("")],
    ["l_in", "m_in"], ["l_out", "m_out"],
    [1 0; 0 1], [sqrt(κ) * a, sqrt(Γ) * b],
    Ω * b' * b - g * (b' + b) * (a' + a)
)

qss = toquadrature(QuantumStateSpace(slh))

# The drift matrix, state order `(x_a, p_a, x_b, p_b)`. Amplitude `x_a` pushes
# the mirror (`p_b`); mirror position `x_b` shifts the phase `p_a`.

qss.A

# ## Parameters

ω_l = 2π * c / λ
x_zpf = sqrt(ħ / (2 * m * Ω_m))
N = P_circ * (2 * L_arm / c) / (ħ * ω_l)   # intracavity photon number

numeric = substitute(qss, Dict(
    κ => T_ITM * c / (2 * L_arm),          # FWHM linewidth
    Ω => Ω_m,
    Γ => Γ_m,
    g => (ω_l / L_arm) * x_zpf * sqrt(N),  # g₀√N [Aspelmeyer_Kippenberg_Marquardt_2014](@cite)
))

# ## Signal and noise
#
# Strain `h` acts on the mirror as a force `F = m ω² L h / 2`. A force enters
# through the momentum quadrature of `m_in` as `p_in = √2 x_zpf F / (ħ √Γ)`.

f = logrange(1.0, 10_000.0, 400)
ω = 2π .* f

G = freqresp(numeric, ω)
lx, lp, mp = 1, 2, 4   # quadrature channels: (l_x, l_p, m_x, m_p)

signal = G[lp, mp, :] .* sqrt(2) .* x_zpf .* (m .* ω .^ 2 .* L_arm ./ 2) ./ (ħ * sqrt(Γ_m))

S_out = spectral_density(numeric, ω)["l_out_p", "l_out_p"]
h_total = sqrt.(real.(S_out)) ./ abs.(signal)
h_shot = sqrt.(abs2.(G[lp, lp, :]) ./ 2) ./ abs.(signal)
h_rpn = sqrt.(abs2.(G[lp, lx, :]) ./ 2) ./ abs.(signal);

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
ax = Axis(fig[1, 1]; xscale = log10, yscale = log10,
    xlabel = "Frequency [Hz]", ylabel = "Strain ASD [1/√Hz]")
lines!(ax, f, h_total; label = "total")
lines!(ax, f, h_shot; label = "shot", linestyle = :dash)
lines!(ax, f, h_rpn; label = "radiation pressure", linestyle = :dash)
lines!(ax, f, h_kimble; label = "closed form", color = :black, linestyle = :dot)
lines!(ax, f, h_SQL; label = "SQL", color = :gray)
axislegend(ax; position = :lb)
fig
