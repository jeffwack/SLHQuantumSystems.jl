using SecondQuantizedAlgebra
using SLHQuantumSystems
using Symbolics
using GLMakie
using PhysicalConstants.CODATA2018: SpeedOfLightInVacuum as c
using ControlSystems

# ## Symbolic construction

hilb = FockSpace(:cav)
a = Destroy(hilb, :a)

@variables ω κ

cav = SLH("cav", [1], [κ * a], ω * a' * a)

cavSS = QuantumStateSpace(cav)

quadcavSS = toquadrature(cavSS)

# ## Calculation of coupling rate

T = 0.00001
L = 400; #meters

# The round-trip travel time of the cavity is 2L/c.
# T is the probability of photon leaving the cavity after one round-trip,
# thus the amplitude coupling rate is

κ_num = T * c.val / (4 * L)

# ## Numeric substitution

paramdict = Dict([ω => 1 * 2π, κ => κ_num])

numcav = substitute(quadcavSS, paramdict)

numeric = toquadrature(numcav)

# ## Frequency grid

freq_hz = collect(logrange(0.1, 100, 500))
freq = 2π .* freq_hz;

# ## Outputs

G = freqresp(numeric, freq) #this is a ControlSystems call

# ## Plotting

fig = Figure()
ax_tf = Axis(
    fig[1, 1];
    xscale = log10, yscale = log10,
    xlabel = "Frequency [Hz]"
)

for ii in 1:2
    for jj in 1:2
        lines!(ax_tf, freq_hz, abs.(G[ii, jj, :]), label = L"q_{%$jj} \rightarrow q_{%$ii}")
    end
end
fig[1, 2] = Legend(fig, ax_tf)

fig
