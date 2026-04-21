using SecondQuantizedAlgebra
using SLHQuantumSystems
using Symbolics
using GLMakie
using PhysicalConstants.CODATA2018: SpeedOfLightInVacuum as c
using ControlSystems

hilb = FockSpace(:cav)
a = Destroy(hilb, :a)

T = 0.00001
L = 400 #meters
κ_num = T * c.val / (4*L)

print(κ_num)

@variables ω κ

cav = SLH("cav",[1],[κ*a],ω*a'*a)

cavSS = QuantumStateSpace(cav)

quadcavSS = toquadrature(cavSS)

paramdict = Dict([ω => 100, κ =>κ_num])

numcav = substitute(quadcavSS,paramdict)

numeric = toquadrature(numcav)

# ============================================================================
# Frequency grid (angular, rad/s)
# ============================================================================

freq_hz = collect(logrange(0.1, 20_000.0, 500))
freq    = 2π .* freq_hz

# ============================================================================
# outputs
# ============================================================================

G = freqresp(numeric,freq)

sd = spectral_density(numeric, freq)

# ============================================================================
# Plotting
# ============================================================================

fig = Figure()
ax_tf = Axis(fig[1, 1];
    xscale = log10, yscale = log10,
    xlabel = "Frequency [Hz]")

for ii in 1:2
    for jj in 1:2
        lines!(ax_tf,freq_hz,abs.(G[ii,jj,:]),label="$ii, $jj")
    end
end
fig[1,2] = Legend(fig,ax_tf)


ax_noise = Axis(fig[2, 1];
    xscale = log10, yscale = log10,
    limits = (nothing,nothing,10e-2,10e2),
    xlabel = "Frequency [Hz]")

for ii in sd.names
    for jj in sd.names
        lines!(ax_noise,freq_hz,abs.(sd[ii,jj]),label="$ii, $jj")
    end
end


fig[2,2] = Legend(fig,ax_noise)
fig
