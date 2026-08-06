using SLHQuantumSystems
using SecondQuantizedAlgebra
using Symbolics
using ControlSystems
using GLMakie

mode = MechanicalMode("")

@variables Ω m Γ
b = Destroy(FockSpace(:mass),operatornames(mode)[1])
paramdict = Dict(zip(nameof.([Ω,m,Γ]),[Ω,m,Γ]))

opdict = Dict(getfield(b,:name)=>b)

slh = SLH("mass",[mode],paramdict,opdict,["in"],["out"],[1],[Γ*b],Ω*b'*b)

ss = QuantumStateSpace(slh)
qss = toquadrature(ss)

paramdict = Dict(
    Ω => 10, #10Hz resonance
    m => 4, #kg
    Γ => 1e-23               
)

numeric = substitute(qss, paramdict)


freq_hz = collect(logrange(0.1, 20_000.0, 500))
freq    = 2π .* freq_hz

# ============================================================================
# output noise spectral density
# ============================================================================

G = freqresp(numeric,freq)

sd = spectral_density(numeric, freq)

# ============================================================================
# Plotting
# ============================================================================

fig = Figure()
ax_tf = Axis(fig[1, 1];
    xscale = log10, yscale = log10,
    xlabel = "Frequency [Hz]",
    title = "transfer functions (unitless)")

for (ii,iname) in enumerate(sd.names)
    for (jj,jname) in enumerate(sd.names)
        lines!(ax_tf,freq_hz,abs.(G[ii,jj,:]),label="$iname, $jname")
    end
end
fig[1,2] = Legend(fig,ax_tf)


ax_noise = Axis(fig[2, 1];
    xscale = log10, yscale = log10,
    limits = (nothing,nothing,10e-3,10e3),
    xlabel = "Frequency [Hz]",
    title = "spectral densities (quanta per root hertz)")

for ii in sd.names
    for jj in sd.names
        lines!(ax_noise,freq_hz,abs.(sd[ii,jj]),label="$ii, $jj")
    end
end


fig[2,2] = Legend(fig,ax_noise)


fig
