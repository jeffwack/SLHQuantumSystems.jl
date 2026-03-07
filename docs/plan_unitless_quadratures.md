# Plan: Unitless Mechanical Quadratures

## Motivation

`MechanicalMode.quadrature_transform` currently embeds `m` and `Ω`, producing a state
vector `(x_phys/(2·p_zpf), p_phys/(2·p_zpf))` with inconsistent scaling relative to the
optical convention. The A matrix therefore mixes rad/s-scale optical entries with entries
that carry `m·Ω` factors. `quadrature_scale` and the explicit `x_zpf` call in
`calibrated_strain.jl` compensate after the fact. This plan eliminates the asymmetry: all
modes use the same `1/√2` unitary transform, the A matrix is uniformly in rad/s, and SI
conversion happens exclusively at the output stage.

---

## Mathematical basis for the formula change

The current and new mechanical momentum quadrature (state 4) relate as follows.

- **Current**: `p_old = 0.5·i·m·Ω·(b†-b) = p_phys/(2·p_zpf)` (dimensionless)
- **New**: `r = i·(b†-b)/√2 = p_phys/(√2·p_zpf)` (same √2 convention as optical)
- Ratio: `r = √2·p_old`

The same scaling applies to the position quadrature (state 3): `q = (b+b†)/√2 = √2·x_old`.
The full state similarity transform is therefore `S = diag(1, 1, √2, √2)`.

**Transfer function scaling**: with `C_old` having zero entries in the mechanical columns
(outputs are optical fields),

```
T_new = C · (iω - A_new)⁻¹ · e₄
      = C · S · (iω - A_old)⁻¹ · S⁻¹ · e₄
      = (1/√2) · T_old
```

**New GW injection factor**: state 4 is now `r = p_phys/(√2·p_zpf)`, so a physical force
drives `ṙ = F_phys/(√2·p_zpf)`. With `F_phys = m·(L_arm/2)·ω²·h`:

```
F_SLH = ω²·L_arm·h / (2√2·Ω·x_zpf)
```

**Invariance check**: the two effects cancel exactly,

```
H_signal_new = T_new · factor_new
             = (1/√2)·T_old · ω²·L / (2√2·Ω·x_zpf)
             = T_old · ω²·L / (4·Ω·x_zpf)
             = H_signal_old   ✓
```

so the physical sensitivity curve is unchanged.

---

## Step 1 — `src/subspace.jl`: simplify `MechanicalMode` transform

### `quadrature_parameter_names`

Change to return `Symbol[]`, matching `OpticalMode`. The transform no longer needs `m` or `Ω`.

```julia
# before
quadrature_parameter_names(subsys::MechanicalMode) = parameternames(subsys)[1:2]

# after
quadrature_parameter_names(subsys::MechanicalMode) = Symbol[]
```

### `quadrature_transform`

Replace the dimensionful matrix with the same `1/√2` unitary as `OpticalMode`.

```julia
# before
function quadrature_transform(subsys::MechanicalMode, params::Dict)
    m = params[param_key(subsys, :m)]
    w = params[param_key(subsys, :Ω)]
    left  = [0.5 0.5; -0.5im*m*w 0.5im*m*w]
    right = [1 im/(m*w); 1 -im/(m*w)]
    return (left, right)
end

# after — identical to OpticalMode
function quadrature_transform(subsys::MechanicalMode, params::Dict)
    c = 1/sqrt(Num(2))
    left  = c*[1 1; -im im]
    right = c*[1 im; 1 -im]
    return (left, right)
end
```

The resulting mechanical state coordinates are `q = (b+b†)/√2` and `r = i(b†-b)/√2` —
dimensionless, normalized consistently with optical X and P. All A-matrix entries are now
in rad/s.

### `quadrature_scale`

Keep the return value `[x_zpf, p_zpf]` but update the docstring. Its role shifts from a
post-hoc normalization correction to the definitive SI conversion:

```
x_phys = √2 · x_zpf · q
p_phys = √2 · p_zpf · r
```

---

## Step 2 — `src/subspace.jl` and `src/abcd.jl`: remove `quadratureblocks` shim

With all modes now using a parameter-free `1/√2` transform, the `quadratureblocks`
compatibility shim in `subspace.jl` is dead code:

```julia
# remove entirely
function quadratureblocks(sys, subsys::Subspace)
    relevant = Dict(k => sys.parameters[k] for k in quadrature_parameter_names(subsys))
    return quadrature_transform(subsys, relevant)
end
```

The shim's only job was to build a `relevant` params dict from `sys.parameters` for the
mechanical mode. After Step 1, `quadrature_parameter_names` returns `[]` for every mode,
so `relevant` is always empty and the shim reduces to `quadrature_transform(subsys, Dict())`.

Update the two call sites in `toquadrature` (`abcd.jl`) to call `quadrature_transform`
directly:

```julia
# before
blockpairs   = [quadratureblocks(sys, mode)        for mode in sys.subspaces]
blockpairsIO = [quadratureblocks(sys, GenericMode("")) for ii in 1:n_ports]

# after
blockpairs   = [quadrature_transform(mode,          Dict()) for mode in sys.subspaces]
blockpairsIO = [quadrature_transform(GenericMode(""), Dict()) for ii in 1:n_ports]
```

Optionally, make the `params` argument default to `Dict()` across all `quadrature_transform`
methods so call sites can omit it entirely:

```julia
quadrature_transform(subsys::Subspace) = quadrature_transform(subsys, Dict())
```

---

## Step 3 — `examples/calibrated_strain.jl`: update signal injection formula

```julia
# before
H_signal = T_mech2opt .* (freq.^2 .* l_arm ./ (4 .* Ω_mech .* x_zpf))

# after
H_signal = T_mech2opt .* (freq.^2 .* l_arm ./ (2*sqrt(2) .* Ω_mech .* x_zpf))
```

Update the comment block to reflect the new derivation:

> In SLH quadrature units, the mechanical momentum state is `r = p_phys/(√2·p_zpf)`,
> so a physical force drives `ṙ = F_phys/(√2·p_zpf)`. With `F_phys = m·(L_arm/2)·(-ω²)·h`:
>
> `F_SLH = -ω²·L_arm·h / (2√2·Ω·x_zpf)`
>
> The signal transfer function is `H_signal(ω) = T_p(ω) · ω²·L_arm / (2√2·Ω·x_zpf)`.

The `T_mech2opt` call (`fresponse_state2output`) requires no change — it automatically
reflects the new A matrix.

---

## Step 4 — Tests: verify and add coverage

The existing `runtests.jl` tests do not cover the quadrature transform or A-matrix entries
numerically, so no existing tests break. Add a new `@testset "optomechanical quadrature"`
that:

1. Builds the optomechanical SLH from `calibrated_strain.jl` (symbolic form).
2. Calls `toquadrature` and substitutes the LIGO-scale numeric values.
3. Asserts all eigenvalues of A are finite and have units consistent with rad/s scale.
4. Ports the vacuum shot-noise sanity check from `calibrated_strain.jl`: with `g = 0`,
   the diagonal of the output spectral density must equal `0.5` to within `1e-8`.

---

## Step 5 — Update docstrings

| Location | Change |
|---|---|
| `QuadratureBasis` struct | Note that all modes use the `1/√2` unitary; mechanical coordinates are `q = (b+b†)/√2`, `r = i(b†-b)/√2` |
| `quadrature_scale(MechanicalMode, ...)` | Clarify `[x_zpf, p_zpf]` are SI conversion factors; physical values recovered as `x_phys = √2·x_zpf·q`, `p_phys = √2·p_zpf·r` |
| `MechanicalMode` struct | Describe the dimensionless quadrature convention |

---

## What does NOT change

- `zpf_length`, `zpf_momentum`, `g0_coupling`, `g_coupling` — operate on physical
  parameters only, unaffected.
- `quadrature_scale` return values — `[x_zpf, p_zpf]` remain correct.
- `toquadrature` logic in `abcd.jl` — the block-diagonal assembly is unchanged; only the
  call site is updated from `quadratureblocks` to `quadrature_transform` directly.
- The `x_zpf` extraction line in `calibrated_strain.jl`:
  `zpf_length(mech_sub, numeric.parameters)` — unchanged, this is exactly the
  "convert at the end" pattern the new design formalises.
