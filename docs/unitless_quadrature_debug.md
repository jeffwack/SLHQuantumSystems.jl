# Unitless Quadrature Implementation Debug Notes

## What was implemented (Steps 1, 2, 3, 5 of the plan)

- `quadrature_parameter_names(MechanicalMode)` → returns `Symbol[]`
- `quadrature_transform(MechanicalMode, params)` → uses same `1/√2` unitary as OpticalMode
- Removed `quadratureblocks` shim; added `quadrature_transform(subsys::Subspace)` default
- Updated `toquadrature` in `abcd.jl` to call `quadrature_transform` directly
- Updated `calibrated_strain.jl` signal factor: `4·Ω·x_zpf` → `2√2·Ω·x_zpf`
- Updated docstrings

## Observed failure

The new code's sensitivity curve beats the SQL, which should not happen for vacuum input.

---

## Root cause analysis

### The plan's mathematical error

The plan claims the old mechanical momentum state is:
```
p_old = 0.5·i·m·Ω·(b†-b) = p_phys / (2·p_zpf)    ← WRONG
```

The actual relation is:
```
p_phys = i·p_zpf·(b†-b) = i·m·Ω·x_zpf·(b†-b)
p_old  = i·m·Ω/2·(b†-b) = p_phys / (2·x_zpf)      ← correct (has units kg/s)
```

The plan confused `x_zpf` with `p_zpf`. The old `p_old` is NOT dimensionless — it has units of kg/s, unlike the new `r = p_phys/(√2·p_zpf)` which is dimensionless.

### The correct similarity transform

The transform from old to new mechanical quadrature coordinates:
```
q_new = √2 · q_old
r     = (√2 / (m·Ω)) · p_old
```

So the full 4×4 similarity: `S_full = diag(1, 1, √2, √2/(m·Ω))`.

The plan assumed `S = diag(1, 1, √2, √2)`, which is wrong by a factor of `1/(m·Ω)` in the 4th diagonal.

### Consequence for T_new vs T_old

From the similarity transform:
```
T_new = C_new · (iω - A_new)^{-1} · e₄
      = C_old · (iω - A_old)^{-1} · S_full^{-1} · e₄
      = (m·Ω/√2) · T_old
```

The plan claimed `T_new = (1/√2) · T_old`. The correct ratio is `m·Ω/√2 ≈ 176` for LIGO.

### What the new signal formula gives

New code:
```
H_signal_new = T_new · ω²·L / (2√2·Ω·x_zpf)
             = (m·Ω/√2) · T_old · ω²·L / (2√2·Ω·x_zpf)
             = m · T_old · ω²·L / (4·x_zpf)
```

### What the CORRECT old signal formula should be

From Heisenberg EOM with physical ℏ, the GW force drives state 4:
```
ṗ_old_gw = F_phys / (2·x_zpf) = m·(L/2)·(-ω²)·h / (2·x_zpf)
```
So `factor_old_correct = m·ω²·L / (4·x_zpf)`, giving:
```
H_signal_old_correct = T_old · m·ω²·L / (4·x_zpf)
```

These match: **H_signal_new = H_signal_old_correct** ✓

So the new formula is internally consistent and correct from first principles.

### The old code's formula was wrong

The old code used `ω²·L / (4·Ω·x_zpf)`, which differs from `m·ω²·L / (4·x_zpf)` by a factor of `m·Ω ≈ 249` for LIGO. The old code **underestimated** the GW signal by m·Ω, making the sensitivity appear 249× worse than reality.

---

## The key unresolved tension

- The new code is **first-principles correct**
- The old code had a wrong signal factor (missing m·Ω), making h_ASD appear m·Ω times worse than truth
- If the new code (correct physics) shows h_ASD < h_SQL → this implies the physical system has sub-SQL sensitivity
- But vacuum input to a standard optomechanical cavity should respect the SQL

### Possible resolutions

1. **The correct physics IS near/at SQL** — the old code's factor error made things appear 249× worse, and the corrected sensitivity is near (but not below) SQL. The "beating SQL" claim might be a misreading of the plot, or a small numerical effect near the pendulum resonance.

2. **The SQL comparison formula is wrong** — the code uses `h_SQL = sqrt(8ℏ/(m·ω²·L²))` for a free mass. For a harmonic oscillator pendulum at Ω=1 Hz, the SQL near mechanical resonance can differ. Perhaps the curve dips below h_SQL_freemass near Ω but respects it in the signal band.

3. **There is still a factor error in the new formula** — perhaps in the derivation of the GW injection force (sign convention, factor of 2, or L vs L/2).

4. **The spectral_density noise is also affected** — `S_out = G·S_in·G†` — if G (the input-to-output transfer function from B-matrix) changes between old and new (it should not for the optical rows, since B_mech ≈ 0 and B_opt is unchanged), the noise would also change. This needs verification.

---

## What to check next

1. Numerically verify that `noise_ASD` (`asd(sd, "l_out_p")`) is unchanged between old and new code (it should be, since the optical output B-matrix rows are unaffected by the mechanical quadrature change).

2. Numerically verify the ratio `T_new / T_old ≈ m·Ω/√2 ≈ 176` by running both codes and comparing `fresponse_state2output`.

3. Determine how far below SQL the new code goes, and whether this occurs at or near the mechanical resonance (1 Hz) vs in the signal band (10–10000 Hz).

4. Check whether the SQL formula `h_SQL = sqrt(8ℏ/(m·ω²·L²))` is the right comparison for this model (single FP cavity with one free mirror).

5. Consider whether the old code's factor `ω²·L/(4·Ω·x_zpf)` was empirically tuned (not first-principles derived), and if so what it was tuned to match.

---

## Files changed so far

- `src/subspace.jl`: quadrature_parameter_names, quadrature_transform, removed quadratureblocks
- `src/abcd.jl`: toquadrature call sites, QuadratureBasis docstring
- `examples/calibrated_strain.jl`: signal factor denominator 4→2√2, updated comment
