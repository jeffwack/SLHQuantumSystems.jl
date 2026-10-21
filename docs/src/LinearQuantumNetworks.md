# Linear Quantum Networks

Linear quantum systems are open quantum systems that consist of
bosonic modes with quadratic Hamiltonians and linear couplings to external
fields. They are extensively studied in quantum optics, where they can describe
systems of optical resonators, beamsplitters, and quadratic optical
nonlinearity such as squeezing.

## ABCD (state space)

The 'state' of a linear system in the Heisenberg picture is specified by the
time evolved creation and annihilation operators for every mode, which we can
collect into the 'state vector' of the system.

```math
\bold{a}(t) = 
\begin{bmatrix} 
a_1(t) \\
a^\dagger_1(t) \\
\vdots \\
a_m(t) \\
a^\dagger_m(t)
\end{bmatrix}
```
This allows us to write the Heisenberg equations of motion for a linear system in vector form,

```math
\dot{\bold{a}}(t) =
i [H,\bold{a}(t)] = A\bold{a}(t),
```                      

for some matrix A.

To include inputs and outputs we introduce vectors of 'bath modes,' ``\bold{a_{in}}(t)`` and ``\bold{a_{out}}(t)``.

Now, the full equations of motion can be written as
```math
\begin{align*}
\dot{\bold{a}}(t) &= A\bold{a}(t) + B\bold{a_{in}}(t) \\
\bold{a_{out}}(t) &= C\bold{a}(t) + D\bold{a_{in}}(t)
\end{align*}
```

## Quadrature Operator Basis

We can also represent the state of the system using quadrature operators,
defined by the relations:

```math
x = \frac{1}{\sqrt{2}}(a^\dagger + a)
\quad \quad
p = \frac{i}{\sqrt{2}}(a^\dagger - a)
```

note that for a single mode,

```math
\begin{bmatrix}
x \\ p
\end{bmatrix}
= 
\frac{1}{\sqrt{2}}
\begin{pmatrix}
1 & 1\\
-i & i
\end{pmatrix}
\begin{bmatrix}
a \\ a^\dagger
\end{bmatrix},
```
giving the vector relationship

```math
\bold{x}(t) = 
\begin{bmatrix}
x_1(t) \\
p_1(t) \\
\vdots \\
x_m(t) \\
p_m(t)
\end{bmatrix}
=
\frac{1}{\sqrt{2}}
\begin{pmatrix}
1 & 1 & \cdots & 0 & 0 \\
-i & i & \cdots & 0 & 0 \\
\vdots & \vdots & \ddots & \vdots & \vdots \\
0 & 0 & \cdots & 1 & 1 \\
0 & 0 & \cdots & -i & i
\end{pmatrix}
\begin{bmatrix}
a_1(t) \\
a^\dagger_1(t) \\
\vdots \\
a_m(t) \\
a^\dagger_m(t)
\end{bmatrix}
= T_m \bold{a}(t) .
```

Where ``T_m = I_m \otimes T`` and we have defined ``T`` to be the single mode
transformation matrix.

The full equations of motion in the quadrature basis are
```math
\begin{align*}
\dot{\bold{x}}(t) &= A\bold{x}(t) + B\bold{u}(t) \\
\bold{y}(t) &= C\bold{x}(t) + D\bold{u}(t)
\end{align*}
```

## Conventions

### Decay rates

A mode that decays into a channel with energy decay rate ``\kappa`` has coupling
operator

```math
L = \sqrt{\kappa}\, a ,
```

a Fabry–Perot cavity of length ``\ell`` with input-mirror power transmission
``T \ll 1``, ``\kappa = Tc/(2\ell)``.

With these conventions, a degenerate parametric amplifier
``H = \frac{i}{2}\epsilon(a^{\dagger 2} - a^2)``, ``L = \sqrt{\kappa}\,a`` has quadrature
drift matrix ``A = \mathrm{diag}(\epsilon - \kappa/2,\, -\epsilon - \kappa/2)`` and
reaches threshold at ``\epsilon = \kappa/2``, matching Collett and Gardiner
[Collett_Gardiner_1984](@cite).

### Units and vacuum normalisation

We set ``\hbar = 1``. All rates and frequencies, and therefore all entries of the
``A`` matrix, are in rad/s.

Quadratures are dimensionless and use the same ``1/\sqrt{2}`` transform for
every mode type, optical or mechanical. Physical units enter through the
zero-point fluctuation amplitudes (`quadrature_scale`):

```math
x_\mathrm{phys} = \sqrt{2}\, x_\mathrm{zpf}\, x,
\qquad
x_\mathrm{zpf} = \sqrt{\frac{\hbar}{2 m \Omega}} .
```

The vacuum state has variance ``\frac{1}{2}`` in each quadrature, so the default
input noise in [`spectral_density`](@ref) is ``\frac{1}{2} I``.

### Frequency responses

Transfer functions are evaluated at ``s = i\omega`` with ``\omega`` in rad/s, as in
ControlSystems.jl.


## Named Inputs and Outputs

Our SLH systems have named input and output ports. When dealing with linear systems, we double the number of inputs and outputs by considering annihilation and creation
operators or both quadratures separately. When specifying a single quadrature, use the syntax [:in][1]. 
