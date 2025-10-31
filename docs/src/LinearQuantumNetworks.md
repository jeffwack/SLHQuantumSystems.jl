# Linear Quantum Networks

Linear quantum systems are open quantum systems that consist of harmonic
bosonic modes with quadratic Hamiltonians and linear couplings to external
fields. They are extensively studied in quantum optics, where they can describe
systems of optical resonators, beamsplitters, and quadratic optical
nonlinearity such as two-mode squeezing.

## ABCD (state space)

The 'state' of a linear system in the Heisenberg picture is specified by the
time evolved creation and annihilation operators for every mode, which we can
collect into the 'state vector' of the system.

```math
\vec{\bold{a}}(t) = 
\begin{bmatrix} 
a_1(t) \\
a^\dagger_1(t) \\
\vdots \\
a_m(t) \\
a^\dagger_m(t)
\end{bmatrix}
```                                    

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
\vec{\bold{x}}(t) = 
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
i & -i & \cdots & 0 & 0 \\
\vdots & \vdots & \ddots & \vdots & \vdots \\
0 & 0 & \cdots & 1 & 1 \\
0 & 0 & \cdots & i & -i
\end{pmatrix}
\begin{bmatrix}
a_1(t) \\
a^\dagger_1(t) \\
\vdots \\
a_m(t) \\
a^\dagger_m(t)
\end{bmatrix}
= T_m \vec{\bold{a}}(t) .

```

Where ``T_m = I_m \otimes T`` and we have defined ``T`` to be the single mode
transformation matrix.


