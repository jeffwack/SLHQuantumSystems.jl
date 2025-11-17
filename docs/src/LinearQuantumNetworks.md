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

### Optical Modes
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

### Mechanical Modes
The operators ``x`` and ``p`` have different units for mechanical systems compared to optical quadratures, and thus the transformation from annihilation and creation operators
to quadratures takes a slightly different form.

```math
\begin{bmatrix}
x \\ p
\end{bmatrix}
=
\frac{1}{\sqrt{2m\omega}}
\begin{pmatrix}
1 & 1 \\
-i m \omega & i m \omega
\end{pmatrix}
\begin{bmatrix}
a \\ a^\dagger
\end{bmatrix}
``` 

## Named Inputs and Outputs

Our SLH systems have named input and output ports. When dealing with linear systems, we double the number of inputs and outputs by considering annihilation and creation
operators or both quadratures separately. When specifying a single quadrature, use the syntax [:in][1]. 
