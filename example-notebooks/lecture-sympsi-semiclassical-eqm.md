---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.11.2
---

# Lecture 7 - Symbolic quantum mechanics using SymPsi - Semiclassical equations of motion

+++

<style>
p {
    font-family: "Liberation Serif", serif;
    font-size: 12pt;
}
</style>

+++

Author: J. R. Johansson (robert@riken.jp), [http://jrjohansson.github.io](http://jrjohansson.github.io), and Eunjong Kim.

Status: Preliminary (work in progress)

This notebook is part of a series of IPython notebooks on symbolic quantum mechanics computations using 
[SymPy](http://sympy.org) and [SymPsi](http://www.github.com/jrjohansson/sympsi). SymPsi is an experimental fork and extension of the [`sympy.physics.quantum`](http://docs.sympy.org/dev/modules/physics/quantum/) module in SymPy. The latest version of this notebook is available at [http://github.com/jrjohansson/sympy-quantum-notebooks](http://github.com/jrjohansson/sympy-quantum-notebooks), and the other notebooks in this lecture series are also indexed at [http://jrjohansson.github.io](http://jrjohansson.github.com).

Requirements: A recent version of SymPy and the latest development version of SymPsi is required to execute this notebook. Instructions for how to install SymPsi is available [here](http://www.github.com/jrjohansson/sympsi).

Disclaimer: The SymPsi module is still under active development and may change in behavior without notice, and the intention is to move some of its features to [`sympy.physics.quantum`](http://docs.sympy.org/dev/modules/physics/quantum/) when they matured and have been tested. However, these notebooks will be kept up-to-date the latest versions of SymPy and SymPsi.

+++

## Setup modules

```{code-cell}
%matplotlib inline
import matplotlib.pyplot as plt
```

```{code-cell}
import numpy as np
```

```{code-cell}
from sympy import *
init_printing()
```

```{code-cell}
from sympsi import *
from sympsi.boson import *
from sympsi.pauli import *
from sympsi.operatorordering import *
from sympsi.expectation import *
from sympsi.operator import OperatorFunction
```

## Semiclassical equations of motion

+++

The dynamics of an open quantum system with a given Hamiltonian, $H$, and some interaction with an environment that acts on the system through the sytem operator $a$, and with rate $\kappa$, can often be described with a Lindblad master equation for the dynamics of the system density matrix $\rho$:

$$
\frac{d}{dt}\rho = -i[H, \rho] + \kappa \mathcal{D}[a]\rho,
$$

where the Lindblad superoperator $\mathcal{D}$ is

$$
\mathcal{D}[a]\rho = a \rho a^\dagger -\frac{1}{2}\rho a^\dagger a - \frac{1}{2}a^\dagger a \rho.
$$

One common approach to solve for the dynamics of this system is to represent the system operators and the density operator as matrices, possibly in a truncated state space, and solve the matrix-valued ODE problem numerically.

Another approach is to use the adjoint master equation for the system operators $X$:

$$
\frac{d}{dt} X =  i [H, X] + \kappa \mathcal{D}[a^\dagger]X
$$

and then solve for dynamics of the expectation values of the relevant system operators. The advantage of this method is that the ODEs are no longer matrix-valued, unlike the ODE for the density matrix. However, from the density matrix we can calculate any same-time expectation values, but with explicit ODEs for expectation values we need to select in advance which operator's expectation values we want to generate equations for. 

We can easily generate an equation for the expectation value of a specific operator by multiplying the master equation for $\rho$ from the left with an operator $X$, and then take the trace over the entire equation. Doing this we obtain:

$$
X\frac{d}{dt}\rho = -iX[H, \rho] + \kappa X\mathcal{D}[a]\rho
$$

and taking the trace:

$$
{\rm Tr}\left(X\frac{d}{dt}\rho\right) = -i{\rm Tr}\left(X[H, \rho]\right) + \kappa {\rm Tr}\left(X\mathcal{D}[a]\rho\right)
$$

using the cyclic permutation properties of traces:

$$
\frac{d}{dt}{\rm Tr}\left(X\rho\right) = -i{\rm Tr}\left([X, H]\rho\right) + \kappa {\rm Tr}\left((\mathcal{D}[a]X) \rho\right)
$$

we end up with an equation for the expectation value of the operator $X$:

$$
\frac{d}{dt}\langle X\rangle 
= 
i\langle [H, X] \rangle + \kappa \langle \mathcal{D}[a]X \rangle
$$

Note that this is a C-number equation, and therefore not as complicated to solve as the master equation for the density matrix. However, the problem with this C-number equation is that the expressions $[H, X]$
 and $\mathcal{D}[a]X$ in general will introduce dependencies on other system operators, so we obtain a system of coupled C-number equations. If this system of equations closes when a finite number of operators are included, then we can use this method to solve the dynamics of these expectation values exactly. If the system of equations do not close, which is often the case for coupled systems, then we can still use this method if we introduce some rule for truncating high-order operator expectation values (for example, by discarding high-order terms or by factoring them in expectation values of lower order). However, in this case the results are no longer exact, and is called a semi-classical equation of motion.

With SymPsi we can automatically generate semiclassical equations of motion for operators in a system described by a given Hamiltonian and a set of collapse operators that describe its coupling to an environment.

+++

## Driven harmonic oscillator

+++

Consider a driven harmonic oscillator, which interaction with an bath at some temperature that corresponds to $N_{\rm th}$ average photons. We begin by setting up symbolic variables for the problem parameters and the system operators in SymPsi: 

```{code-cell}
w, t, Nth, Ad, kappa = symbols(r"\omega, t, n_{th}, A_d, kappa", positive=True)
```

```{code-cell}
a = BosonOp("a")
rho = Operator(r"\rho")
rho_t = OperatorFunction(rho, t)
```

```{code-cell}
H = w * Dagger(a) * a + Ad * (a + Dagger(a))

Eq(Symbol("H"), H)
```

The master equation for this system can be generated using the `master_equation` function:

```{code-cell}
c_ops = [sqrt(kappa * (Nth + 1)) * a, sqrt(kappa * Nth) * Dagger(a)]
```

```{code-cell}
me = master_equation(rho_t, t, H, c_ops)

me
```

Equation for the system operators can be generated using the function `operator_master_equation`, and for the specific case of the cavity operator $a$ we obtain:

```{code-cell}
# first setup time-dependent operators
a_t = OperatorFunction(a, t)
a_to_a_t = {a: a_t, Dagger(a): Dagger(a_t)}
H_t = H.subs(a_to_a_t)
c_ops_t = [c.subs(a_to_a_t) for c in c_ops]
```

```{code-cell}
# operator master equation for a
ome_a = operator_master_equation(a_t, t, H_t, c_ops_t)

Eq(ome_a.lhs, normal_ordered_form(ome_a.rhs.doit().expand()))
```

```{code-cell}
# operator master equation for n = Dagger(a) * a
ome_n = operator_master_equation(Dagger(a_t) * a_t, t, H_t, c_ops_t)

Eq(ome_n.lhs, normal_ordered_form(ome_n.rhs.doit().expand()))
```

From these operator equations we see that the equation for $a$ depends only on the operator $a$, while the equation for $n$ depends on $n$, $a$, $a^\dagger$. So to solve the latter equation we therefore also have to generate an equation for $a^\dagger$.

+++

### System of semiclassical equations

```{code-cell}
ops, op_eqm, sc_eqm, sc_ode, ofm, oim = semi_classical_eqm(H, c_ops)
```

```{code-cell}
html_table([[Eq(Expectation(key), ofm[key]), sc_ode[key]] for key in operator_sort_by_order(sc_ode)])
```

Since this is a system of linear ODEs, we can write it on matrix form:

```{code-cell}
A_eq, A, M, b = semi_classical_eqm_matrix_form(sc_ode, t, ofm)

A_eq
```

### Steadystate

+++

We can solve for the steadystate by setting the left-hand-side of the ODE to zero, and solve the linear system of equations:

```{code-cell}
A_sol = M.LUsolve(-b)
```

The solution for the three system operators are:

```{code-cell}
A_sol[ops.index(Dagger(a)*a)]
```

```{code-cell}
A_sol[ops.index(a)]
```

```{code-cell}
A_sol[ops.index(Dagger(a))]
```

We can also solve for the steadystate directly from the ODE by settings its right-hand-side to zero, and using the SymPy `solve` function: 

```{code-cell}
solve([eq.rhs for eq in sc_ode.values()], list(ofm.values()))
```

### Solve in the ODEs

+++

For systems with a small number of dependent operators we can solve the resulting system of ODEs directly:

```{code-cell}
sols = dsolve(list(sc_ode.values())); sols
```

```{code-cell}
# hack
tt = [s for s in sols[0].rhs.free_symbols if s.name == 't'][0]
```

We also need o specify the initial conditions: Here the initial conditions are $\langle a(0) \rangle = \langle a^\dagger(0) \rangle = 2$ and $\langle a^\dagger(0)a(0) \rangle = 4$.

```{code-cell}
ics = {ofm[Dagger(a)].subs(tt, 0): 2, 
       ofm[a].subs(tt, 0): 2,
       ofm[Dagger(a)*a].subs(tt, 0): 4}; ics
```

```{code-cell}
constants = set(sum([[s for s in sol.free_symbols if (str(s)[0] == 'C')] for sol in sols], [])); constants
```

```{code-cell}
C_sols = solve([sol.subs(tt, 0).subs(ics) for sol in sols], constants); C_sols
```

```{code-cell}
sols_with_ics = [sol.subs(C_sols) for sol in sols]; sols_with_ics
```

Now let's insert numerical values for the system parameters so we can plot the solution:

```{code-cell}
values = {w: 1.0, Ad: 0.0, kappa: 0.1, Nth: 0.0}
```

```{code-cell}
sols_funcs = [sol.rhs.subs(values) for sol in sols_with_ics]; sols_funcs
```

```{code-cell}
times = np.linspace(0, 50, 500)

y_funcs = [lambdify([tt], sol_func, 'numpy') for sol_func in sols_funcs]

fig, axes = plt.subplots(len(y_funcs), 1, figsize=(12, 6))

for n, y_func in enumerate(y_funcs):
    axes[n].plot(times, np.real(y_func(times)), 'r')
    
axes[2].set_ylim(0, 5);
```

## Driven dissipative two-level system

```{code-cell}
sx, sy, sz, sm, sp = SigmaX(), SigmaY(), SigmaZ(), SigmaMinus(), SigmaPlus()
```

```{code-cell}
Omega, gamma_0, N, t = symbols("\Omega, \gamma_0, N, t", positive=True)

values = {Omega: 1.0, gamma_0: 0.5, N: 1.75}
```

```{code-cell}
H = -Omega/2 * sx
H
```

```{code-cell}
c_ops = [sqrt(gamma_0 * (N + 1)) * pauli_represent_x_y(sm), 
         sqrt(gamma_0 * N) * pauli_represent_x_y(sp)]
```

```{code-cell}
ops, op_eqm, sc_eqm, sc_ode, ofm, oim = semi_classical_eqm(H, c_ops)
```

```{code-cell}
html_table([[Eq(Expectation(key), ofm[key]), sc_ode[key]] for key in operator_sort_by_order(sc_ode)])
```

```{code-cell}
A_eq, A, M, b = semi_classical_eqm_matrix_form(sc_ode, t, ofm)

A_eq
```

### Steadystate

```{code-cell}
A_sol = M.LUsolve(-b)
```

The steadystate expectation value of $\sigma_x$:

```{code-cell}
A_sol[ops.index(sx)]
```

The steadystate expectation value of $\sigma_y$:

```{code-cell}
A_sol[ops.index(sy)]
```

The steadystate expectation value of $\sigma_z$:

```{code-cell}
A_sol[ops.index(sz)]
```

Steadystate of $\sigma_+$:

```{code-cell}
pauli_represent_x_y(sp).subs({sx: A_sol[ops.index(sx)], sy: A_sol[ops.index(sy)]})
```

Steadystate of $\sigma_-$:

```{code-cell}
pauli_represent_x_y(sm).subs({sx: A_sol[ops.index(sx)], sy: A_sol[ops.index(sy)]})
```

Alternatively we can also use the SymPy `solve` function to find the steadystate solutions:

```{code-cell}
solve([eq.rhs for eq in sc_ode.values()], list(ofm.values()))
```

At zero temperature:

```{code-cell}
solve([eq.subs(N, 0).rhs for eq in sc_ode.values()], list(ofm.values()))
```

## Versions

```{code-cell}
%reload_ext version_information

%version_information sympy, sympsi
```
