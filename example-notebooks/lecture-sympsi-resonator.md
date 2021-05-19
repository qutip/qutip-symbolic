---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.11.2
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# Lecture 3 - Symbolic quantum mechanics using SymPsi - Resonators and cavities

+++

Author: J. R. Johansson (robert@riken.jp), [http://jrjohansson.github.io](http://jrjohansson.github.io).

Status: Preliminary (work in progress)

This notebook is part of a series of IPython notebooks on symbolic quantum mechanics computations using 
[SymPy](http://sympy.org) and [SymPsi](http://www.github.com/jrjohansson/sympsi). SymPsi is an experimental fork and extension of the [`sympy.physics.quantum`](http://docs.sympy.org/dev/modules/physics/quantum/) module in SymPy. The latest version of this notebook is available at [http://github.com/jrjohansson/sympy-quantum-notebooks](http://github.com/jrjohansson/sympy-quantum-notebooks), and the other notebooks in this lecture series are also indexed at [http://jrjohansson.github.io](http://jrjohansson.github.com).

Requirements: A recent version of SymPy and the latest development version of SymPsi is required to execute this notebook. Instructions for how to install SymPsi is available [here](http://www.github.com/jrjohansson/sympsi).

Disclaimer: The SymPsi module is still under active development and may change in behavior without notice, and the intention is to move some of its features to [`sympy.physics.quantum`](http://docs.sympy.org/dev/modules/physics/quantum/) when they matured and have been tested. However, these notebooks will be kept up-to-date the latest versions of SymPy and SymPsi.

+++

## Setup modules

```{code-cell} ipython3
from sympy import *
init_printing()
```

```{code-cell} ipython3
from sympsi import *
from sympsi.boson import *
from sympsi.operatorordering import *
```

## Introduction

+++

In this notebook we will work with cavities and resonators. A single mode of a resonator can be modelled with a quantum harmonic oscillator. Here we look at the effect of classical driving fields, transformation to different rotating frames, and various types of coupling between different modes in a resonator.

+++

## Quantum Harmonic Oscillator

+++

First we consider the Hamiltonian for a single harmonic oscillator, which describes a single mode of for example a cavity or a waveguide resonator:

$$
H = \hbar \omega_r a^\dagger a
$$

To represent this Hamiltonian in sympy we create a symbol `omega_r` and an instance of the class `BosonOp`:

```{code-cell} ipython3
omega_r = symbols("omega_r", positive=True)
Hsym = symbols("H")
a = BosonOp("a")
```

```{code-cell} ipython3
H0 = omega_r * Dagger(a) * a

Eq(Hsym, H0)
```

### Classical drive signal

+++

A classical driving field of frequency $\omega_d$ and amplitude $A$ and phase $\phi_0$ can be modelled by including an additional term in the Hamiltonian:

$$
H = \hbar \omega_r a^\dagger a + A \cos(\omega_dt  + \phi_0) (a + a^\dagger)
$$

where we have assumed that the driving field couples to the quadrature of the resonator $a + a^\dagger$, which is the canonical case. It is convenient to rewrite the $\cos$ factor

$$
H = \hbar \omega_r a^\dagger a + (Ae^{-i\omega_dt}+ A^*e^{i\omega_dt})(a + a^\dagger)
$$

where we have redefined $A \rightarrow Ae^{-i\phi_0}$. In Sympy we can represent this Hamiltonian as:

```{code-cell} ipython3
omega_d, t = symbols("omega_d, t")
A = symbols("A")
```

```{code-cell} ipython3
Hdrive = (A * exp(-I * omega_d * t) + conjugate(A) * exp(I * omega_d * t)) * (a + Dagger(a))

Hdrive
```

```{code-cell} ipython3
H = H0 + Hdrive

Eq(Hsym, H)
```

When working with Hamiltonians like this one, one common operation is to move to different rotating frames (by performing unitary transformations) where the Hamiltonian takes a simplier form.

Here we want to transform this Hamiltonian to a rotating frame in which the drive term (and all other terms) are no longer explicitly depening on time `t`. We can accomplish this by performing the unitary transformation 

$$
U = \exp(i\omega_d a^\dagger a t)
$$

```{code-cell} ipython3
U = exp(I * omega_d * t * Dagger(a) * a)

U
```

When doing a unitary basis transformation, the Hamiltonian is transformed as

$$
H \rightarrow UHU^\dagger -i U \frac{d}{dt}U^\dagger
$$

and we can carry out this transformation using the function `hamiltonian_transformation`:

```{code-cell} ipython3
H2 = hamiltonian_transformation(U, H.expand())

H2
```

Now let's introduce the detuning $\Delta = \omega_r - \omega_d$:

```{code-cell} ipython3
Delta = symbols("Delta", positive=True)
```

```{code-cell} ipython3
H3 = collect(H2, Dagger(a) * a).subs(omega_r - omega_d, Delta)

H3
```

### Rotating-wave approximation (RWA)

+++

Now we invoke the rotating-wave approximation, under which we assume that $\omega_d$ is much larger than $\Delta$, so that the we can neglect the two fast rotating terms which contains factors $\exp(\pm 2i \omega_d t)$:

```{code-cell} ipython3
H3 = drop_terms_containing(H3, [exp( 2 * I * omega_d * t),
                                exp(-2 * I * omega_d * t)])

Eq(Hsym, H3)
```

This is a time-independent hamiltonian describing the resonator in a rotating frame, where fast rotating terms have been dropped. We can no apply a displacement transformation by applying the unitary transformation:

```{code-cell} ipython3
alpha = symbols("alpha")
H = Dagger(a) * alpha - conjugate(alpha) * a
U = exp(H)

U
```

```{code-cell} ipython3
H4 = hamiltonian_transformation(U, H3)

H4 = collect(H4.expand(), [Dagger(a)*a, a, Dagger(a)])

H4
```

If we choose $\alpha = A/\Delta$ and drop c-numbers from the hamiltonian, we obtain a particularly simple form:

```{code-cell} ipython3
H5 = H4.subs(alpha, A/Delta)

H5 = drop_c_number_terms(H5)

H5
```

So a driven harmonic oscillator can be described by the hamiltonian of an undriven harmonic oscillator in a displaced frame.

+++

## Optical parametric oscillator

+++

An optical parametric oscillator (OPO) is a two-mode system with a particular nonlinear interaction. The Hamiltonian of the OPO is

$$
H = \omega_r a^\dagger a + \omega_p b^\dagger b + (\kappa {a^\dagger}^2 b + \kappa^* a^2 b^\dagger) + (A e^{i\omega_pt} + A^* e^{-i\omega_pt})(b + b^\dagger)
$$

where the operators of the two modes are$a$ and $b$, and the nonlinear interaction strength is $\kappa$, and where we have also included a classical drive field applied to mode $b$.

To model this system in SymPy we create two instances of `BosonOp`, one for each mode, and construct the Hamiltonian:

```{code-cell} ipython3
a, b = BosonOp("a"), BosonOp("b")
```

```{code-cell} ipython3
kappa, omega_p, omega_a, omega_b, t = symbols("kappa, omega_p, omega_a, omega_b, t")
Delta_a, Delta_b = symbols("Delta_a, Delta_b", positive=True)
```

```{code-cell} ipython3
Hdrive = (A * exp(-I * omega_p * t) + conjugate(A) * exp(I * omega_p * t)) * (b + Dagger(b))
H1 = omega_a * Dagger(a) * a + omega_b * Dagger(b) * b  + kappa * (a ** 2 * Dagger(b) + Dagger(a) ** 2 * b) + Hdrive

Eq(Hsym, H1)
```

We first move to a frame rotating with frequency $\omega_p$ with respect to mode $b$:

```{code-cell} ipython3
U = exp(I * omega_p * t * Dagger(b) * b)

U
```

```{code-cell} ipython3
H2 = hamiltonian_transformation(U, H1.expand(), independent=True)

H2
```

and we can perform the rotating wave approximation and drop the fast rotating terms $\exp{\left(\pm i \omega_p t\right)}$:

```{code-cell} ipython3
H3 = drop_terms_containing(H2, [exp(2 * I * omega_p * t),
                                exp(-2 * I * omega_p * t)])

H3 = collect(H3, kappa)

H3
```

Introduce a new variable for the detuning $\Delta_b = \omega_b - \omega_p$:

```{code-cell} ipython3
H4 = H3.subs(omega_b, Delta_b + omega_p).expand()

H4
```

Next we want to displace the mode $b$ so that the drive terms are eliminated:

```{code-cell} ipython3
beta = symbols("beta")
H = Dagger(b) * beta - conjugate(beta) * b
U = exp(H)

U
```

```{code-cell} ipython3
H5 = hamiltonian_transformation(U, H4.expand(), independent=True)

H5
```

```{code-cell} ipython3
H5 = collect(H5.expand(), [Dagger(a) * a, Dagger(b) * b, Dagger(a) ** 2 * b, a ** 2 * Dagger(b), b, Dagger(b)])

H5
```

The choice $\beta = A/\Delta_b$ eliminates the drive terms. After dropping c numbers we have:

```{code-cell} ipython3
H6 = H5.subs(beta, A/Delta_b)

H6 = drop_c_number_terms(H6)

H6 = collect(H6, [exp( I * omega_p * t) * kappa * a ** 2,
                  exp(-I * omega_p * t) * kappa * Dagger(a) ** 2])

H6
```

If we now assume that the dynamics of the mode $b$ is dominated by the classical drive field, we can neglect the $b$ operators in the hamiltonian (in this displaced frame, which describes deviations from the dynamics induced by the classical driving field):

```{code-cell} ipython3
H7 = drop_terms_containing(H6.expand(), [b, Dagger(b)])

H7 = collect(H7, kappa / Delta_b)

H7
```

To simpify the notation we redefine $- \kappa A / \Delta_b \rightarrow \kappa$, and assume that $A$ is real:

```{code-cell} ipython3
H = omega_a * Dagger(a) * a  + kappa * (a ** 2 * exp(I * omega_p * t) + Dagger(a) ** 2 * exp(-I * omega_p * t))

Eq(Hsym, H)
```

To eliminate the time-dependence in the interaction term we move to a rotating frame using the unitary transformation:

```{code-cell} ipython3
U = exp(I * omega_d * t * Dagger(a) * a)

U
```

```{code-cell} ipython3
H2 = hamiltonian_transformation(U, H)

H2
```

Now consider the case $\omega_p = 2 \omega_d$:

```{code-cell} ipython3
H3 = H2.subs(omega_d, omega_p/2)

H3 = collect(H3, Dagger(a) * a)

H3
```

Introduce $\Delta_a = \omega_a - \omega_p / 2$:

```{code-cell} ipython3
H4 = H3.subs(omega_a, Delta_a + omega_p/2)

H4
```

We can diagonalize this Hamiltonian by introducing the squeezing transformation:

```{code-cell} ipython3
chi = symbols("chi")
```

```{code-cell} ipython3
U = exp(chi/2 * a **2 - chi/2 * Dagger(a) ** 2)

U
```

which tranforms the operators $a$ and $a^\dagger$ according to the well-known relations:

```{code-cell} ipython3
hamiltonian_transformation(U, a)
```

```{code-cell} ipython3
hamiltonian_transformation(U, Dagger(a))
```

and in this frame the Hamiltonian takes the form:

```{code-cell} ipython3
H5 = hamiltonian_transformation(U, H4)

H5
```

```{code-cell} ipython3
H6 = normal_ordered_form(H5.expand(), independent=True)

H6 = drop_c_number_terms(H6)

H6 = collect(H6, [Dagger(a) * a, Dagger(a)**2, a**2])

H6
```

```{code-cell} ipython3
H7 = collect(H6, H6.args[1].args[0])

H7
```

```{code-cell} ipython3
# Trick to simplify the coefficients for the quantum operators
H8 = Add(*(simplify(arg.args[0]) * Mul(*(arg.args[1:])) for arg in H7.args))

H8
```

Now if we choose $\chi$ such that the coefficient of $a^2 + {a^\dagger}^2$ is zero:

```{code-cell} ipython3
chi_eq = H8.args[0].args[0]
```

```{code-cell} ipython3
Eq(chi_eq, 0)
```

```{code-cell} ipython3
chi_sol = simplify(solve(chi_eq, chi)[3])

Eq(chi, chi_sol)
```

we obtain a diagonal Hamiltonian:

```{code-cell} ipython3
H9 = H8.subs(chi_eq.args[0], -chi_eq.args[1])

H9
```

with frequency:

```{code-cell} ipython3
Eq(Symbol("omega"), H9.args[0])
```

### Versions

```{code-cell} ipython3
%reload_ext version_information

%version_information sympy, sympsi
```
