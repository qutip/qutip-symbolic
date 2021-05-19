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

# Lecture 4 - Symbolic quantum mechanics using SymPsi - Atom and cavity

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

```{code-cell} ipython3
from sympy import collect, symbols, simplify, Add, Mul, Eq, I, exp, powsimp, init_printing
from sympy.physics.quantum import Dagger
from sympy.physics.quantum.boson import BosonOp
init_printing()
```

```{code-cell} ipython3
from qutip_symbolic.compat.pauli import SigmaX, SigmaY, SigmaZ, SigmaMinus, SigmaPlus
```

```{code-cell} ipython3
from qutip_symbolic.transformations import hamiltonian_transformation, unitary_transformation
from qutip_symbolic.operator_utilities import drop_terms_containing, drop_c_number_terms
from qutip_symbolic.qsimplify import qsimplify
```

## The Jaynes-Cummings model

+++

The [Jaynes-Cummings model](http://en.wikipedia.org/wiki/Jaynes%E2%80%93Cummings_model) is one of the most elementary quantum mechanical models light-matter interaction. It describes a single two-level atom that interacts with a single harmonic-oscillator mode of a electromagnetic cavity.

The Hamiltonian for a two-level system in its eigenbasis (see [Two-level systems](lecture-sympy-quantum-two-level-system.ipynb)) can be written as

$$
H = \frac{1}{2}\Omega \sigma_z
$$

and the Hamiltonian of a quantum harmonic oscillator (see [Resonators and cavities](lecture-sympy-quantum-resonators.ipynb)) is

$$
H = \hbar\omega_r (a^\dagger a + 1/2)
$$

The atom interacts with the electromagnetic field produced by the cavity mode $a + a^\dagger$ through its dipole moment. The dipole-transition operators is $\sigma_x$ (which cause a transition from the two dipole states of the atom). The combined atom-cavity Hamiltonian can therefore be written in the form

$$
H = 
\hbar\omega_r (a^\dagger a + 1/2)
+ \frac{1}{2}\hbar\Omega\sigma_z 
+
\hbar
g\sigma_x(a + a^\dagger)
$$

+++

To obtain the Jaynes-Cumming Hamiltonian 

$$
H = 
\hbar\omega_r (a^\dagger a + 1/2)
%-\frac{1}{2}\Delta\sigma_x 
+ \frac{1}{2}\hbar\Omega\sigma_z 
+
\hbar
g(\sigma_+ a + \sigma_- a^\dagger)
$$

we also need to perform a rotating-wave approximation which simplifies the interaction part of the Hamiltonian. In the following we will begin with looking at how these two Hamiltonians are related.

To represent the atom-cavity Hamiltonian in SymPy we creates an instances of the operator classes `BosonOp` and `SigmaX`, `SigmaY`, and `SigmaZ`, and use these to construct the Hamiltonian (we work in units where $\hbar = 1$).

```{code-cell} ipython3
omega_r, Omega, g, Delta, t, x, Hsym = symbols("omega_r, Omega, g, Delta, t, x, H")
```

```{code-cell} ipython3
sx, sy, sz, sm, sp = SigmaX(), SigmaY(), SigmaZ(), SigmaMinus(), SigmaPlus()
a = BosonOp("a")
```

```{code-cell} ipython3
H = omega_r * Dagger(a) * a + Omega/2 * sz + g * sx * (a + Dagger(a))

Eq(Hsym, H)
```

To simplify the interaction term we carry out two unitary transformations that corresponds to moving to the interaction picture:

```{code-cell} ipython3
U = exp(I * omega_r * t * Dagger(a) * a)

U
```

```{code-cell} ipython3
H2 = hamiltonian_transformation(U, H.expand())

H2
```

```{code-cell} ipython3
U = exp(I * Omega * t * sp * sm)

U
```

```{code-cell} ipython3
H3 = hamiltonian_transformation(U, H2.expand())

H3 = H3.subs(sx, sm + sp).expand()

H3 = powsimp(H3)

H3
```

We introduce the detuning parameter $\Delta = \Omega - \omega_r$ and substitute into this expression

```{code-cell} ipython3
# trick to simplify exponents
def simplify_exp(e):
    if isinstance(e, exp):
        return exp(simplify(e.exp.expand()))

    if isinstance(e, (Add, Mul)):
        return type(e)(*(simplify_exp(arg) for arg in e.args)) 

    return e
```

```{code-cell} ipython3
H4 = simplify_exp(H3).subs(-omega_r + Omega, Delta)

H4
```

Now, in the rotating-wave approximation we can drop the fast oscillating terms containing the factors $e^{\pm i(\Omega + \omega_r)t}$

```{code-cell} ipython3
H5 = drop_terms_containing(H4, [exp( I * (Omega + omega_r) * t),
                                exp(-I * (Omega + omega_r) * t)])

H5 = drop_c_number_terms(H5.expand())

Eq(Hsym, H5)
```

This is the interaction term of in the Jaynes-Cumming model in the interaction picture. If we transform back to the Schrödinger picture we have:

```{code-cell} ipython3
U = exp(-I * omega_r * t * Dagger(a) * a)
H6 = hamiltonian_transformation(U, H5.expand())
```

```{code-cell} ipython3
U = exp(-I * Omega * t * sp * sm)
H7 = hamiltonian_transformation(U, H6.expand())
```

```{code-cell} ipython3
H8 = simplify_exp(H7).subs(Delta, Omega - omega_r)

H8 = simplify_exp(powsimp(H8)).expand()

H8 = drop_c_number_terms(H8)

H = collect(H8, g)

Eq(Hsym, H)
```

This is the Jaynes-Cumming model give above, and we have now seen that it is obtained to the dipole interaction Hamiltonian through the rotating wave approximation.

+++

## Dispersive regime

+++

In the dispersive regime, where the two-level system is detuned from the cavity by much more than the interaction strength, $\Delta \gg g$, an effective Hamiltonian can be dervied which describes the Stark shift of the two-level system (which depends on the number of photons in the cavity) and the frequency shift of the cavity (which depend on the state of the two-level system).

This effective Hamiltonian, which is correct up to second order in the small paramter $g/\Delta$, is obtained by performing the unitary transformation

$$
U = e^{\frac{g}{\Delta}(a \sigma_- - a^\dagger \sigma_+)}
$$

```{code-cell} ipython3
U = exp((x * (a * sp - Dagger(a) * sm)).expand())

U
```

```{code-cell} ipython3
#H1 = unitary_transformation(U, H, allinone=True, expansion_search=False, N=3).expand()
#H1 = qsimplify(H1)
#H1
```

```{code-cell} ipython3
H1 = hamiltonian_transformation(U, H, expansion_search=False, N=3).expand()

H1 = qsimplify(H1)

H1
```

```{code-cell} ipython3
H2 = drop_terms_containing(H1.expand(), [x**3, x**4])

H2
```

```{code-cell} ipython3
H3 = H2.subs(x, g/Delta)

H3
```

```{code-cell} ipython3
H4 = drop_c_number_terms(H3)

H4
```

```{code-cell} ipython3
H5 = collect(H4, [Dagger(a) * a, sz])

H5
```

Now move to a frame co-rotating with the qubit and oscillator frequencies:

```{code-cell} ipython3
H5.expand()
```

```{code-cell} ipython3
U = exp(I * omega_r * t * Dagger(a) * a)
```

```{code-cell} ipython3
H6 = hamiltonian_transformation(U, H5.expand()); H6
```

```{code-cell} ipython3
U = exp(I * Omega * t * Dagger(sm) * sm)
```

```{code-cell} ipython3
H7 = hamiltonian_transformation(U, H6.expand()); H7
```

Now, since we are in the dispersive regime $|\Omega-\omega_r| \gg g$, we can do a rotating-wave approximation and drop all the fast rotating terms in the Hamiltonian above:

```{code-cell} ipython3
H8 = drop_terms_containing(H7, [exp(I * omega_r * t), exp(-I * omega_r * t),
                                exp(I * Omega * t), exp(-I * Omega * t)])

H8
```

```{code-cell} ipython3
H9 = qsimplify(H8)

H9 = collect(H9, [Dagger(a) * a, sz])

H9
```

Now move back to the lab frame:

```{code-cell} ipython3
U = exp(-I * omega_r * t * Dagger(a) * a)
```

```{code-cell} ipython3
H10 = hamiltonian_transformation(U, H9.expand()); H10
```

```{code-cell} ipython3
U = exp(-I * Omega * t * Dagger(sm) * sm)
```

```{code-cell} ipython3
H11 = hamiltonian_transformation(U, H10.expand()); H11
```

```{code-cell} ipython3
H12 = qsimplify(H11)

H12 = collect(H12, [Dagger(a) * a, sz])

H12 = H12.subs(omega_r, Omega-Delta).expand().collect([Dagger(a)*a, sz]).subs(Omega-Delta,omega_r)

Eq(Hsym, H12)
```

This is the Hamiltonian of the Jaynes-Cummings model in the the dispersive regime. It can be interpreted as the resonator having a qubit-state-dependent frequency shift, or alternatively that the qubit is feeling a resonator-photon-number dependent Stark-shift.

+++

## Versions

```{code-cell} ipython3
import sympy
import qutip_symbolic
print("sympy:", sympy.__version__)
print("qutip_symbolic", qutip_symbolic.__version__)
```

```{code-cell} ipython3

```
