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

# Lecture 2 - Symbolic quantum mechanics using SymPsi - Two-level systems

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
from sympy import collect, expand, simplify, symbols, Eq, I, exp, atan, tan, sqrt, init_printing
init_printing()
```

```{code-cell} ipython3
from qutip_symbolic.compat.pauli import SigmaX, SigmaY, SigmaZ
```

```{code-cell} ipython3
from qutip_symbolic.transformations import hamiltonian_transformation
from qutip_symbolic.operator_utilities import split_coeff_operator
```

## Introduction: Hamiltonian and Pauli matrices

+++

One of the simplest possible quantum system is the two-level system, such as a spin-$1/2$, or an atom or artificial-atom (nano-device) with (effectively) only two quantum states. 

The Hamiltonian for a two-level quantum system is conveniently written in terms of the [Pauli matrices](http://en.wikipedia.org/wiki/Pauli_matrices) $\sigma_x$, $\sigma_y$ and $\sigma_z$, and their annihilation and creation operators $\sigma_-$ and $\sigma_+$. In terms of these operators we can write the Hamiltonian on the form:

$$
H = -\frac{1}{2}\Delta\sigma_x - \frac{1}{2}\epsilon\sigma_z 
$$

where $\epsilon$ is the bare energy splitting and $\Delta$ is the tunneling rate between the two states of the two-level system.

In SymPy we can represent this Hamiltonian as:

```{code-cell} ipython3
theta, t = symbols("theta, t")
eps, Delta, Omega = symbols("epsilon, Delta, Omega", positive=True)
Hsym = symbols("H")
```

```{code-cell} ipython3
sx, sy, sz = SigmaX(), SigmaY(), SigmaZ()
```

```{code-cell} ipython3
H = -eps/2 * sz - Delta/2 * sx

Eq(Hsym, H)
```

## Instantaneous eigenbasis

+++

Is often convenient to perform basis transformation that simplifies the Hamiltonian. For example, we can transform the Hamiltonian to the eigenbasis (where the Hamiltonian is diagonal, that is only containing a $\sigma_z$ term) by applying the unitary tranformation:

```{code-cell} ipython3
U = exp(I * theta/2 * sy); U
```

This unitary tranformation transforms the operators in the Hamiltonian according to these well-known relations:

```{code-cell} ipython3
hamiltonian_transformation(U, sx)
```

```{code-cell} ipython3
hamiltonian_transformation(U, sz)
```

so the Hamiltonian after this transformation takes the form

```{code-cell} ipython3
H1 = hamiltonian_transformation(U, H)

H1
```

```{code-cell} ipython3
H2 = collect(H1.expand(), (sx, sz))

H2
```

In the eigenbasis we require the coefficient of $\sigma_x$ to be zero, so we have the condition:

```{code-cell} ipython3
c, o = split_coeff_operator(H2.args[0])

Eq(c, 0)
```

with the solution

```{code-cell} ipython3
Eq(tan(theta), Delta/eps)
```

Substituting this into the Hamiltonian results in

```{code-cell} ipython3
H3 = simplify(H2.subs(theta, atan(Delta/eps)))

H3
```

Now introduce $\Omega = \sqrt{\Delta^2 + \epsilon^2}$, which is the eigenenergies of the two-level system:

```{code-cell} ipython3
H3.subs(Delta, sqrt(Omega ** 2 - eps ** 2))
```

In summary, to reach this basis, we have transformed $\sigma_x$ and $\sigma_z$ as:

```{code-cell} ipython3
hamiltonian_transformation(U, sx)
```

```{code-cell} ipython3
hamiltonian_transformation(U, sz)
```

and chosen $\theta = \arctan(\Delta/\epsilon)$.

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
