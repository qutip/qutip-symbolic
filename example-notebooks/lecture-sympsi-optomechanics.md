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

# Lecture 5 - Symbolic quantum mechanics using SymPsi - Optomechanics

+++

<style>
p {
    font-family: "Liberation Serif", serif;
    font-size: 12pt;
}
</style>

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
from sympy import collect, symbols, exp, conjugate, Eq, I, Symbol, init_printing
from sympy.physics.quantum import Dagger
from sympy.physics.quantum.boson import BosonOp
from sympy.physics.quantum.operatorordering import normal_ordered_form
init_printing()
```

```{code-cell} ipython3
from qutip_symbolic.transformations import hamiltonian_transformation
from qutip_symbolic.operator_utilities import drop_terms_containing, drop_c_number_terms
```

## Optomechanical system

+++

Consider the standard Hamiltonian for an optomechanical system, including a classical driving signal on the optical mode:

$$
H = \hbar\omega_a a^\dagger a + \hbar \omega_b b^\dagger b - \hbar g a^\dagger a (b + b^\dagger) + (A e^{-i\omega_d t} + A^* e^{i\omega_d t})(a + a^\dagger)
$$

```{code-cell} ipython3
omega_a, omega_b, g, A, Delta, t = symbols("omega_a, omega_b, g, A, Delta, t", positive=True)
Hsym, omega_d = symbols("H, omega_d")
```

```{code-cell} ipython3
a, b = BosonOp("a"), BosonOp("b")
```

```{code-cell} ipython3
H0 = omega_a * Dagger(a) * a + omega_b * Dagger(b) * b - g * Dagger(a) * a * (b + Dagger(b))

Hdrive = (A * exp(-I * omega_d * t) + conjugate(A) * exp(I * omega_d * t)) * (a + Dagger(a))

H = H0 + Hdrive

Eq(Hsym, H)
```

### Linearized interaction

+++

First we apply the unitary transformation $U = e^{i \omega_d a^\dagger a t}$:

```{code-cell} ipython3
U = exp(I * Dagger(a) * a * omega_d * t)

U
```

```{code-cell} ipython3
H1 = hamiltonian_transformation(U, H, independent=True)

H1
```

We can now perform a rotating-wave approximation (RWA) by eliminating all terms that rotate with frequencies $2\omega_d$:

```{code-cell} ipython3
H2 = drop_terms_containing(H1.expand(), [exp(-2*I*omega_d*t), exp(2*I*omega_d*t)])

Eq(Symbol("H_{rwa}"), H2)
```

Introduce the detuning $\Delta = \omega_a - \omega_d$:

```{code-cell} ipython3
H3 = H2.subs(omega_a, Delta + omega_d).expand()

H3
```

To eliminate the coherent part of the state of the cavity mode we apply the unitary displacement operator $U = e^{\alpha a^\dagger - \alpha^*a}$:

```{code-cell} ipython3
alpha = symbols("alpha")
```

```{code-cell} ipython3
UH = Dagger(a) * alpha - conjugate(alpha) * a
U = exp(UH)

U
```

```{code-cell} ipython3
H4 = hamiltonian_transformation(U, H3, independent=True)

H4
```

Now want to cancel out the drivng terms so we set $A - \Delta \alpha = 0$, i.e. $\alpha = A/\Delta$:

```{code-cell} ipython3
H5 = H4.expand().subs({A: alpha * Delta, conjugate(alpha): alpha}) 

H5 = collect(H5, [g * Dagger(a) * a, - alpha * g])

H5
```

Drop C-numbers from the Hamiltonian:

```{code-cell} ipython3
H6 = drop_c_number_terms(H5)

H6
```

Now, if driving strength is large, so that $\alpha \gg 1$, we can drop the nonlinear interaction term, and we have an linear effective coupling:

```{code-cell} ipython3
H7 = H6.subs(g * Dagger(a) * a, 0)

e = (a + Dagger(a)) * (b + Dagger(b))
H7 = H7.subs(e.expand(), e)

Eq(Hsym, H7)
```

This linearlized optomechanical Hamiltonian has at least two interesting regimes:

+++

### Red sideband

+++

The red sideband regime occurs when the detuning is $\Delta = \omega_b$. In this case, if we move to a frame rotating with the driving field, we obtain:

```{code-cell} ipython3
H1 = H7
H1
```

```{code-cell} ipython3
U = exp(I * Dagger(a) * a * Delta * t)
U
```

```{code-cell} ipython3
H2 = hamiltonian_transformation(U, H1, independent=True)
H2
```

```{code-cell} ipython3
U = exp(I * Dagger(b) * b * omega_b * t)
U
```

```{code-cell} ipython3
H3 = hamiltonian_transformation(U, H2, independent=True)
H3
```

If we substitute $\Delta = \omega_b$ we get:

```{code-cell} ipython3
H4 = H3.expand().subs(Delta, omega_b)
H4
```

Now we can do a rotating-wave approximation and get rid of terms rotating at angular frequencies $2\omega_b$, and then transform back to the original frame:

```{code-cell} ipython3
H5 = drop_terms_containing(H4, [exp(+2 * I * omega_b * t), exp(-2 * I * omega_b * t)])
H5
```

```{code-cell} ipython3
U = exp(-I * Dagger(a) * a * omega_b * t)  # Delta = omega_b
H6 = hamiltonian_transformation(U, H5, independent=True)
U = exp(-I * Dagger(b) * b * omega_b * t)
H7 = hamiltonian_transformation(U, H6, independent=True)
H7 = collect(H7, [alpha**2, g])
H7
```

Now the interaction term is $a^\dagger b + ab^\dagger$, which is a swapping interaction that can be used for state transfer or energy transfer in for example cooling application (side-band cooling).

+++

### Blue sideband

+++

If, instead, we choose a driving frequency such that $\Delta = -\omega_b$, we obtain:

```{code-cell} ipython3
H4 = H3.expand().subs(Delta, -omega_b)
H4
```

As before, we do a rotating-wave approximation to get rid of fast oscillating terms:

```{code-cell} ipython3
H5 = drop_terms_containing(H4, [exp(+2 * I * omega_b * t), exp(-2 * I * omega_b * t)])
H5
```

and moving back to the original frame results in:

```{code-cell} ipython3
U = exp( I * Dagger(a) * a * omega_b * t)  # Delta = -omega_b
H6 = hamiltonian_transformation(U, H5, independent=True)
U = exp(-I * Dagger(b) * b * omega_b * t)
H7 = hamiltonian_transformation(U, H6, independent=True)
H7 = collect(H7, [alpha**2, g])
H7
```

Here, instead of a swap interaction we have obtained an interaction on the form $a^\dagger b^\dagger + a b$, which is the parametric amplification Hamiltonian. It can be used to parametrically amplify the states of the optical and mechanical modes, and generated interesting nonclassical states like Schrodinger-cat states.

+++

## Nonlinear regime: effective Kerr nonlinearity

+++

In the regime where $g \sim \omega_b$, the effect of the coupling to the mechanical mode $b$ on the optical mode $a$ is an effective Kerr-nonlinearity, i.e., a term on the form $(a^\dagger a)^2$ in the Hamiltonian. To see this we can perform the so-called polariton defined by the unitary

$$
U = \exp\left(-\frac{g}{\omega_b} a^\dagger a (b^\dagger - b)\right)
$$

```{code-cell} ipython3
H0
```

```{code-cell} ipython3
x = symbols("x")
```

```{code-cell} ipython3
U = exp(- x * Dagger(a) * a * (Dagger(b) - b))

U
```

```{code-cell} ipython3
H1 = H0.subs(A, 0)

H1
```

```{code-cell} ipython3
H2 = hamiltonian_transformation(U, H1, independent=True, expansion_search=False, N=2)

H3 = normal_ordered_form(H2.expand(), independent=True)

H3
```

```{code-cell} ipython3
H4 = H3.subs({x**2: 0, x**3: 0, x: g/omega_b})  # neglect higher order terms

H4 = collect(H4, Dagger(a)*a)

H4
```

In this Hamiltonian, the mechanical and optical mode is effectively decoupled, but the influence of the mechanical mode on the optical mode is described by the induced Kerr nonlinearity for the optical mode.

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
