---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.11.2
---

# Lectures - Symbolic quantum mechanics with SymPsi - Table of content

+++

Author: J. R. Johansson (robert@riken.jp), [http://jrjohansson.github.io](http://jrjohansson.github.io).

Status: Preliminary (work in progress)

This notebook is part of a series of IPython notebooks on symbolic quantum mechanics computations using 
[SymPy](http://sympy.org) and [SymPsi](http://www.github.com/jrjohansson/sympsi). SymPsi is an experimental fork and extension of the [`sympy.physics.quantum`](http://docs.sympy.org/dev/modules/physics/quantum/) module in SymPy. The latest version of this notebook is available at [http://github.com/jrjohansson/sympy-quantum-notebooks](http://github.com/jrjohansson/sympy-quantum-notebooks), and the other notebooks in this lecture series are also indexed at [http://jrjohansson.github.io](http://jrjohansson.github.com).

Requirements: A recent version of SymPy and the latest development version of SymPsi is required to execute this notebook. Instructions for how to install SymPsi is available [here](http://www.github.com/jrjohansson/sympsi).

Disclaimer: The SymPsi module is still under active development and may change in behavior without notice, and the intention is to move some of its features to [`sympy.physics.quantum`](http://docs.sympy.org/dev/modules/physics/quantum/) when they matured and have been tested. However, these notebooks will be kept up-to-date the latest versions of SymPy and SymPsi.

+++

## About these notebooks

+++

These notebooks are vaguely written in the form of lecture notes. However, the point of these notebooks is not to give an introduction to quantum mechanics or even the topics covered by the notebooks, but rather demonstrate how analytical calculations that are often necessary when working with these quantum mechanical systems can be carried out symbolically. 

+++

## Table of content

+++

 * [Basics](lecture-sympsi-basics.ipynb)
   * About SymPy and the SymPsi module
   * Getting started
 * [Two-level systems](lecture-sympsi-two-level-system.ipynb)
   * Eigenbasis
 * [Resonators and cavities](lecture-sympsi-resonator.ipynb)
   * Classical driving fields
   * Rotating frames and the rotating-wave approximation
   * Optical parametric oscillator
 * [Atom and cavity](lecture-sympsi-atom-cavity.ipynb)
   * The Jaynes-Cumming model
   * Dispersive regime
 * [Optomechanics](lecture-sympsi-optomechanics.ipynb)
   * Linear regime
     * Red sideband
     * Blue sideband
   * Kerr effect and photon blockade
 * Input/Output theory
   * Quantum Langevin equation
   * One-sided cavity
 * Master equations
   * Adjoint master equation
 * [Semiclassical equations of motion](lecture-sympsi-semiclassical-eqm.ipynb)
 * Composite systems
