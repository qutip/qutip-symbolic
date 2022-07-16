from sympy.core.symbol import symbols
from sympy.physics.quantum.commutator import Commutator as Comm

from qutip_symbolic.commutators import recursive_commutator

A, B = symbols('A, B', commutative=False)


def test_recursive_commutators_n_eq_1():
    C = recursive_commutator(A, B, n=1)
    assert C == Comm(A, B)


def test_recursive_commutators_n_eq_2():
    C = recursive_commutator(A, B, n=2)
    assert C == Comm(A, Comm(A, B))
