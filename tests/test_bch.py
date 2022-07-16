from sympy.core.symbol import symbols
from sympy.physics.quantum.commutator import Commutator as Comm

from qutip_symbolic.bch import _bch_expansion

A, B = symbols('A, B', commutative=False)


def test_simple_expand():
    C = _bch_expansion(A, B, N=4)
    assert C == (
        B
        + Comm(A, B)
        + Comm(A, Comm(A, B)) / 2
        + Comm(A, Comm(A, Comm(A, B))) / 6
    )
