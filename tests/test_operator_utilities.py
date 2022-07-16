import pytest
from sympy.core.symbol import symbols, Expr
from sympy.physics.quantum import Operator

from qutip_symbolic.operator_utilities import split_coeff_operator

a, b = symbols('a, b')
A, B = symbols('A, B', commutative=False)
Op = Operator()
E = Expr()


def test_split_coeff_operator():
    assert split_coeff_operator(A) == (1, A)
    assert split_coeff_operator(Op) == (1, Op)
    assert split_coeff_operator(a) == (a, 1)
    assert split_coeff_operator(A * B) == (1, A * B)
    assert split_coeff_operator(a * A) == (a, A)
    assert split_coeff_operator(a * A * b * B) == (a * b, A * B)
    assert split_coeff_operator(a * A + b * B) == [(a, A), (b, B)]
    assert split_coeff_operator(A + B) == [(1, A), (1, B)]
    assert split_coeff_operator(a + b) == [(a, 1), (b, 1)]
    assert split_coeff_operator(a * (A + B)) == (a, A + B)
    with pytest.raises(TypeError) as err:
        split_coeff_operator(E)
    assert str(err.value) == (
        "split_coeff_operator: Expr() has unsupported type"
    )
