
#from .compat.commutator import Commutator
from sympy.physics.quantum.commutator import Commutator


def recursive_commutator(a, b, n=1):
    """
    Generate a recursive commutator of order n:

    [a, b]_1 = [a, b]
    [a, b]_2 = [a, [a, b]]
    [a, b]_3 = [a, [a, b]_2] = [a, [a, [a, b]]]
    ...

    """
    if n == 1:
        return Commutator(a, b)
    else:
        return Commutator(a, recursive_commutator(a, b, n-1))
