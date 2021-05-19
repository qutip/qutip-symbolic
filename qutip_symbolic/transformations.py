from sympy import I, diff, exp
from .bch import bch_expansion
from .operator_utilities import extract_operators, subs_single

debug = False  # TODO: Replace with logging


def unitary_transformation(U, O, N=6, collect_operators=None,
                           independent=False, allinone=False,
                           expansion_search=True):
    """
    Perform a unitary transformation

        O = U O U^\\dagger

    and automatically try to identify series expansions in the resulting
    operator expression.
    """
    if not isinstance(U, exp):
        raise ValueError("U must be a unitary operator on the form "
                         "U = exp(A)")

    A = U.exp

    if debug:
        print("unitary_transformation: using A = ", A)

    if allinone:
        return bch_expansion(A, O, N=N, collect_operators=collect_operators,
                             independent=independent,
                             expansion_search=expansion_search)
    else:
        ops = extract_operators(O.expand())
        ops_subs = {op: bch_expansion(A, op, N=N,
                                      collect_operators=collect_operators,
                                      independent=independent,
                                      expansion_search=expansion_search)
                    for op in ops}

        #  return O.subs(ops_subs, simultaneous=True) # XXX: this this
        return subs_single(O, ops_subs)


def hamiltonian_transformation(U, H, N=6, collect_operators=None,
                               independent=False, expansion_search=True):
    """
    Apply an unitary basis transformation to the Hamiltonian H:

        H = U H U^\\dagger -i U d/dt(U^\\dagger)

    """
    t = [s for s in U.exp.free_symbols if str(s) == 't']
    if t:
        t = t[0]
        H_td = - I * U * diff(exp(-U.exp), t)
    else:
        H_td = 0

    # H_td = I * diff(U, t) * exp(- U.exp)  # hack: Dagger(U) = exp(-U.exp)
    H_st = unitary_transformation(U, H, N=N,
                                  collect_operators=collect_operators,
                                  independent=independent,
                                  expansion_search=expansion_search)
    return H_st + H_td
