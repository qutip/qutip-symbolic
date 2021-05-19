from sympy import (
    Add, collect, Dummy, cos, cosh, sin, sinh, exp, factorial, I, Mul, Pow,
    simplify
)
from sympy.physics.quantum.operatorordering import normal_ordered_form
from .commutators import recursive_commutator
from .operator_utilities import extract_operator_products, split_coeff_operator
from .qsimplify import qsimplify

debug = False  # TODO: replace with logging


def _order(e):
    fs = list(e.free_symbols)
    if isinstance(e, Pow) and e.base == fs[0]:
        return e.exp
    elif isinstance(e, Mul):
        o = sum([_order(arg) for arg in e.args])
        return o
    elif isinstance(e, Add):
        o = max([_order(arg) for arg in e.args])
        return o
    elif e.is_Symbol:
        return 1
    else:
        return 0


def _lowest_order_term(e):

    if isinstance(e, Add):
        min_order = _order(e.args[0])
        min_expr = e.args[0]
        for arg in e.args:
            arg_order = _order(arg)
            if arg_order < min_order:
                min_order = arg_order
                min_expr = arg
        return min_expr, min_order
    else:
        return e, _order(e)


def _expansion_search(e, rep_list, N):
    """
    Search for and substitute terms that match a series expansion of
    fundamental math functions.

    e: expression

    rep_list: list containing dummy variables

    """
    if e.find(I):
        is_complex = True
    else:
        is_complex = False

    if debug:
        print("_expansion_search: ", e)

    try:
        dummy = Dummy()

        flist0 = [exp, lambda x: exp(-x), cos, cosh]

        flist1 = [lambda x: (exp(x) - 1) / x,
                  lambda x: (1 - exp(-x)) / x,
                  lambda x: sin(x) / x,
                  lambda x: sinh(x) / x]

        flist2 = [lambda x: (1 - cos(x))/(x**2/2),
                  lambda x: (cosh(x)-1)/(x**2/2)]

        if is_complex:
            iflist0 = [lambda x: exp(I*x),
                       lambda x: exp(-I*x)]
            iflist1 = [lambda x: (exp(I*x) - 1) / (I*x),
                       lambda x: (1 - exp(-I*x)) / (I*x)]

            flist0 = iflist0 + flist0
            flist1 = iflist1 + flist1

        flist = [flist0, flist1, flist2]
        fseries = {}

        if isinstance(e, Mul):
            e_args = [e]
        elif isinstance(e, Add):
            e_args = e.args
        else:
            return e

        newargs = []
        for e in e_args:
            if isinstance(e, Mul):
                c, nc = e.args_cnc()
                if nc and c:
                    c_expr = Mul(*c).expand()
                    d, d_order = _lowest_order_term(c_expr)
                    c_expr_normal = (c_expr / d).expand()
                    c_expr_subs = c_expr_normal

                    for alpha in rep_list:
                        if alpha not in c_expr_subs.free_symbols:
                            continue

                        for f in flist[d_order]:
                            if f not in fseries.keys():
                                fseries[f] = f(dummy).series(
                                    dummy, n=N-d_order).removeO()
                            c_expr_subs = c_expr_subs.subs(
                                fseries[f].subs(dummy, alpha), f(alpha))
                            if c_expr_subs != c_expr_normal:
                                break
                        if c_expr_subs != c_expr_normal:
                            break
                    newargs.append(d * c_expr_subs * Mul(*nc))
                else:
                    newargs.append(e)
            else:
                newargs.append(e)

        return Add(*newargs)

    except Exception as e:
        print("Failed to identify series expansions: " + str(e))
        raise


def _bch_expansion(A, B, N=10):
    """
    Baker–Campbell–Hausdorff formula:

    e^{A} B e^{-A} = B + 1/(1!)[A, B] +
                     1/(2!)[A, [A, B]] + 1/(3!)[A, [A, [A, B]]] + ...
                   = B + Sum_n^N 1/(n!)[A, B]^n

    Truncate the sum at N terms.
    """
    e = B
    for n in range(1, N):
        e += recursive_commutator(A, B, n=n) / factorial(n)

    return e


def bch_expansion(A, B, N=6, collect_operators=None, independent=False,
                  expansion_search=True):

    # Use BCH expansion of order N

    if debug:
        print("bch_expansion: ", A, B)

    cno = split_coeff_operator(A)
    if isinstance(cno, list):
        nvar = len(cno)
        c_list = []
        o_list = []
        for n in range(nvar):
            c_list.append(cno[n][0])
            o_list.append(cno[n][1])
    else:
        nvar = 1
        c_list, o_list = [cno[0]], [cno[1]]

    if debug:
        print("A coefficient: ", c_list)

    rep_list = []
    var_list = []

    for n in range(nvar):
        rep_list.append(Dummy())

        coeff, sym = c_list[n].as_coeff_Mul()
        if isinstance(sym, Mul):
            sym_ = simplify(sym)
            if I in sym_.args:
                var_list.append(sym_/I)
            elif any([isinstance(arg, exp) for arg in sym_.args]):
                nexps = Mul(*[arg for arg in sym_.args
                              if not isinstance(arg, exp)])
                exps = Mul(*[arg for arg in sym_.args if isinstance(arg, exp)])

                if I in simplify(exps).exp.args:
                    var_list.append(nexps)
                else:
                    var_list.append(sym_)
            else:
                var_list.append(sym_)
        else:
            var_list.append(sym)

    A_rep = A.subs({var_list[n]: rep_list[n] for n in range(nvar)})

    e_bch_rep = _bch_expansion(A_rep, B, N=N).doit(independent=independent)

    if debug:
        print("simplify: ")

    e = qsimplify(normal_ordered_form(e_bch_rep.expand(),
                                      recursive_limit=25,
                                      independent=independent).expand())
    if debug:
        print("extract operators: ")

    ops = extract_operator_products(e, independent=independent)

    # make sure that product operators comes first in the list
    ops = list(reversed(sorted(ops, key=lambda x: len(str(x)))))

    if debug:
        print("operators in expression: ", ops)

    if collect_operators:
        e_collected = collect(e, collect_operators)
    else:
        e_collected = collect(e, ops)

    if debug:
        print("search for series expansions: ", expansion_search)

    if debug:
        print("e_collected: ", e_collected)

    if expansion_search and c_list:
        e_collected = _expansion_search(e_collected, rep_list, N)
        e_collected = e_collected.subs({rep_list[n]: var_list[n]
                                        for n in range(nvar)})

        return e_collected
    else:
        return e_collected.subs(
            {rep_list[n]: var_list[n] for n in range(nvar)})
