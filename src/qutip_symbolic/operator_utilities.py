from sympy import Add, Mul, Pow, Symbol, Expr
#from .compat.operator import Operator
from sympy.physics.quantum import Operator
from sympy.physics.quantum.operatorordering import normal_ordered_form


def split_coeff_operator(e):
    """
    Split a product of coefficients, commuting variables and quantum
    operators into two factors containing the commuting factors and the
    quantum operators, resepectively.

    Returns:
    c_factor, o_factors:
        Commuting factors and noncommuting (operator) factors
    """
    if isinstance(e, Symbol):
        if e.is_commutative:
            return e, 1
        return 1, e

    if isinstance(e, Pow):
        c, o = split_coeff_operator(e.base)
        return (c ** e.exp, o ** e.exp)

    if isinstance(e, Mul):
        c_args = []
        o_args = []

        for arg in e.args:
            if isinstance(arg, Add):
                if arg.is_commutative:
                    c_args.append(arg)
                else:
                    o_args.append(arg)
            else:
                c, o = split_coeff_operator(arg)
                if c != 1:
                    c_args.append(c)
                if o != 1:
                    o_args.append(o)

        return Mul(*c_args), Mul(*o_args)

    if isinstance(e, Add):
        return [split_coeff_operator(arg) for arg in e.args]

    if isinstance(e, Expr):  # also covers the Operator case
        # is_commuative may also be None
        if e.is_commutative is True:
            return e, 1
        elif e.is_commutative is False:
            return 1, e

    raise TypeError(f"split_coeff_operator: {e!r} has unsupported type")


def extract_operators(e, independent=False):
    """
    Return a list of unique quantum operator products in the
    expression e.
    """
    ops = []

    if isinstance(e, Operator):
        ops.append(e)

    elif isinstance(e, Add):
        for arg in e.args:
            ops += extract_operators(arg, independent=independent)

    elif isinstance(e, Mul):
        for arg in e.args:
            ops += extract_operators(arg, independent=independent)
    else:
        if debug:
            print("Unrecongized type: %s: %s" % (type(e), str(e)))

    return list(set(ops))


def extract_operator_products(e, independent=False):
    """
    Return a list of unique normal-ordered quantum operator products in the
    expression e.
    """
    ops = []

    if isinstance(e, Operator):
        ops.append(e)

    elif isinstance(e, Add):
        for arg in e.args:
            ops += extract_operator_products(arg, independent=independent)

    elif isinstance(e, Mul):
        c, o = split_coeff_operator(e)
        if o != 1:
            ops.append(o)
    else:
        if debug:
            print("Unrecongized type: %s: %s" % (type(e), str(e)))

    no_ops = []
    for op in ops:
        no_op = normal_ordered_form(op.expand(), independent=independent)
        if isinstance(no_op, (Mul, Operator, Pow)):
            no_ops.append(no_op)
        elif isinstance(no_op, Add):
            for sub_no_op in extract_operator_products(no_op, independent=independent):
                no_ops.append(sub_no_op)
        else:
            raise ValueError("Unsupported type in loop over ops: %s: %s" %
                             (type(no_op), no_op))

    return list(set(no_ops))


def drop_terms_containing(e, e_drops):
    """
    Drop terms contaning factors in the list e_drops
    """
    if isinstance(e, Add):
        # fix this
        # e = Add(*(arg for arg in e.args if not any([e_drop in arg.args
        #                                            for e_drop in e_drops])))

        new_args = []

        for term in e.args:

            keep = True
            for e_drop in e_drops:
                if e_drop in term.args:
                    keep = False

                if isinstance(e_drop, Mul):
                    if all([(f in term.args) for f in e_drop.args]):
                        keep = False

            if keep:
                # new_args.append(arg)
                new_args.append(term)
        e = Add(*new_args)
        # e = Add(*(arg.subs({key: 0 for key in e_drops}) for arg in e.args))

    return e


def drop_c_number_terms(e):
    """
    Drop commuting terms from the expression e
    """
    if isinstance(e, Add):
        return Add(*(arg for arg in e.args if not arg.is_commutative))

    return e


def subs_single(O, subs_map):

    if isinstance(O, Operator):
        if O in subs_map:
            return subs_map[O]
        else:
            print("warning: unresolved operator: ", O)
            return O
    elif isinstance(O, Add):
        new_args = []
        for arg in O.args:
            new_args.append(subs_single(arg, subs_map))
        return Add(*new_args)

    elif isinstance(O, Mul):
        new_args = []
        for arg in O.args:
            new_args.append(subs_single(arg, subs_map))
        return Mul(*new_args)

    elif isinstance(O, Pow):
        return Pow(subs_single(O.base, subs_map), O.exp)

    else:
        return O
