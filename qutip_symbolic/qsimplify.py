import warnings
from sympy import Add, exp, Mul, Pow, simplify
from sympy.physics.quantum.operatorordering import normal_ordered_form


def qsimplify(e_orig, _n=0):
    """
    Simplify an expression containing operators.
    """
    if _n > 15:
        warnings.warn("Too high level or recursion, aborting")
        return e_orig

    e = normal_ordered_form(e_orig)

    if isinstance(e, Add):
        return Add(*(qsimplify(arg, _n=_n+1) for arg in e.args))

    elif isinstance(e, Pow):
        return Pow(*(qsimplify(arg, _n=_n+1) for arg in e.args))

    elif isinstance(e, exp):
        return exp(*(qsimplify(arg, _n=_n+1) for arg in e.args))

    elif isinstance(e, Mul):
        args1 = tuple(arg for arg in e.args if arg.is_commutative)
        args2 = tuple(arg for arg in e.args if not arg.is_commutative)
        #x = 1
        #for y in args2:
        #    x = x * y

        x = 1
        for y in reversed(args2):
            x = y * x

        if isinstance(x, Mul):
            args2 = x.args
            x = 1
            for y in args2:
                x = x * y

        e_new = simplify(Mul(*args1)) * x

        if e_new == e:
            return e
        else:
            return qsimplify(e_new.expand(), _n=_n+1)

    if e == e_orig:
        return e
    else:
        return qsimplify(e, _n=_n+1).expand()
