from __init__ import *
from symbols import *
from sympy.matrices import *


def Cocharacter(size=dimension):  # Generic cocharacter
    r"""
    Generic cocharacter of the form \lambda^\vee(t) \in G(K).
    :param size: The size of the cocharacter
    :return: A matrix for a cocharacter
    """
    λs = create_symbols('λ', size)

    return diag(*[t ** λs[i] for i in range(size)])
