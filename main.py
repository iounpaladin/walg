from sympy import *
from sympy.interactive.printing import init_printing
from affine import *
from finite import *
from symbols import *

init_printing(use_latex=False)


# Code goal: given cocharacter lambda, determine if i lambda(t) g in N(k) has a solution with i in I_p, g in G(O)
# with given bounds on the powers of t appearing in i and g

def iwahori(size=dimension):
    pass


def g_o(size=dimension):
    pass


e = Matrix([[0, 1, 0, 0], [0, 0, 0, 0], [0, 0, 0, 1], [0, 0, 0, 0]])

pprint(parabolic(e))
