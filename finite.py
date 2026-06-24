from __init__ import *
from sympy.matrices import *
from sympy import pprint


# The goal of this library is to compute a generic element of a parabolic subgroup given a nilpotent f.
# This also only works in type A lol
def bracket(x, y):
    return x * y - y * x


def generate_h(size):
    return diag(*[
        x for x in range(size - 1, -size, -2)
    ])


def sl2_triple(e):
    e_normal = e.jordan_form()[1]
    _, e_cells = e.jordan_cells()

    h_cells = []

    for cell in e_cells:
        size = cell.shape
        assert size[0] == size[1], "Non-square Jordan block???"

        h_cells.append(generate_h(size[0]))

    h = diag(*h_cells)
    assert bracket(h, e_normal) == 2 * e_normal, "sl2 triple h failed consistency"

    f = None  # todo

    return e_normal, f, h


def parabolic(e: Matrix):
    e, f, h = sl2_triple(e)

    # h-weights >= 0

    return h
