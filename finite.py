import itertools

from __init__ import *
from sympy.matrices import *
from sympy import *


# The goal of this library is to compute a generic element of a parabolic subgroup given a nilpotent f.
# This also only works in type A lol
# More specifically it only works for GL_n


# [I_p \lambda^\vee(t) G(O) \cap U_P(K) G(O)]/G(O) \neq \emptyset?


# Proposal:
# For g \in GL_n(K), figure out how to quickly compute its Smith normal = Cartan decomposition
# Then refine using Bruhat decomposition
# That says if you have a matrix, you get its I,G(O) double coset
# Then can read I_P,G(O) off of that
# Then sample "random"ish matrices in U_P(K)
def bracket(x, y):
    return x * y - y * x


def generate_h(size):
    return diag(*[
        x for x in range(size - 1, -size, -2)
    ])


def sl2_triple(e: Matrix):
    e_normal: Matrix
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


def parabolic(e: Matrix, name='n'):
    e, f, h = sl2_triple(e)

    # h weights >= 0

    assert e.shape[0] == e.shape[1]
    size = e.shape[0]

    ad_h = []

    for i, j in itertools.product(range(size), range(size)):  # Choice: [h, -] (standard ad)
        E_ij = zeros(size, size)
        E_ij[i, j] = 1
        result = bracket(h, E_ij)

        ad_h.append(result.reshape(size * size, 1).transpose())

    ad_h = Matrix(ad_h).transpose()  # needs to be tested lol

    eigenvectors = ad_h.eigenvects()
    parabolic_basis = dict()
    for λ, count, vecs in eigenvectors:
        if λ >= 0:
            parabolic_basis[λ] = [vec.reshape(size, size) for vec in vecs]

    return _formal_span(parabolic_basis, name=name)


def _formal_span(eigenbasis, name='N'):
    total = []
    size = 0
    for λ, basis in eigenbasis.items():
        elts = symbols(f'{name}_{λ}(1:{len(basis) + 1})')
        total.extend([
            sym * matrix for sym, matrix in zip(elts, basis)
        ])
        size = basis[0].shape[0]

    return sum(total, start=zeros(size, size))


def Parabolic(e: Matrix, name='N'):
    return parabolic(e, name)  # fuck baker campbell hausdorff
    # P is positive eigenspace of h acting on G by conjugation
