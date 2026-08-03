from sympy import *
from sympy.interactive.printing import init_printing
from sympy.matrices.normalforms import *
from sympy.polys.puiseux import puiseux_ring

init_printing(use_latex=False)

# Proposal:
# For g \in GL_n(K), figure out how to quickly compute its Smith normal = Cartan decomposition
# Then refine using Bruhat decomposition
# That says if you have a matrix, you get its I,G(O) double coset
# Then can read I_P,G(O) off of that
# Then sample "random"ish matrices in U_P(K)

# Important note: we are taking bounded terms of t so we take K = C(t) not C((t))
# For some reason this code doesn't work not over Q
base = QQ
t = symbols('t')
O = PolynomialRing(base, symbols=[t])
K = base.frac_field(t)

# g = Matrix([
#     [1, t, 1/t],
#     [4, t**2, 1 - t],
#     [1/t**2, 3, 2 * t]
# ])

g2 = Matrix([
    [t ** 2, t ** 3, t],
    [4 * (t ** 2), t ** 4, t ** 2 - t ** 3],
    [1, 3 * (t ** 2), 2 * (t ** 3)]
])


# Cartan decomposition says
# G(K) = G(O) T(K) G(O)
# and T(K) = X_*(T) T(O) so G(O)\G(K)/G(O) = X_*(T)
# But in particular from G(K) we want the elts g', g'' so g = g' diag g''
# But SNF only works if g \in G(O)

def cartan(g):
    """
    g \in G(K) \mapsto (g', t, g'') with g'tg'' = g; g', g'' \in G(O), t \in T(K)
    """

    # clear denominator, do smith normal form, then divide t's diagonal by the original denominator

    t_k, g_o_left, g_o_right = smith_normal_decomp(g, domain=O)

    return g_o_left.inv(), t_k, g_o_right.inv()


def mod_t(matrix: Matrix) -> Matrix:
    return matrix.applyfunc(lambda x: x.subs(t, 0))

def elementary_swap(size, row_i, row_j):
    out = eye(size)
    out[row_i, row_i] = 0
    out[row_j, row_j] = 0
    out[row_i, row_j] = 1
    out[row_j, row_i] = 1

    return out

def elementary_scale(size, row_i, scalar):
    return diag(*[
        1 if i != row_i else scalar for i in range(size)
    ])

def elementary_add(size, row_inp, row_outp, scalar):
    out = eye(size)
    out[row_outp, row_inp] = scalar
    return out

def bruhat(matrix_fd, e):
    # currently only works for borel
    pivot_row = 0
    pivot_col = 0

    size = matrix_fd.shape[0]
    left = eye(size)
    w = eye(size)
    right = matrix_fd

    while pivot_row < size and pivot_col < size:
        i_max = max(
            [(i, right[i, pivot_col]) for i in range(pivot_row, size)], key=lambda x: x[1]
        )[0]
        if right[i_max, pivot_col] == 0:
            pivot_col += 1
        else:
            swap = elementary_swap(size, pivot_row, i_max)
            w *= swap
            right = swap * right

            for i in range(pivot_row + 1, size):
                quotient = right[i, pivot_col] / right[pivot_row, pivot_col]
                subtr = elementary_add(size, pivot_row, i, -quotient)
                right = subtr * right
                left *= subtr.inv()

            pivot_row += 1
            pivot_col += 1

    return left, w, right


def iwahori_bruhat(matrix: Matrix, e):
    """
    G(O) = P W I_P
    :param matrix:
    :return:
    """
    matrix_fd = mod_t(matrix)
    rest_of_matrix = matrix - matrix_fd

    left_p, w, right_p = bruhat(matrix_fd, e)

    right_i = right_p + w.inv() * left_p.inv() * rest_of_matrix

    return left_p, w, right_i

def right_iwahori_decomposition(matrix, e):
    # Step 1: Cartan
    cartan_decomposition = cartan(matrix)
    t_k: Matrix
    g_o_left, t_k, g_o_right = cartan_decomposition
    assert g_o_left * t_k * g_o_right == matrix, "Cartan failed"
    t_k_inverse = t_k.inv()
    # Step 2: Iwahori-Bruhat on right G(O)
    p_left, w, i_right = iwahori_bruhat(g_o_right, e)
    assert p_left * w * i_right == g_o_right, "Iwahori-Bruhat on right factor failed"
    # Step 3: Pass P, W through T
    p_left = t_k * p_left * t_k_inverse
    w = t_k * w * t_k_inverse
    g_o_left *= (p_left * w)
    assert (g_o_left * t_k * i_right).applyfunc(simplify) == matrix, "Final sanity check failed"

    return g_o_left, t_k, i_right


# 1. Cartan decomposition of arbitrary g
# 2. Bruhat decomposition of left G(O) factor
# 3. Sample matrices in U_P(K) and see what Bruhat says

# Note that G(O)/(first congruence subgroup) = G , so G(O)/I = G/B on the nose, so G(O) = B*W*I, and if you have an element g \in G(O), you just reduce it mod t and take the corresponding Bruhat decomposition to get this decomposition. There are similar statements for parabolics, of course.
#
#
# 
# From here probably it's clear: G(K) = G(O)T(K)G(O) = G(O)*T(K)*B*W*I = G(O)*B*W*T(K)*I= G(O) T(K) I.
g = Matrix([
    [1 + t, 2, 3, 4],
    [0, 1, 0, 4],
    [0, 5, 1, 1],
    [0, 0, 0, 1]
])

e = Matrix([[0, 1, 0, 0], [0, 0, 0, 0], [0, 0, 0, 1], [0, 0, 0, 0]])

pprint(right_iwahori_decomposition(g, e))