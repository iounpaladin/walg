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

    return g_o_left, t_k, g_o_right


def iwahori_bruhat(matrix: Matrix) -> Tuple(Matrix, Matrix, Matrix):
    """
    G(O) = BWI
    :param matrix:
    :return:
    """
    m = matrix.rref
    # for testing, just bruhat
    return eye(3), eye(3), eye(3)


def right_iwahori_decomposition(matrix):
    # Step 1: Cartan
    cartan_decomposition = cartan(matrix)
    t_k: Matrix
    g_o_left, t_k, g_o_right = cartan_decomposition
    t_k_inverse = t_k.inv()
    # Step 2: Iwahori-Bruhat on right G(O)
    b_left, w, i_right = iwahori_bruhat(g_o_right)
    # Step 3: Pass B, W through T
    b_left = t_k * b_left * t_k_inverse
    w = t_k * w * t_k_inverse
    g_o_left *= (b_left * w)

    assert g_o_left * t_k * i_right == matrix
    return g_o_left, t_k, i_right


# TODO:
# 1. Cartan decomposition of arbitrary g
# 2. Bruhat decomposition of left G(O) factor
# 3. Sample matrices in U_P(K) and see what Bruhat says

# Note that G(O)/(first congruence subgroup) = G , so G(O)/I = G/B on the nose, so G(O) = B*W*I, and if you have an element g \in G(O), you just reduce it mod t and take the corresponding Bruhat decomposition to get this decomposition. There are similar statements for parabolics, of course.
#
#
# 
# From here probably it's clear: G(K) = G(O)T(K)G(O) = G(O)*T(K)*B*W*I = G(O)*B*W*T(K)*I= G(O) T(K) I.
g = Matrix([
    [1, 2, 3],
    [0, 1, 4],
    [0, 5, 1]
])

pprint(right_iwahori_decomposition(g))