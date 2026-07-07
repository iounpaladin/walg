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
    [sympify(t**2), t ** 3, t],
    [4 * (t ** 2), t**4, t ** 2 - t ** 3],
    [1, 3 * (t ** 2), 2 * (t ** 3)]
])

# Cartan decomposition says
# G(K) = G(O) T(K) G(O)
# and T(K) = X_*(T) T(O) so G(O)\G(K)/G(O) = X_*(T)
# But in particular from G(K) we want the elts g', g'' so g = g' diag g''
# But SNF only works if g \in G(O)

pprint(smith_normal_decomp(g2, domain=O))

# TODO:
# 1. Cartan decomposition of arbitrary g
# 2. Bruhat decomposition of left G(O) factor
# 3. Sample matrices in U_P(K) and see what Bruhat says