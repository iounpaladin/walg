from sympy import *
from sympy.interactive.printing import init_printing
from sympy.matrices import *

init_printing(use_latex=False)

t = symbols('t')  # Field = C((t))
dimension = 2
bounded_t_power = 3
bounded_t_power_below = bounded_t_power


def create_symbols(name, number):
    return create_symbols_list(name, map(lambda x: x + 1, range(number)))


def create_symbols_list(name, list_of_numbers):
    return symbols(' '.join([
        f'{name}_{i}' for i in list_of_numbers
    ]))


def formal_power_series(coeff_name, t_power=bounded_t_power):
    return formal_laurent_series(coeff_name, t_power, 0)


def formal_laurent_series(coeff_name, t_power=bounded_t_power, t_power_below=bounded_t_power_below):
    numbers = list(range(-t_power_below, t_power + 1))
    coefficients = zip(numbers, create_symbols_list(coeff_name, numbers))

    return sum([
        coeff * (t ** i) for i, coeff in coefficients
    ])


# Code goal: given cocharacter lambda, determine if i lambda(t) g in N(k) has a solution with i in I_p, g in G(O)
# with given bounds on the powers of t appearing in i and g

def cocharacter(size=dimension):  # Generic cocharacter
    ns = create_symbols('λ', size)

    return diag(*[t ** ns[i] for i in range(size)])


def iwahori(size=dimension):
    pass


def parabolic(size=dimension):
    pass


def g_o(size=dimension):
    pass

pprint(cocharacter())
pprint(formal_laurent_series('N'))
