from __init__ import *
from sympy import symbols

t = symbols('t')  # Field = C((t))


def create_symbols(name, number):
    """
    Creates sympy symbols for a range of variables to be indexed in the form n_1, n_2, ...

    @:param name: The name of the variables
    @:param number: The number of variables to create
    @:return: A list of sympy symbols corresponding to the variables
    """
    return create_symbols_list(name, map(lambda x: x + 1, range(number)))


def create_symbols_list(name, list_of_numbers):
    """
    Creates sympy symbols for a range of variables to be indexed in a more general form

    @:param name: The name of the variables
    @:param list_of_numbers: The list of numbered variables
    @:return: A list of sympy symbols corresponding to the variables
        """
    return symbols(' '.join([
        f'{name}_{i}' for i in list_of_numbers
    ]))


def formal_power_series(coeff_name, t_power=bounded_t_power):
    """
    Creates a general formal power series with named coefficients and a bounded degree.
    Calls :ref:`formal_laurent_series` internally.
    :param coeff_name: Name of coefficients
    :param t_power: Bounded degree
    :return: Expression of a formal power series
    """
    return formal_laurent_series(coeff_name, t_power, 0)


def formal_laurent_series(coeff_name, t_power=bounded_t_power, t_power_below=bounded_t_power_below):
    """
    Creates a general formal Laurent series with named coefficients and a bounded degree
    :param coeff_name: Name of coefficients
    :param t_power: Bounded degree
    :param t_power_below: Bounded negative degree
    :return: Expression of a formal Laurent series
    """
    numbers = list(range(-t_power_below, t_power + 1))
    coefficients = zip(numbers, create_symbols_list(coeff_name, numbers))

    return sum([
        coeff * (t ** i) for i, coeff in coefficients
    ])
