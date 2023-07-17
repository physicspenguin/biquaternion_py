"""Implementation of polynomial class and associated functions."""

# from sympy.core.symbol import Symbol
from sympy import Poly as sympoly

# from sympy import symbols


def _coeff_gen(poly, indets):
    """Generate coefficient array for a polynomial with given indeterminate.

    Parameters
    ----------
    poly : sympy expression
    Polynomial given.

    indets : list of sympy.core.symbol.Symbol
    Indeterminate for which to generate the coefficients

    """
    print(poly)
    if isinstance(poly, list):
        for i, val in enumerate(poly):
            poly[i] = _coeff_gen(val, indets[1:])
    else:
        return sympoly(poly, indets[0])


class Poly:
    """Class implementing arbitrary polynomials."""

    def __init__(self, indets, poly):
        """Initialize an instance of the polynomial class.

        Parameters
        ----------
        indets : tuple of sympy.core.symbol.Symbol, or strings
        indeterminates used in the polynomial
        poly: sympy.Expr, tuple

        """
        pass

        # # Handling of indeterminates
        # if isinstance(indets, Symbol):
        #     self.indets = [indets]
        # elif isinstance(indets, list):
        #     self.indet = [0 for i, _ in enumerate(indets)]
        #     if isinstance(indets[0], Symbol):
        #         for i, val in enumerate(indets):
        #             self.indet[i] = val
        #     else:
        #         for i, val in enumerate(indets):
        #             sym = Symbol(str(val))
        #             self.indet[i] = sym
        # elif isinstance(indets, str):
        #     self.indets = symbols(indets)

        # else:
        #     raise ValueError(
        #         "indets must be of type, list(sympy Symbols), list(strings), string"
        #     )

        # # Assume that this is now given as a list of coefficients.
        # # Each coefficient can then again be a list, such that multivariate cases
        # # are handled
        # if isinstance(poly, list):
        #     self.coeff = poly
        # # This assumes a sympy expression
        # # (Allthough most likely it will then be read as a biquaternion class)
        # else:
        #     coeffs = sympoly(poly, self.indets[0]).all_coeffs()
        #     self.coeff = _coeff_gen(coeffs, indets[1:])
