"""Implementation of polynomial class and associated functions."""

# from sympy.core.symbol import Symbol
from sympy import Poly as sympoly
from sympy import expand, Pow, Expr, sympify

# from .biquaternion import BiQuaternion


def _coeff_gen(poly, indets):
    """
    Generate coefficient array for a polynomial with given indeterminate.

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


def _max_pow(expr, indet):
    """Find the maximal power of indet in expr.

    Parameters
    ----------
    expr : sympy expression
        Expression of which to find the maximal power of indet.
    indet : sympy symbol
        Indeterminate of which to find the maximal power in expr

    Returns
    -------
    int
        Maximal power of indet in expr
    """
    deg = max(
        (
            z.as_base_exp()[1] if z.as_base_exp()[0] == indet else 0
            for z in expand(expr).atoms(Pow)
        ),
        default=1,
    )
    if deg <= 1:
        if expr.coeff(indet, 1) != 0:
            return 1
        elif expr.coeff(indet, 0) != 0:
            return 0
        else:
            return deg
    return deg


def _all_indet_coeffs(expr, indet):
    """Generate list of all coefficients of expr with respect to indet.

    Parameters
    ----------
    expr : sympy.Expr
        Expression of which to generate list of coefficients
    indet : sympy.Symbol
        Indeterminate of which to find list of coefficients

    Returns
    -------
    List
        List containing coefficients of powers of indet in descending order.
    """
    deg = _max_pow(expr, indet)

    coeffs = [expr.coeff(indet, n) for n in range(deg + 1)]
    return coeffs


def _all_coeffs(expr, indets):
    """Generate list of all coefficients of expr ordered same as indets.

    Parameters
    ----------
    expr : sympy.Expr
        Expression of which to generate list of coefficients
    indets : sympy.Symbol
        Indeterminate of which to find list of coefficients

    Returns
    -------
    List
        List containing coefficients of powers of indet in descending order.
    """
    pass


class Poly(Expr):
    """Class implementing arbitrary polynomials."""

    _op_priority = 12.1

    def __new__(cls, *args):
        """Create an instance of the polynomial class.

        Parameters
        ----------
        indets : tuple of sympy.core.symbol.Symbol, or strings
        indeterminates used in the polynomial
        poly: sympy.Expr, tuple

        """
        poly, *indets = args
        poly, *indets = map(sympify, (poly, *indets))
        obj = Expr.__new__(cls, poly, *indets)
        obj._poly = poly
        obj._indets = indets
        return obj

    @property
    def poly(self):
        return self._poly

    @property
    def indets(self):
        return self._indets

    def __pos__(self):
        return Poly(self.poly, self.indets)

    def __neg__(self):
        return Poly(-self.poly, self.indets)

    def __mul__(self, other):
        if isinstance(other, Poly):
            return Poly(
                expand(self.poly * other.poly),
                *self.indets,
                *set(other.indets).difference(set(self.indets)),
            )
        else:
            return Poly(expand(self.poly * sympify(other)), self.indets)

    def __rmul__(self, other):
        if isinstance(other, Poly):
            return Poly(
                expand(other.poly * self.poly),
                *self.indets,
                *set(other.indets).difference(set(self.indets)),
            )
        else:
            return Poly(expand(sympify(other) * self.poly), self.indets)

    def __add__(self, other):
        if isinstance(other, Poly):
            return Poly(
                expand(self.poly + other.poly),
                *self.indets,
                *set(other.indets).difference(set(self.indets)),
            )
        else:
            return Poly(expand(self.poly + sympify(other)), self.indets)

    __radd__ = __add__

    def __sub__(self, other):
        return self + (-other)

    def __rsub__(self, other):
        return other + (-self)

    def __repr__(self):
        return f"Poly({repr(self.poly)},{repr(self.indets)})"

    def __str__(self):
        return f"Poly({self.poly},{self.indets})"

    def __eq__(self, other):
        return self.poly == other.poly and set(self.indets) == set(other.indets)

    __hash__ = super.__hash__

    def coeff(self, var, power=1, right=False, _first=True):
        return expand(self.poly).coeff(var, power, right, _first)

    def lcoeff(self, var):
        return expand(self.poly).coeff(var, _max_pow(self.poly, var))

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

    def all_coeffs(self, var):
        """Compute all coefficients with respect to var.

        Parameters
        ----------
        var : sympy.Symbol
            Variable with respect which to computa all coefficients.

        Returns
        -------
        list of BiQuaternions
            List of coefficients for powers of var in descending order.
        """
