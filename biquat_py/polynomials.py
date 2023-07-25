"""Implementation of polynomial class and associated functions."""

# from sympy.core.symbol import Symbol
from sympy import expand, Pow, Expr, sympify

# from .biquaternion import BiQuaternion


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

    def deg(self, var):
        return _max_pow(self.poly, var)

    def __pos__(self):
        return Poly(self.poly, *self.indets)

    def __neg__(self):
        return Poly(-self.poly, *self.indets)

    def __mul__(self, other):
        if isinstance(other, Poly):
            return Poly(
                expand(self.poly * other.poly),
                *self.indets,
                *set(other.indets).difference(set(self.indets)),
            )
        else:
            return Poly(expand(self.poly * sympify(other)), *self.indets)

    def __rmul__(self, other):
        if isinstance(other, Poly):
            return Poly(
                expand(other.poly * self.poly),
                *self.indets,
                *set(other.indets).difference(set(self.indets)),
            )
        else:
            return Poly(expand(sympify(other) * self.poly), *self.indets)

    def __add__(self, other):
        if isinstance(other, Poly):
            return Poly(
                expand(self.poly + other.poly),
                *self.indets,
                *set(other.indets).difference(set(self.indets)),
            )
        else:
            return Poly(expand(self.poly + sympify(other)), *self.indets)

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


def poly_div(poly_1, poly_2, var, right=True):
    """Polynomial division with remainder of poly_1 and poly_2 with respect to var.

    Parameters
    ----------
    poly_1 : Poly
        Polynomial which should be divided.
    poly_2 : Poly
        Polynomial by which should be divided.
    var : sympy.Symbol
        Variable with respect to which to divide.
    right : (optional, default = True) bool
        Should right division be used.

    Returns
    -------
    quotient : Poly
        Result of division
    remainder : Poly
        Remainder of division

    Notes
    -----
    This function produces polynomials `quotient` and `remainder` such that
    poly_1 = poly_2 * quotient + remainder for `right = True`
    poly_1 = quotient * poly_2 + remainder for `right = False`
    """

    init_lead_coeff = poly_2.lcoeff(var)

    if right:
        f0 = init_lead_coeff * poly_1
        g0 = init_lead_coeff * poly_2
    else:
        f0 = poly_1 * init_lead_coeff
        g0 = poly_2 * init_lead_coeff

    quotient = Poly(0, var)
    remainder = f0

    m = remainder.deg(var)
    n = poly_2.deg(var)

    while m >= n:
        lead_coeff = remainder.lcoeff(var)
        quotient = quotient + (lead_coeff * (var ** (m - n)))
        if right:
            remainder = expand(remainder - (g0 * lead_coeff * var ** (m - n)))
        else:
            remainder = expand(remainder - (lead_coeff * var ** (m - n) * g0))

        m = remainder.deg(var)

    if right:
        return quotient, init_lead_coeff * remainder
    else:
        return quotient, remainder * init_lead_coeff
