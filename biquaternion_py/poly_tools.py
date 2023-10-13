"""Extra functions for polynomials."""

from .biquaternion import BiQuaternion
from .polynomials import poly_div, Poly
import sympy as sy


def max_real_poly_fact(poly):
    """Calculate maximal real polynomial factor of the BiQuaternionpolynomial `poly`."""
    if len(poly.indets) != 1:
        raise ValueError("Only univariate polynomials are supported.")
    if not isinstance(poly.poly, BiQuaternion):
        raise ValueError(
            "Only polynomials with coefficients in "
            + "the BiQuaternions are supported."
        )

    polys = [val for val in poly.poly.coeffs]

    gcd = 0

    for val in polys:
        gcd = sy.gcd(gcd, val)
    return gcd


def gcd_conj_pd(poly):
    """Real gcd of c, primal*dual.conjugate(), primal.conjugate()*dual."""
    c = max_real_poly_fact(poly.primal())
    primal = (poly * (1 / c)).poly.primal()
    dual = poly.poly.dual()

    primal_dual_conj = (primal * dual.conjugate()).coeffs[0:4]
    primal_conj_dual = (primal.conjugate() * dual).coeffs[0:4]

    gcd = sy.gcd(c, primal_dual_conj[0])

    for i in range(3):
        gcd = sy.gcd(gcd, primal_dual_conj[i + 1])

    for val in primal_conj_dual:
        gcd = sy.gcd(gcd, val)

    return gcd


def is_poly_reduced(poly):
    """Check if polynomial is reduced.

    Parameters
    ----------
    poly : Poly
        Polynomial which to check for reducedness.

    Returns
    -------
    bool
        True if the polynomial is reduced.

    Notes
    -----
    A polynomial is called reduced, if the primal and dual part have no common
    real factor. [1]_

    .. [1] Z. Li, J. Schicho, H.-P. Schröcker,
       The rational motion of minimal dual quaternion degree with prescribed trajectory,
       Computer Aided Geometric Design,
       Volume 41,
       2016,
       Pages 1-9,
       ISSN 0167-8396,
       https://doi.org/10.1016/j.cagd.2015.10.002.
    """
    return (
        sy.gcd(max_real_poly_fact(poly.primal()), max_real_poly_fact(poly.dual())) == 1
    )


# def is_poly_real(poly):
#     return np.allclose(poly.poly.coeffs[1:], np.zeros(7))


def irreducible_factors(poly):
    """Calculate the irreducible factors of a polynomial."""
    # if not is_poly_real(poly):
    #     raise ValueError("Polynomial must have real coefficients.")
    var = poly.indets[0]
    t = sy.Symbol(var.name, real=True)
    poly1 = Poly(poly.poly.subs({var: t}), t)
    factors = sy.polys.polyroots.root_factors(poly1.poly)
    out = []
    for i, val in enumerate(factors):
        if val.is_real:
            out = out + [val]
        else:
            for j, val2 in enumerate(factors[i:]):
                if val == val2.conjugate():
                    out = out + [sy.expand(val * val2)]
                    factors.pop(i + j)

    for i, val in enumerate(out):
        out[i] = Poly(val.subs({t: var}), var)
    return poly.lcoeff(var), out


def split_lin_factor(poly, norm):
    """Split off linear factor with norm `norm` from poly."""
    if len(poly.indets) != 1:
        raise ValueError("Only univariate polynomials supported.")
    if poly.indets != norm.indets:
        raise ValueError("Poly and norm must have the same indeterminates.")

    indet = poly.indets[0]
    _, rem = poly_div(poly, norm, indet, False)
    root = -(rem.poly.coeff(indet, 1).inv()) * rem.poly.coeff(indet, 0)
    lin_fact = Poly(indet - root, indet)
    quot, _ = poly_div(poly, lin_fact, indet, False)
    return quot, lin_fact


def factorize_from_list(poly, factors):
    """Factorize polynomial given a list of factors of the norm polynomial."""
    if len(poly.indets) != 1:
        raise ValueError("Only univariate polynomials supported.")
    out = []
    poly0 = poly
    for i, val in enumerate(factors[::-1]):
        if poly.indets != val.indets:
            raise ValueError("Poly and factor must have the same indeterminates.")
        poly0, lin_fact = split_lin_factor(poly0, val)
        out = [lin_fact] + out
    return out


def factorize_bq_poly(poly):
    """Factorize Biquaternion polynomial into linear factors."""
    norm = poly.norm()
    # if not is_poly_real(norm):
    #     raise ValueError("Norm must be a real polynomial.")
    norm = Poly(norm.poly.scal, *norm.indets)
    _, factors = irreducible_factors(norm)
    return factorize_from_list(poly, factors)
