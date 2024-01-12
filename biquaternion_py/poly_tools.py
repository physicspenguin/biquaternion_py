"""Extra functions for polynomials."""

from .biquaternion import BiQuaternion
from .polynomials import poly_div, Poly
import sympy as sy


def max_real_poly_fact(poly):
    """Calculate maximal real polynomial factor of the BiQuaternionpolynomial `poly`.

    Parameters
    ----------
    poly : Poly
        Polynomial of which to find the maximal real factor.

    Returns
    -------
    gcd : Poly
        Maximal real factor of `poly`
    """
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
    """Real gcd of c, primal*dual.conjugate(), primal.conjugate()*dual.

    Parameters
    ----------
    poly : Poly
        Polynomial of which to find the gcd of maximal real poly factor of the
    primal part, primal*dual.conjugate(), primal.conjugate()*dual

    Returns
    -------
    gcd : Poly

    """
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


def irreducible_factors(poly, domain=None):
    """Calculate the irreducible factors of a polynomial.

    Parameters
    ----------
    poly : Poly
        Polynomial of which to calculate the irreducible factors.
    domain : string, optional
        Domain over which to calculate the irreducible factors.
        (Default None lets sympy decide which domain to use.)

    Returns
    -------
    out : array of Poly
        List of irreducible factors.

    Notes
    -----
    If the irreducible factors are not calculated correctly this might be an issue
    of sympy assuming to much about the domain. Chaning this to "QQ" or "RR"
    """
    var = poly.indets[0]
    t = sy.Symbol(var.name, real=True)
    poly1 = Poly(poly.poly.subs({var: t}), t)
    if domain:
        factors = sy.polys.polyroots.root_factors(poly1.poly, t, domain=domain)
    else:
        factors = sy.polys.polyroots.root_factors(poly1.poly, t)
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
    """Split off linear factor with norm `norm` from poly.

    Parameters
    ----------
    poly : Poly
        Polynomial of which to split of a linear factor
    norm : Poly
        Norm of the linear factor to split off

    Returns
    -------
    quot : Poly
        Quotient of the polynomial left division of `poly` by `norm`
    lin_fact : Poly
        Linearfactor that was split off corresponding to `norm`
    """
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
    """Factorize polynomial given a list of factors of the norm polynomial.

    Parameters
    ----------
    poly : Poly
        Polynomial which to factorize
    factors : array of Poly
        Array of irreducible factors of the norm polynomial of `poly`

    Returns
    -------
    out : array of Poly
        Array of linear factors of Poly corresponding to the order of factors
    """
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


def factorize_bq_poly(poly, domain=None):
    """Factorize Biquaternion polynomial into linear factors.

    Parameters
    ----------
    poly : Poly
        Polynomial which to factorize
    domain : string, optional
        Domain over which to calculate the irreducible factors.
        (Default None lets sympy decide which domain to use.)
    Returns
    -------
    factors : array of Poly
        Array of linear factors of `poly` associated to the order of the factors
    given by irreducible_factors.
    """
    norm = poly.norm()
    # if not is_poly_real(norm):
    #     raise ValueError("Norm must be a real polynomial.")
    norm = Poly(norm.poly.scal, *norm.indets)
    _, factors = irreducible_factors(norm, domain)
    return factorize_from_list(poly, factors)
