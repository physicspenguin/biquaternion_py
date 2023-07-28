"""Extra functions for polynomials."""

from .biquaternion import BiQuaternion
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
