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

    for i in range(4):
        gcd = sy.gcd(gcd, polys[i])
        gcd = sy.gcd(gcd, polys[i + 4])
    return gcd
