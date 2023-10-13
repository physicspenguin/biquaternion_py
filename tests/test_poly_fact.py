import biquaternion_py as bq
import sympy as sy


def test_factorization():
    t = sy.symbols("t")
    h1 = bq.rand_rational() + bq.rand_line()
    h2 = bq.rand_rational() + bq.rand_line()
    h3 = bq.rand_rational() + bq.rand_line()

    poly = bq.Poly((t - h1) * (t - h2) * (t - h3), t)
    norm_poly = bq.Poly(poly.norm().poly.scal, *poly.indets)
    _, facts = bq.irreducible_factors(norm_poly)
    fact1 = bq.factorize_from_list(poly, facts)
    fact2 = bq.factorize_from_list(poly, facts[::-1])
    poly1 = 1
    for fac in fact1:
        poly1 *= fac
    poly2 = 1
    for fac in fact2:
        poly2 *= fac

    print(poly == poly1 == poly2)
    assert poly == poly1 == poly2
