import biquaternion_py.polynomials as bp
from biquaternion_py import Poly
import sympy as sy

t, s = sy.symbols("t s")
a = sy.symbols("a:8")
b = sy.symbols("b:8")


def test_max_power():
    assert 0 == bp._max_pow(3.2, t)


def test_all_indet_coeffs():
    assert a == tuple(
        Poly(sum([b * t**i for i, b in enumerate(a)]), t).all_indet_coeffs(t)
    )


def test_all_coeffs():
    p = Poly(
        a[0] * t**2 * s + a[1] * t**1 + a[2] * s**2 + a[3] * s + a[4], [t, s]
    )
    q = Poly(sum([b * t**i for i, b in enumerate(a)]), t)
    assert [[a[4], a[3], a[2]], [a[1]], [0, a[0]]] == bp._all_coeffs(p, [t, s])
    assert a == tuple(q.all_coeffs())
    assert [
        [q.poly],
        [a[1] * t**1 + a[4], a[3] + a[0] * t**2, a[2]],
    ] == bp._all_coeffs([q, p], [t, s])


def test_terms():
    q = Poly(sum([b * t**i for i, b in enumerate(a)]), t)
    p = Poly(
        a[0] * t**2 * s + a[1] * t**1 + a[2] * s**2 + a[3] * s + a[4], [t, s]
    )
    out = [
        ((7,), a[7]),
        ((6,), a[6]),
        ((5,), a[5]),
        ((4,), a[4]),
        ((3,), a[3]),
        ((2,), a[2]),
        ((1,), a[1]),
        ((0,), a[0]),
    ]
    out2 = [
        ((2, 1), a[0]),
        ((1, 0), a[1]),
        ((0, 2), a[2]),
        ((0, 1), a[3]),
        ((0, 0), a[4]),
    ]

    assert out == bp._terms(q)
    assert out2 == bp._terms(p)


def test_eval_poly():
    q = Poly(sum([x * t**i for i, x in enumerate(a)]), t)
    p = Poly(sum([x * t**i for i, x in enumerate(b)]), t)
    assert sum(a) == q.eval(1, True)
    assert sum(a) == q.eval(1, False)
    assert sum(a) == q.eval((1), True)
    assert sum(a) + sum(b) == (p + q).eval([1, 1], False)
    assert sum(a) + sum(b) == (p + q).eval((1, 1), False)
