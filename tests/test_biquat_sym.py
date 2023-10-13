import biquaternion_py as bq
import numpy.testing as nt
from biquaternion_py import II, JJ, EE
import sympy as sy

bq.define_algebra()

x1, x2, x3, x4, x5, x6, x7, x8 = sy.symbols("x1 x2 x3 x4 x5 x6 x7 x8")
y1, y2, y3, y4, y5, y6, y7, y8 = sy.symbols("y1 y2 y3 y4 y5 y6 y7 y8")
x = bq.BiQuaternion([x1, x2, x3, x4, x5, x6, x7, x8])
y = bq.BiQuaternion([y1, y2, y3, y4, y5, y6, y7, y8])


def test_multiplication():
    assert (x * y).coeffs == [
        x1 * y1 - x2 * y2 - x3 * y3 - x4 * y4,
        x1 * y2 + x2 * y1 + x3 * y4 - x4 * y3,
        x1 * y3 - x2 * y4 + x3 * y1 + x4 * y2,
        x1 * y4 + x2 * y3 - x3 * y2 + x4 * y1,
        x1 * y5 - x2 * y6 - x3 * y7 - x4 * y8 + x5 * y1 - x6 * y2 - x7 * y3 - x8 * y4,
        x1 * y6 + x2 * y5 + x3 * y8 - x4 * y7 + x5 * y2 + x6 * y1 + x7 * y4 - x8 * y3,
        x1 * y7 - x2 * y8 + x3 * y5 + x4 * y6 + x5 * y3 - x6 * y4 + x7 * y1 + x8 * y2,
        x1 * y8 + x2 * y7 - x3 * y6 + x4 * y5 + x5 * y4 + x6 * y3 - x7 * y2 + x8 * y1,
    ]


def test_division():
    assert sy.simplify(str(x / x)) == 1


def test_addiditon():
    assert x + y == y + x


def test_subtraction():
    assert x - x == 0
    assert x - y == -(y - x)


def test_conjugation():
    assert ~x == bq.BiQuaternion([x1, -x2, -x3, -x4, x5, -x6, -x7, -x8])
    assert x.conjugate() == bq.BiQuaternion([x1, -x2, -x3, -x4, x5, -x6, -x7, -x8])


def test_epsilon_conjugation():
    assert x.eps_conjugate() == bq.BiQuaternion([x1, x2, x3, x4, -x5, -x6, -x7, -x8])


def test_inverse():
    assert sy.simplify(x * x.inv()) == 1
    assert sy.simplify(x.inv() * x) == 1


def test_non_invertible():
    with nt.assert_raises(ValueError):
        bq.BiQuaternion([0, 0, 0, 0, x1, x2, x3, x4]).inv()


def test_scalars():
    assert bq.BiQuaternion(x1) == x1


def test_define_algebra():
    bq.define_algebra(y1, y2, y3)

    assert II * II == y1
    assert JJ * JJ == y2
    assert EE * EE == y3
    bq.define_algebra()
