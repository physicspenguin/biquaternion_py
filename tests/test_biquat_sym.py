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


def test_too_many_args():
    with nt.assert_raises(ValueError):
        bq.BiQuaternion([1, 2, 3, 4, 5, 6, 7, 8, 9])


def test_non_commutative_args():
    with nt.assert_raises(ValueError):
        bq.BiQuaternion([1, bq.BiQuaternion([1, 3, 4])])


def test_scal_setter():
    b = bq.BiQuaternion()
    b.scal = 4
    assert b.scal == 4


def test_i_setter():
    b = bq.BiQuaternion()
    b.i = 4
    assert b.i == 4


def test_j_setter():
    b = bq.BiQuaternion()
    b.j = 4
    assert b.j == 4


def test_k_setter():
    b = bq.BiQuaternion()
    b.k = 4
    assert b.k == 4


def test_eps_setter():
    b = bq.BiQuaternion()
    b.eps = 4
    assert b.eps == 4


def test_ei_setter():
    b = bq.BiQuaternion()
    b.ei = 4
    assert b.ei == 4


def test_ej_setter():
    b = bq.BiQuaternion()
    b.ej = 4
    assert b.ej == 4


def test_ek_setter():
    b = bq.BiQuaternion()
    b.ek = 4
    assert b.ek == 4


def test_coeffs_setter():
    b = bq.BiQuaternion()
    b.coeffs = [1, 2, 3, 4, 5, 6, 7, 8]
    assert b.coeffs == [1, 2, 3, 4, 5, 6, 7, 8]


def test_pos():
    b = bq.BiQuaternion()
    b.coeffs = [1, 2, 3, 4, 5, 6, 7, 8]
    assert b == +b


def test_poly_bq_add():
    t = sy.symbols("t")
    p = bq.Poly((t**2 + bq.II * t + bq.KK), t)
    b = bq.BiQuaternion([1, 2, 3, 4, 5, 6, 7, 8])
    assert b + p == bq.Poly((t**2 + bq.II * t + bq.KK + b), [t])


def test_repr():
    b = bq.BiQuaternion([1, 2, 3, 4, 5, 6, 7, 8])

    assert (
        repr(b)
        == "(( 1 ) + ( 2 ) * II + ( 3 ) * JJ + ( 4 ) * KK) +"
        + " EE * (( 5 ) + ( 6 ) * II + ( 7 ) * JJ + ( 8 ) * KK)"
    )


def test_powers():
    b = bq.BiQuaternion([1, 2, 3, 4, 5, 6, 7, 8])
    assert b**3 == b * b * b
    assert b ** (-3) == 1 / b / b / b
    with nt.assert_raises(TypeError):
        b ** (0.34)


def test_quadrance():
    b = bq.BiQuaternion([1, 2, 3, 4, 5, 6, 7, 8])
    assert b.quadrance() == 140 * bq.EE + 30
    assert b.norm() == 140 * bq.EE + 30


def test_reverse_add():
    a = bq.BiQuaternion([1, 1, 1, 1, 1, 1, 1, 1])
    b = a + 4

    assert bq.BiQuaternion([5, 1, 1, 1, 1, 1, 1, 1]) == b


def test_reverse_div():
    a = bq.BiQuaternion([1, 1, 1, 1, 1, 1, 1, 1])
    b = a * 4

    assert b / 4 == a


def test_parts():
    b = bq.BiQuaternion([1, 2, 3, 4, 5, 6, 7, 8])

    assert b.primal() == bq.BiQuaternion([1, 2, 3, 4])
    assert b.dual() == bq.BiQuaternion([5, 6, 7, 8])
    assert b.scalar_part() == bq.BiQuaternion([1, 0, 0, 0, 5, 0, 0, 0])
    assert b.vector_part() == bq.BiQuaternion([0, 2, 3, 4, 0, 6, 7, 8])


def test_elementwise():
    b = bq.BiQuaternion([1, 2, 3, 4, 5, 6, 7, 8])
    b = b.apply_elementwise(lambda x: x**2)
    assert b == bq.BiQuaternion([1, 4, 9, 16, 25, 36, 49, 64])
