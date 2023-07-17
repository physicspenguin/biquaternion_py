"""
This module implements the biquaternion class, which takes care of
 the basic biquaternion arithmetic.
"""

import numpy as np
from sympy.core.expr import Expr
from sympy import sympify

_BQ_I = -1
_BQ_J = -1
_BQ_E = 0
_BACKEND = "general"


def define_backend(backend="general"):
    """Defines the backend to be used.
    You can choose to use a general implementation that supports
    symbolic computations, or numpy for fast numerical computations

    Parameters
    ----------

    backend : string, optional
        values "general" and "numpy" are supported

    Returns
    -------
    None

    """
    global _BACKEND
    _BACKEND = backend


def define_algebra(i_square=-1, j_square=-1, e_square=0):
    """ Defines the algebra

    Parameters
    ----------
    i_square : float
        value of II**2
    j_square : float
        value of JJ**2
    e_square : float
        value of EE**2

    Returns
    -------
    None

    Notes
    -----
    The algebra of bi-quaternions is fully defined by the following conditions:
    .. math:: i^2 = a, j^2 = b, e^2 = c,\\
    ij = k, ji = -k,\\
    ei = ie, ej = je, ek = ke
    All objects commute with elements of the chosen base field.
    Mostly this is :math:`\mathbb{R}` or :math:`\mathbb{C}`.
    """
    global _BQ_I, _BQ_J, _BQ_E, II, JJ, KK, EE
    _BQ_I = i_square
    _BQ_J = j_square
    _BQ_E = e_square
    II = BiQuaternion(*[0, 1, 0, 0, 0, 0, 0, 0])
    JJ = BiQuaternion(*[0, 0, 1, 0, 0, 0, 0, 0])
    KK = BiQuaternion(*[0, 0, 0, 1, 0, 0, 0, 0])
    EE = BiQuaternion(*[0, 0, 0, 0, 1, 0, 0, 0])


def _sanitize_args(*args):
    """Sanitizes the input of the __new__ method of biquaternion"""
    coeffs = [0, 0, 0, 0, 0, 0, 0, 0]
    if len(args) == 1:
        gen = args[0]
        if isinstance(gen, BiQuaternion):
            cof = gen.coeffs
        elif isinstance(gen, (list, tuple, np.ndarray)):
            if len(gen) >= 9:
                raise ValueError("Maximum array length is 8")
            cof = gen
        else:
            cof = [gen]
    else:
        cof = args

    for i, val in enumerate(cof):
        coeffs[i] = val

    return coeffs


class BiQuaternion(Expr):
    """Biquaternions as $a + II b + JJ c + KK d + EE (w + II x + JJ y + KK z)$"""

    is_commutative = False
    _op_priority = 11.5

    def __new__(cls, *args):
        # Sanitize the arguments

        scal, i, j, k, eps, ei, ej, ek = _sanitize_args(*args)
        scal, i, j, k, eps, ei, ej, ek = map(sympify, (scal, i, j, k, eps, ei, ej, ek))

        if any(i.is_commutative is False for i in [scal, i, j, k, eps, ei, ej, ek]):
            raise ValueError("arguments have to be commutative")
        else:
            obj = Expr.__new__(cls, scal, i, j, k, eps, ei, ej, ek)
            obj._scal = scal
            obj._i = i
            obj._j = j
            obj._k = k
            obj._eps = eps
            obj._ei = ei
            obj._ej = ej
            obj._ek = ek
            return obj

    # def __init__(self, *args):
    #     # if _BACKEND == "numpy":
    #     #     self._coeff = np.zeros(8)
    #     # else:
    #     #     self._coeff = [0, 0, 0, 0, 0, 0, 0, 0]
    #     # if isinstance(gen, BiQuaternion):
    #     #     cof = gen._coeff
    #     # elif isinstance(gen, (list, tuple, np.ndarray)):
    #     #     if len(gen) >= 9:
    #     #         raise ValueError("Maximum array length is 8")
    #     #     cof = gen
    #     # else:
    #     #     cof = [gen]

    #     # for i, val in enumerate(cof):
    #     #     self._coeff[i] = val

    @property
    def scal(self):
        return self._scal

    @property
    def i(self):
        return self._i

    @property
    def j(self):
        return self._j

    @property
    def k(self):
        return self._k

    @property
    def eps(self):
        return self._eps

    @property
    def ei(self):
        return self._ei

    @property
    def ej(self):
        return self._ej

    @property
    def ek(self):
        return self._ek

    @property
    def coeffs(self):
        return [
            self.scal,
            self.i,
            self.j,
            self.k,
            self.eps,
            self.ei,
            self.ej,
            self.ek,
        ]

    def __mul__(self, other):
        if isinstance(other, BiQuaternion):
            out = [
                -_BQ_E * _BQ_I * _BQ_J * other.coeffs[7] * self.coeffs[7]
                + _BQ_E * _BQ_I * other.coeffs[5] * self.coeffs[5]
                + _BQ_E * _BQ_J * other.coeffs[6] * self.coeffs[6]
                - _BQ_I * _BQ_J * other.coeffs[3] * self.coeffs[3]
                + _BQ_E * other.coeffs[4] * self.coeffs[4]
                + _BQ_I * other.coeffs[1] * self.coeffs[1]
                + _BQ_J * other.coeffs[2] * self.coeffs[2]
                + other.coeffs[0] * self.coeffs[0],
                _BQ_E * _BQ_J * other.coeffs[6] * self.coeffs[7]
                - _BQ_E * _BQ_J * other.coeffs[7] * self.coeffs[6]
                + _BQ_E * other.coeffs[4] * self.coeffs[5]
                + _BQ_E * other.coeffs[5] * self.coeffs[4]
                + _BQ_J * other.coeffs[2] * self.coeffs[3]
                - _BQ_J * other.coeffs[3] * self.coeffs[2]
                + other.coeffs[0] * self.coeffs[1]
                + other.coeffs[1] * self.coeffs[0],
                -_BQ_E * _BQ_I * other.coeffs[5] * self.coeffs[7]
                + _BQ_E * _BQ_I * other.coeffs[7] * self.coeffs[5]
                + _BQ_E * other.coeffs[4] * self.coeffs[6]
                + _BQ_E * other.coeffs[6] * self.coeffs[4]
                - _BQ_I * other.coeffs[1] * self.coeffs[3]
                + _BQ_I * other.coeffs[3] * self.coeffs[1]
                + other.coeffs[0] * self.coeffs[2]
                + other.coeffs[2] * self.coeffs[0],
                _BQ_E * other.coeffs[4] * self.coeffs[7]
                - _BQ_E * other.coeffs[5] * self.coeffs[6]
                + _BQ_E * other.coeffs[6] * self.coeffs[5]
                + _BQ_E * other.coeffs[7] * self.coeffs[4]
                + other.coeffs[0] * self.coeffs[3]
                - other.coeffs[1] * self.coeffs[2]
                + other.coeffs[2] * self.coeffs[1]
                + other.coeffs[3] * self.coeffs[0],
                -_BQ_I * _BQ_J * other.coeffs[3] * self.coeffs[7]
                - _BQ_I * _BQ_J * other.coeffs[7] * self.coeffs[3]
                + _BQ_I * other.coeffs[1] * self.coeffs[5]
                + _BQ_I * other.coeffs[5] * self.coeffs[1]
                + _BQ_J * other.coeffs[2] * self.coeffs[6]
                + _BQ_J * other.coeffs[6] * self.coeffs[2]
                + other.coeffs[0] * self.coeffs[4]
                + other.coeffs[4] * self.coeffs[0],
                _BQ_J * other.coeffs[2] * self.coeffs[7]
                - _BQ_J * other.coeffs[3] * self.coeffs[6]
                + _BQ_J * other.coeffs[6] * self.coeffs[3]
                - _BQ_J * other.coeffs[7] * self.coeffs[2]
                + other.coeffs[0] * self.coeffs[5]
                + other.coeffs[1] * self.coeffs[4]
                + other.coeffs[4] * self.coeffs[1]
                + other.coeffs[5] * self.coeffs[0],
                -_BQ_I * other.coeffs[1] * self.coeffs[7]
                + _BQ_I * other.coeffs[3] * self.coeffs[5]
                - _BQ_I * other.coeffs[5] * self.coeffs[3]
                + _BQ_I * other.coeffs[7] * self.coeffs[1]
                + other.coeffs[0] * self.coeffs[6]
                + other.coeffs[2] * self.coeffs[4]
                + other.coeffs[4] * self.coeffs[2]
                + other.coeffs[6] * self.coeffs[0],
                other.coeffs[0] * self.coeffs[7]
                - other.coeffs[1] * self.coeffs[6]
                + other.coeffs[2] * self.coeffs[5]
                + other.coeffs[3] * self.coeffs[4]
                + other.coeffs[4] * self.coeffs[3]
                - other.coeffs[5] * self.coeffs[2]
                + other.coeffs[6] * self.coeffs[1]
                + other.coeffs[7] * self.coeffs[0],
            ]
            return BiQuaternion(*out)
            # return BiQuaternion(
            #     out[0], out[1], out[2], out[3], out[4], out[5], out[6], out[7]
            # )

        return self * BiQuaternion(other)

    def __pos__(self):
        """Positive of itsself
        Returns
        -------
        BiQuaternion
            self
        """
        return BiQuaternion(self)

    def __neg__(self):
        """Negative of itsself
        Returns
        -------
        BiQuaternion
            -self
        """
        return BiQuaternion(*[-self.coeffs[i] for i in range(8)])

    def __add__(self, other):
        """Addition of two biquaternions
        Parameters
        ----------
        self: BiQuaternion

        other: BiQuaternion, float, symbolic

        Returns
        -------
        BiQuaternion
            Sum of self and input parameter
        """

        if isinstance(other, BiQuaternion):
            out = [self.coeffs[i] + other.coeffs[i] for i in range(8)]
            return BiQuaternion(*out)

        return self + BiQuaternion(other)

    def __sub__(self, other):
        """subtraction of two biquaternions
        Parameters
        ----------
        self: BiQuaternion

        other: BiQuaternion, float, symbolic

        Returns
        -------
        BiQuaternion
            Difference of self and input parameter
        """
        return self + (-other)

    def __rsub__(self, other):
        """Addition of two biquaternions
        Parameters
        ----------
        self: BiQuaternion

        other: BiQuaternion, float, symbolic

        Returns
        -------
        BiQuaternion
            Difference of input parameter and self
        """
        return other + (-self)

    __radd__ = __add__
    __rmul__ = __mul__

    def __eq__(self, other):
        """Test equality of two biquaternions"""
        othercoeff = BiQuaternion(other).coeffs
        for i, val in enumerate(self.coeffs):
            if val != othercoeff[i]:
                return False
        return True

    __hash__ = super.__hash__

    def __repr__(self):
        """Representation function.

        Parameters
        ----------
        self: self

        Returns
        -------
        result : string
            String representation of biquaternion that can be used to reproduce the
            exact object.
        """

        result = (
            "(" + "( " + repr(self.coeffs[0]) + " )" + " + "
            "( " + repr(self.coeffs[1]) + " )" + " * II" + " + "
            "( " + repr(self.coeffs[2]) + " )" + " * JJ" + " + "
            "( " + repr(self.coeffs[3]) + " )" + " * KK"
            ") + EE * ("
            "( " + repr(self.coeffs[4]) + " )" + " + "
            "( " + repr(self.coeffs[5]) + " )" + " * II" + " + "
            "( " + repr(self.coeffs[6]) + " )" + " * JJ" + " + "
            "( " + repr(self.coeffs[7]) + " )" + " * KK)"
        )
        return result

    def __str__(self):
        """String function.

        Parameters
        ----------
        self: self

        Returns
        -------
        result : string
        Human readable string representation of biquaternion.
        """
        result = (
            "(" + "( " + repr(self.coeffs[0]) + " )" + " + "
            "( " + repr(self.coeffs[1]) + " )" + " * i" + " + "
            "( " + repr(self.coeffs[2]) + " )" + " * j" + " + "
            "( " + repr(self.coeffs[3]) + " )" + " * k"
            ") + eps * ("
            "( " + repr(self.coeffs[4]) + " )" + " + "
            "( " + repr(self.coeffs[5]) + " )" + " * i" + " + "
            "( " + repr(self.coeffs[6]) + " )" + " * j" + " + "
            "( " + repr(self.coeffs[7]) + " )" + " * k)"
        )
        return result

    def __pow__(self, other):
        if isinstance(other, int):
            if other >= 0:
                pw = 1
                for i in range(other):
                    pw = pw * self
            else:
                pw = 1
                for i in range(-other):
                    pw = pw / self
            return pw
        else:
            raise TypeError(
                "unsupported operand type(s) for ** or pow(): "
                + str(type(self))
                + " and "
                + str(type(other))
            )

    def conjugate(self):
        """(Bi)-Quaternion conjugate of the quaternion.

        Conjugation of a quaternion inverts the sign of the non-scalar part of
        a (bi-)quaternion.
        This happens in the same fashion as for complex numbers.
        """
        return BiQuaternion(
            *[
                self.coeffs[0],
                -self.coeffs[1],
                -self.coeffs[2],
                -self.coeffs[3],
                self.coeffs[4],
                -self.coeffs[5],
                -self.coeffs[6],
                -self.coeffs[7],
            ]
        )

    def eps_conjugate(self):
        """Epsilon conjugation of the biquaternion

        Epsilon conjugation inverts the sign of the dual part of a quaternion
        """
        return BiQuaternion(
            *[
                self.coeffs[0],
                self.coeffs[1],
                self.coeffs[2],
                self.coeffs[3],
                -self.coeffs[4],
                -self.coeffs[5],
                -self.coeffs[6],
                -self.coeffs[7],
            ]
        )

    def quadrance(self):
        """Quadrance of a quaternion

        The quadrance of a biquaternion is its norm. It is defined as
        the product of a biquaternion and its conjugate.
        """
        return self * (self.conjugate())

    def __invert__(self):
        """(Bi)-Quaternion conjugate of the quaternion.

        Conjugation of a quaternion changes the sign of the non-scalar part of
        a (bi-)quaternion.
        This happens in the same fashion as for complex numbers.
        """
        return self.conjugate()

    def inv(self):
        """Calculates the inverse of the biquaternion"""
        quad = self.quadrance()
        primal = quad.coeffs[0]
        dual = quad.coeffs[5]
        s = primal * primal - _BQ_E * dual * dual
        if s == 0:
            raise ValueError("Object is not invertible")
            return
        return (quad.eps_conjugate() * (1 / s)) * (~self)

    def __truediv__(self, other):
        """BiQuaternion division"""
        if isinstance(other, BiQuaternion):
            return self * other.inv()
        return self * (1 / other)

    def __rtruediv__(self, other):
        """BiQuaternion division"""
        return other * self.inv()

    def primal(self):
        """Calculates the primal part of the dual quaternion.

        Returns
        -------
        BiQuaternion

        Notes
        -----
        The primal part of a dual quaternion is defined as the part not containing a
        factor epsilon.

        """
        return BiQuaternion(*self.coeffs[0:4])

    def dual(self):
        """Calculates the dual part of the dual quaternion

        Returns
        -------
        BiQuaternion

        Notes
        -----
        The dual part of a dual quaternion is defined as the part only containing
        factors epsilon.

        """
        return BiQuaternion(*self.coeffs[4:])

    def scalar_part(self):
        """Calculates the scalar part of the dual quaternion

        Returns
        -------
        BiQuaternion

        Notes
        -----
        The scalar part of a dual quaternion is defined as the part only not containing
        any of the numbers i, j, or k

        """
        return BiQuaternion([self.coeffs[0], 0, 0, 0, self.coeffs[5], 0, 0, 0])

    def vector_part(self):
        """Calculates the vector part of the dual quaternion

        Returns
        -------
        BiQuaternion

        Notes
        -----
        The dual part of a dual quaternion is defined as the part only containing
        the numbers i, j, k.

        """
        return BiQuaternion(*([0] + self.coeffs[1:4] + [0] + self.coeffs[5:]))


II = BiQuaternion(0, 1, 0, 0, 0, 0, 0, 0)
JJ = BiQuaternion(0, 0, 1, 0, 0, 0, 0, 0)
KK = BiQuaternion(0, 0, 0, 1, 0, 0, 0, 0)
EE = BiQuaternion(0, 0, 0, 0, 1, 0, 0, 0)
