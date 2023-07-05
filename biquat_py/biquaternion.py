"""
This module implements the biquaternion class, which takes care of
 the basic biquaternion arithmetic.
"""

import numpy as np

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
    II = BiQuaternion([0, 1, 0, 0, 0, 0, 0, 0])
    JJ = BiQuaternion([0, 0, 1, 0, 0, 0, 0, 0])
    KK = BiQuaternion([0, 0, 0, 1, 0, 0, 0, 0])
    EE = BiQuaternion([0, 0, 0, 0, 1, 0, 0, 0])


class BiQuaternion:
    """Biquaternions as $a + II b + JJ c + KK d + EE (w + II x + JJ y + KK z)$"""

    def __init__(self, gen):
        if _BACKEND == "numpy":
            self.coeff = np.zeros(8)
        else:
            self.coeff = [0, 0, 0, 0, 0, 0, 0, 0]

        if isinstance(gen, BiQuaternion):
            cof = gen.coeff
        elif isinstance(gen, (list, tuple, np.ndarray)):
            if len(gen) >= 9:
                raise ValueError("Maximum array length is 8")
            cof = gen
        else:
            cof = [gen]

        for i, val in enumerate(cof):
            self.coeff[i] = val

    def __mul__(self, other):
        if isinstance(other, BiQuaternion):
            out = [
                -_BQ_E * _BQ_I * _BQ_J * other.coeff[7] * self.coeff[7]
                + _BQ_E * _BQ_I * other.coeff[5] * self.coeff[5]
                + _BQ_E * _BQ_J * other.coeff[6] * self.coeff[6]
                - _BQ_I * _BQ_J * other.coeff[3] * self.coeff[3]
                + _BQ_E * other.coeff[4] * self.coeff[4]
                + _BQ_I * other.coeff[1] * self.coeff[1]
                + _BQ_J * other.coeff[2] * self.coeff[2]
                + other.coeff[0] * self.coeff[0],
                _BQ_E * _BQ_J * other.coeff[6] * self.coeff[7]
                - _BQ_E * _BQ_J * other.coeff[7] * self.coeff[6]
                + _BQ_E * other.coeff[4] * self.coeff[5]
                + _BQ_E * other.coeff[5] * self.coeff[4]
                + _BQ_J * other.coeff[2] * self.coeff[3]
                - _BQ_J * other.coeff[3] * self.coeff[2]
                + other.coeff[0] * self.coeff[1]
                + other.coeff[1] * self.coeff[0],
                -_BQ_E * _BQ_I * other.coeff[5] * self.coeff[7]
                + _BQ_E * _BQ_I * other.coeff[7] * self.coeff[5]
                + _BQ_E * other.coeff[4] * self.coeff[6]
                + _BQ_E * other.coeff[6] * self.coeff[4]
                - _BQ_I * other.coeff[1] * self.coeff[3]
                + _BQ_I * other.coeff[3] * self.coeff[1]
                + other.coeff[0] * self.coeff[2]
                + other.coeff[2] * self.coeff[0],
                _BQ_E * other.coeff[4] * self.coeff[7]
                - _BQ_E * other.coeff[5] * self.coeff[6]
                + _BQ_E * other.coeff[6] * self.coeff[5]
                + _BQ_E * other.coeff[7] * self.coeff[4]
                + other.coeff[0] * self.coeff[3]
                - other.coeff[1] * self.coeff[2]
                + other.coeff[2] * self.coeff[1]
                + other.coeff[3] * self.coeff[0],
                -_BQ_I * _BQ_J * other.coeff[3] * self.coeff[7]
                - _BQ_I * _BQ_J * other.coeff[7] * self.coeff[3]
                + _BQ_I * other.coeff[1] * self.coeff[5]
                + _BQ_I * other.coeff[5] * self.coeff[1]
                + _BQ_J * other.coeff[2] * self.coeff[6]
                + _BQ_J * other.coeff[6] * self.coeff[2]
                + other.coeff[0] * self.coeff[4]
                + other.coeff[4] * self.coeff[0],
                _BQ_J * other.coeff[2] * self.coeff[7]
                - _BQ_J * other.coeff[3] * self.coeff[6]
                + _BQ_J * other.coeff[6] * self.coeff[3]
                - _BQ_J * other.coeff[7] * self.coeff[2]
                + other.coeff[0] * self.coeff[5]
                + other.coeff[1] * self.coeff[4]
                + other.coeff[4] * self.coeff[1]
                + other.coeff[5] * self.coeff[0],
                -_BQ_I * other.coeff[1] * self.coeff[7]
                + _BQ_I * other.coeff[3] * self.coeff[5]
                - _BQ_I * other.coeff[5] * self.coeff[3]
                + _BQ_I * other.coeff[7] * self.coeff[1]
                + other.coeff[0] * self.coeff[6]
                + other.coeff[2] * self.coeff[4]
                + other.coeff[4] * self.coeff[2]
                + other.coeff[6] * self.coeff[0],
                other.coeff[0] * self.coeff[7]
                - other.coeff[1] * self.coeff[6]
                + other.coeff[2] * self.coeff[5]
                + other.coeff[3] * self.coeff[4]
                + other.coeff[4] * self.coeff[3]
                - other.coeff[5] * self.coeff[2]
                + other.coeff[6] * self.coeff[1]
                + other.coeff[7] * self.coeff[0],
            ]
            return BiQuaternion(out)

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
        return BiQuaternion([-self.coeff[i] for i in range(8)])

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
            out = [self.coeff[i] + other.coeff[i] for i in range(8)]
            return BiQuaternion(out)

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
        other_coeff = BiQuaternion(other).coeff
        for i, val in enumerate(self.coeff):
            if val != other_coeff[i]:
                return False
        return True

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
            "(" + "( " + repr(self.coeff[0]) + " )" + " + "
            "( " + repr(self.coeff[1]) + " )" + " * II" + " + "
            "( " + repr(self.coeff[2]) + " )" + " * JJ" + " + "
            "( " + repr(self.coeff[3]) + " )" + " * KK"
            ") + EE * ("
            "( " + repr(self.coeff[4]) + " )" + " + "
            "( " + repr(self.coeff[5]) + " )" + " * II" + " + "
            "( " + repr(self.coeff[6]) + " )" + " * JJ" + " + "
            "( " + repr(self.coeff[7]) + " )" + " * KK)"
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
            "(" + "( " + repr(self.coeff[0]) + " )" + " + "
            "( " + repr(self.coeff[1]) + " )" + " * i" + " + "
            "( " + repr(self.coeff[2]) + " )" + " * j" + " + "
            "( " + repr(self.coeff[3]) + " )" + " * k"
            ") + eps * ("
            "( " + repr(self.coeff[4]) + " )" + " + "
            "( " + repr(self.coeff[5]) + " )" + " * i" + " + "
            "( " + repr(self.coeff[6]) + " )" + " * j" + " + "
            "( " + repr(self.coeff[7]) + " )" + " * k)"
        )
        return result

    def conjugate(self):
        """(Bi)-Quaternion conjugate of the quaternion.

        Conjugation of a quaternion inverts the sign of the non-scalar part of
        a (bi-)quaternion.
        This happens in the same fashion as for complex numbers.
        """
        return BiQuaternion(
            [
                self.coeff[0],
                -self.coeff[1],
                -self.coeff[2],
                -self.coeff[3],
                self.coeff[4],
                -self.coeff[5],
                -self.coeff[6],
                -self.coeff[7],
            ]
        )

    def eps_conjugate(self):
        """Epsilon conjugation of the biquaternion

        Epsilon conjugation inverts the sign of the dual part of a quaternion
        """
        return BiQuaternion(
            [
                self.coeff[0],
                self.coeff[1],
                self.coeff[2],
                self.coeff[3],
                -self.coeff[4],
                -self.coeff[5],
                -self.coeff[6],
                -self.coeff[7],
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
        primal = quad.coeff[0]
        dual = quad.coeff[5]
        s = primal * primal - _BQ_E * dual * dual
        if s == 0:
            raise ValueError("Object is not invertible")
            return
        return (quad.eps_conjugate() * (1 / s)) * (~self)


II = BiQuaternion([0, 1, 0, 0, 0, 0, 0, 0])
JJ = BiQuaternion([0, 0, 1, 0, 0, 0, 0, 0])
KK = BiQuaternion([0, 0, 0, 1, 0, 0, 0, 0])
EE = BiQuaternion([0, 0, 0, 0, 1, 0, 0, 0])
