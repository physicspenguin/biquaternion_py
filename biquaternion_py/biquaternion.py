"""Calculate with BiQuaternions.

This module implements BiQuaternions as a class for general calculations.

Classes:

    BiQuaternion

Misc variables:

    II
    JJ
    KK
    EE
"""

import numpy as np
from sympy.core.expr import Expr
from sympy import sympify, expand
from .polynomials import Poly

_BQ_I = -1
_BQ_J = -1
_BQ_E = 0


def define_algebra(i_square=-1, j_square=-1, e_square=0):
    r"""Define the algebra.

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
    """Sanitizes the input of the __new__ method of BiQuaternion."""
    coeffs = [0, 0, 0, 0, 0, 0, 0, 0]
    if len(args) == 1:
        gen = args[0]
        cof = [gen]
        if isinstance(gen, BiQuaternion):
            cof = gen.coeffs
        elif isinstance(gen, (list, tuple, np.ndarray)):
            if len(gen) >= 9:
                raise ValueError("Maximum array length is 8")
            cof = gen
    else:
        cof = args

    for i, val in enumerate(cof):
        coeffs[i] = val

    return coeffs


class BiQuaternion(Expr):
    """
    Class implementing Bi-Quaternions.

    Bi-Quaternions are represented as
    $a + II b + JJ c + KK d + EE (w + II x + JJ y + KK z)$.

    Attributes
    ----------
    coeffs : list
        coefficients of the quaternion as a list in the canonical order.
    scal : sympy.Expr, numeric
        scalar value of the Bi-Quaternion
    i : sympy.Expr, numeric
        II value of the Bi-Quaternion
    j : sympy.Expr, numeric
        JJ value of the Bi-Quaternion
    k : sympy.Expr, numeric
        KK value of the Bi-Quaternion
    eps : sympy.Expr, numeric
        dual scalar value of the Bi-Quaternion
    ei : sympy.Expr, numeric
        dual II value of the Bi-Quaternion
    ej : sympy.Expr, numeric
        dual JJ value of the Bi-Quaternion
    ek : sympy.Expr, numeric
        dual KK value of the Bi-Quaternion
    conjugate : BiQuaternion
        Conjugate of this instance of BiQuaternion.
    eps_conjugate : BiQuaternion
        Epsilon conjugation of the biquaternion.
    quadrance : BiQuaternion
        Quadrance of a quaternion.
    inv : BiQuaternion
        Inverse of the biquaternion.
    primal : BiQuaternion
        Primal part of the dual quaternion.
    dual : BiQuaternion
        Dual part of the dual quaternion.
    scalar_part : BiQuaternion
        Scalar part of the dual quaternion.
    vector_part : BiQuaternion
        Vector part of the dual quaternion.

    Methods
    -------
    __new__(cls, *args):
        Create new instance of the BiQuaternion class.
    __mul__(other):
        Multiply BiQuaternion with other.
    __pos__(self):
        Positive of itsself.
    __neg__(self):
        Negative of itsself.
    __add__(other):
        Add BiQuaternion to other.
    __sub__(other):
        Subtract other from BiQuaternion.
    __rsub__(other):
        Subtract BiQuaternion from other.
    __radd__ = __add__
    __rmul__ = __mul__
    __eq__(other):
        Test equality of two biquaternions.
    __hash__ = super.__hash__
    __repr__():
        Convert BiQuaternion to a readable format in shell.
    __str__():
        Converst BiQuaternion to string.
    __pow__(other):
        Power function of BiQuaternion.
    __invert__():
        (Bi)-Quaternion conjugate of the quaternion.
    __truediv__(other):
        Division of BiQuaternion by other.
    __rtruediv__(other):
        Divide other by BiQuaternion.
    coeff(var, power):
        Rewriting of Expr.coeff to work for BiQuaternions
    """

    is_commutative = False
    _op_priority = 11.1

    def __new__(cls, *args):
        """Create new instance of BiQuaternion."""
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

    @property
    def scal(self):
        """Value of the scalar part of the instance of BiQuaternion."""
        return self._scal

    @scal.setter
    def scal(self, val):
        self._scal = val

    @property
    def i(self):
        """Value of the II part of the instance of BiQuaternion."""
        return self._i

    @i.setter
    def i(self, val):
        self._i = val

    @property
    def j(self):
        """Value of the JJ part of the instance of BiQuaternion."""
        return self._j

    @j.setter
    def j(self, val):
        self._j = val

    @property
    def k(self):
        """Value of the KK part of the instance of BiQuaternion."""
        return self._k

    @k.setter
    def k(self, val):
        self._k = val

    @property
    def eps(self):
        """Value of the eps part of the instance of BiQuaternion."""
        return self._eps

    @eps.setter
    def eps(self, val):
        self._eps = val

    @property
    def ei(self):
        """Value of the eps*II part of the instance of BiQuaternion."""
        return self._ei

    @ei.setter
    def ei(self, val):
        self._ei = val

    @property
    def ej(self):
        """Value of the eps*JJ part of the instance of BiQuaternion."""
        return self._ej

    @ej.setter
    def ej(self, val):
        self._ej = val

    @property
    def ek(self):
        """Value of the eps*KK part of the instance of BiQuaternion."""
        return self._ek

    @ek.setter
    def ek(self, val):
        self._ek = val

    @property
    def coeffs(self):
        """Coefficients describing an instance of BiQuaternion."""
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

    @coeffs.setter
    def coeffs(self, val):
        self.scal = val[0]
        self.i = val[1]
        self.j = val[2]
        self.k = val[3]
        self.eps = val[4]
        self.ei = val[5]
        self.ej = val[6]
        self.ek = val[7]

    def __mul__(self, other):
        """Multiply BiQuaternion with other."""
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
        elif isinstance(other, Poly):
            return other.__rmul__(self)

        return self * BiQuaternion(other)

    def __pos__(self):
        """Positive of itsself.

        Returns
        -------
        BiQuaternion
            self
        """
        return BiQuaternion(self)

    def __neg__(self):
        """Negative of itsself.

        Returns
        -------
        BiQuaternion
            -self
        """
        return BiQuaternion(*[-self.coeffs[i] for i in range(8)])

    def __add__(self, other):
        """Add BiQuaternion to other.

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

        elif isinstance(other, Poly):
            return other.__radd__(self)

        return self + BiQuaternion(other)

    def __sub__(self, other):
        """Subtract other from BiQuaternion.

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
        """Subtract BiQuaternion from other.

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
        """Test equality of two biquaternions."""
        othercoeff = BiQuaternion(other).coeffs
        for i, val in enumerate(self.coeffs):
            if val != othercoeff[i]:
                return False
        return True

    __hash__ = super.__hash__

    def __repr__(self):
        """Convert BiQuaternion to a readable format in shell.

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
        """Converst BiQuaternion to string.

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
        """Power function of BiQuaternion."""
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
        """Conjugate of this instance of BiQuaternion.

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
        """Epsilon conjugation of the biquaternion.

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
        """Quadrance of a quaternion.

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
        """Inverse of the biquaternion."""
        quad = self.quadrance()
        primal = quad.coeffs[0]
        dual = quad.coeffs[5]
        s = primal * primal - _BQ_E * dual * dual
        if s == 0:
            raise ValueError("Object is not invertible")
            return
        return (quad.eps_conjugate() * (1 / s)) * (~self)

    def __truediv__(self, other):
        """Division of BiQuaternion by other."""
        if isinstance(other, BiQuaternion):
            return self * other.inv()
        return self * (1 / other)

    def __rtruediv__(self, other):
        """Divide other by BiQuaternion."""
        return other * self.inv()

    def primal(self):
        """Primal part of the dual quaternion.

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
        """Dual part of the dual quaternion.

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
        """Scalar part of the dual quaternion.

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
        """Vector part of the dual quaternion.

        Returns
        -------
        BiQuaternion

        Notes
        -----
        The dual part of a dual quaternion is defined as the part only containing
        the numbers i, j, k.

        """
        return BiQuaternion(*([0] + self.coeffs[1:4] + [0] + self.coeffs[5:]))

    def coeff(self, var, power=1, right=False, _first=True):
        """Rewriting of Expr.coeff to work for BiQuaternions."""
        # TODO: Get every coeff after each other. Run expr.coeff on it and
        # construct new quaternion describing the coeffs.
        coeffed = [0, 0, 0, 0, 0, 0, 0, 0]
        for i, val in enumerate(self.coeffs):
            coeffed[i] = expand(val).coeff(var, power, right, _first)
        return BiQuaternion(coeffed)

    def apply_elementwise(self, func, *args):
        """Apply a function with specified arguments elementwise.

        Parameters
        ----------
        func : function
            Function to be applied elementwise
        *args : unknown
            Arguments to be passed to the function

        Returns
        -------
        BiQuaternion
            Biquaternion with function applied to each coefficient individually.

        """
        coeffed = [0, 0, 0, 0, 0, 0, 0, 0]
        for i, val in enumerate(self.coeffs):
            coeffed[i] = func(val, *args)
        return BiQuaternion(coeffed)


II = BiQuaternion(0, 1, 0, 0, 0, 0, 0, 0)
JJ = BiQuaternion(0, 0, 1, 0, 0, 0, 0, 0)
KK = BiQuaternion(0, 0, 0, 1, 0, 0, 0, 0)
EE = BiQuaternion(0, 0, 0, 0, 1, 0, 0, 0)
