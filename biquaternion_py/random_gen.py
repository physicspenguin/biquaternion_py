from .biquaternion import BiQuaternion
import numpy.random as rand
from sympy import Rational


def rand_rational(maximum=10):
    """Generate random rational number up to n

    Parameters
    ----------
    maximum : int (optional)
        Maximum of random number generated
         (Default value = 10)

    Returns
    -------
    (float)
        Random rational number in Range maximum

    """
    return (-1) ** (rand.randint(2)) * (
        Rational(rand.randint(maximum), (rand.randint(maximum) + 1))
    )


def rand_quat(maximum=10):
    """Random quaternion with rational coefficients."""
    return BiQuaternion([rand_rational(maximum) for i in range(4)])


def rand_bq(maximum=10):
    """Random biquaternion with rational coefficients."""
    return BiQuaternion([rand_rational(maximum) for i in range(8)])


def rand_line(maximum=10):
    p = [rand_rational(maximum) for i in range(4)]
    q = [rand_rational(maximum) for i in range(4)]
    return BiQuaternion(
        [
            0,
            (p[0] * q[1] - p[1] * q[0]),
            (p[0] * q[2] - p[2] * q[0]),
            (p[0] * q[3] - p[3] * q[0]),
            0,
            (p[2] * q[3] - p[3] * q[2]),
            (p[3] * q[1] - p[1] * q[3]),
            (p[1] * q[2] - p[2] * q[1]),
        ]
    )
