global _bq_i, _bq_j, _bq_e
import numpy as np

_bq_i = -1
_bq_j = -1
_bq_e = 0


def define_algebra(a=-1, b=-1, c=0):
    _bq_i = a
    _bq_j = b
    _bq_e = c
    II = biquaternion([0, 1, 0, 0, 0, 0, 0, 0])
    JJ = biquaternion([0, 0, 1, 0, 0, 0, 0, 0])
    KK = biquaternion([0, 0, 0, 1, 0, 0, 0, 0])
    EE = biquaternion([0, 0, 0, 0, 1, 0, 0, 0])


class biquaternion:

    """Biquaternions as $a + II b + JJ c + KK d + EE (w + II x + JJ y + KK z)$"""

    coeff = [0, 0, 0, 0, 0, 0, 0, 0]

    def __init__(self, gen):
        if isinstance(gen, biquaternion):
            cof = gen.coeff
        elif isinstance(gen, (list, tuple, np.ndarray)):
            if len(gen) >= 9:
                raise ValueError("Maximum array length is 8")
            cof = gen
        else:
            cof = [gen]

        self.coeff = [0, 0, 0, 0, 0, 0, 0, 0]

        for i in range(len(cof)):
            self.coeff[i] = cof[i]

    def __mul__(self, other):
        if isinstance(other, biquaternion):
            out = [
                -_bq_e * _bq_i * _bq_j * other.coeff[7] * self.coeff[7]
                + _bq_e * _bq_i * other.coeff[5] * self.coeff[5]
                + _bq_e * _bq_j * other.coeff[6] * self.coeff[6]
                - _bq_i * _bq_j * other.coeff[3] * self.coeff[3]
                + _bq_e * other.coeff[4] * self.coeff[4]
                + _bq_i * other.coeff[1] * self.coeff[1]
                + _bq_j * other.coeff[2] * self.coeff[2]
                + other.coeff[0] * self.coeff[0],
                _bq_e * _bq_j * other.coeff[6] * self.coeff[7]
                - _bq_e * _bq_j * other.coeff[7] * self.coeff[6]
                + _bq_e * other.coeff[4] * self.coeff[5]
                + _bq_e * other.coeff[5] * self.coeff[4]
                + _bq_j * other.coeff[2] * self.coeff[3]
                - _bq_j * other.coeff[3] * self.coeff[2]
                + other.coeff[0] * self.coeff[1]
                + other.coeff[1] * self.coeff[0],
                -_bq_e * _bq_i * other.coeff[5] * self.coeff[7]
                + _bq_e * _bq_i * other.coeff[7] * self.coeff[5]
                + _bq_e * other.coeff[4] * self.coeff[6]
                + _bq_e * other.coeff[6] * self.coeff[4]
                - _bq_i * other.coeff[1] * self.coeff[3]
                + _bq_i * other.coeff[3] * self.coeff[1]
                + other.coeff[0] * self.coeff[2]
                + other.coeff[2] * self.coeff[0],
                _bq_e * other.coeff[4] * self.coeff[7]
                - _bq_e * other.coeff[5] * self.coeff[6]
                + _bq_e * other.coeff[6] * self.coeff[5]
                + _bq_e * other.coeff[7] * self.coeff[4]
                + other.coeff[0] * self.coeff[3]
                - other.coeff[1] * self.coeff[2]
                + other.coeff[2] * self.coeff[1]
                + other.coeff[3] * self.coeff[0],
                -_bq_i * _bq_j * other.coeff[3] * self.coeff[7]
                - _bq_i * _bq_j * other.coeff[7] * self.coeff[3]
                + _bq_i * other.coeff[1] * self.coeff[5]
                + _bq_i * other.coeff[5] * self.coeff[1]
                + _bq_j * other.coeff[2] * self.coeff[6]
                + _bq_j * other.coeff[6] * self.coeff[2]
                + other.coeff[0] * self.coeff[4]
                + other.coeff[4] * self.coeff[0],
                _bq_j * other.coeff[2] * self.coeff[7]
                - _bq_j * other.coeff[3] * self.coeff[6]
                + _bq_j * other.coeff[6] * self.coeff[3]
                - _bq_j * other.coeff[7] * self.coeff[2]
                + other.coeff[0] * self.coeff[5]
                + other.coeff[1] * self.coeff[4]
                + other.coeff[4] * self.coeff[1]
                + other.coeff[5] * self.coeff[0],
                -_bq_i * other.coeff[1] * self.coeff[7]
                + _bq_i * other.coeff[3] * self.coeff[5]
                - _bq_i * other.coeff[5] * self.coeff[3]
                + _bq_i * other.coeff[7] * self.coeff[1]
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
            return biquaternion(out)
        else:
            return self * biquaternion(other)

    def __pos__(self):
        return biquaternion(self)

    def __neg__(self):
        return biquaternion([-self.coeff[i] for i in range(8)])

    def __add__(self, other):
        if isinstance(other, biquaternion):
            out = [self.coeff[i] + other.coeff[i] for i in range(8)]
            return biquaternion(out)
        else:
            return self + biquaternion(other)

    def __sub__(self, other):
        return self + (-other)

    def __rsub__(self, other):
        return other + (-self)

    __radd__ = __add__
    __rsub__ = __sub__
    __rmul__ = __mul__

    def __repr__(self):
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


II = biquaternion([0, 1, 0, 0, 0, 0, 0, 0])
JJ = biquaternion([0, 0, 1, 0, 0, 0, 0, 0])
KK = biquaternion([0, 0, 0, 1, 0, 0, 0, 0])
EE = biquaternion([0, 0, 0, 0, 1, 0, 0, 0])
