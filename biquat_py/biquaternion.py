global _bq_i, _bq_j, _bq_e

_bq_i = -1
_bq_j = -1
_bq_e = 0

def define_algebra(a = -1, b = -1, c = 0):
    _bq_i = a
    _bq_j = b
    _bq_e = c

class biquaternion(a = 0, ai = 0, aj = 0, ak = 0, ae = 0, aei = 0, aej = 0, aek = 0):

    """Docstring for biquaternion. """

    coeff = [0,0,0,0,0,0,0,0]

    def __init__(self):
        self.coeff[0] = a
        self.coeff[1] = ai
        self.coeff[2] = aj
        self.coeff[3] = ak
        self.coeff[4] = ae
        self.coeff[5] = aei
        self.coeff[6] = aej
        self.coeff[7] = aek

    def mul(self,other):
        self.coeff = [-_bq_e*_bq_i*_bq_j*other.coeff[7]*self.coeff[7] +
        _bq_e*_bq_i*other.coeff[5]*self.coeff[5] +
        _bq_e*_bq_j*other.coeff[6]*self.coeff[6] -
        _bq_i*_bq_j*other.coeff[3]*self.coeff[3] +
        _bq_e*other.coeff[4]*self.coeff[4] + _bq_i*other.coeff[1]*self.coeff[1]
        + _bq_j*other.coeff[2]*self.coeff[2] + other.coeff[0]*self.coeff[0],
        _bq_e*_bq_j*other.coeff[6]*self.coeff[7] -
        _bq_e*_bq_j*other.coeff[7]*self.coeff[6] +
        _bq_e*other.coeff[4]*self.coeff[5] + _bq_e*other.coeff[5]*self.coeff[4]
        + _bq_j*other.coeff[2]*self.coeff[3] -
        _bq_j*other.coeff[3]*self.coeff[2] + other.coeff[0]*self.coeff[1] +
        other.coeff[1]*self.coeff[0], -_bq_e*_bq_i*other.coeff[5]*self.coeff[7]
        + _bq_e*_bq_i*other.coeff[7]*self.coeff[5] +
        _bq_e*other.coeff[4]*self.coeff[6] + _bq_e*other.coeff[6]*self.coeff[4]
        - _bq_i*other.coeff[1]*self.coeff[3] +
        _bq_i*other.coeff[3]*self.coeff[1] + other.coeff[0]*self.coeff[2] +
        other.coeff[2]*self.coeff[0], _bq_e*other.coeff[4]*self.coeff[7] -
        _bq_e*other.coeff[5]*self.coeff[6] + _bq_e*other.coeff[6]*self.coeff[5]
        + _bq_e*other.coeff[7]*self.coeff[4] + other.coeff[0]*self.coeff[3] -
        other.coeff[1]*self.coeff[2] + other.coeff[2]*self.coeff[1] +
        other.coeff[3]*self.coeff[0], -_bq_i*_bq_j*other.coeff[3]*self.coeff[7]
        - _bq_i*_bq_j*other.coeff[7]*self.coeff[3] +
        _bq_i*other.coeff[1]*self.coeff[5] + _bq_i*other.coeff[5]*self.coeff[1]
        + _bq_j*other.coeff[2]*self.coeff[6] +
        _bq_j*other.coeff[6]*self.coeff[2] + other.coeff[0]*self.coeff[4] +
        other.coeff[4]*self.coeff[0], _bq_j*other.coeff[2]*self.coeff[7] -
        _bq_j*other.coeff[3]*self.coeff[6] + _bq_j*other.coeff[6]*self.coeff[3]
        - _bq_j*other.coeff[7]*self.coeff[2] + other.coeff[0]*self.coeff[5] +
        other.coeff[1]*self.coeff[4] + other.coeff[4]*self.coeff[1] +
        other.coeff[5]*self.coeff[0], -_bq_i*other.coeff[1]*self.coeff[7] +
        _bq_i*other.coeff[3]*self.coeff[5] - _bq_i*other.coeff[5]*self.coeff[3]
        + _bq_i*other.coeff[7]*self.coeff[1] + other.coeff[0]*self.coeff[6] +
        other.coeff[2]*self.coeff[4] + other.coeff[4]*self.coeff[2] +
        other.coeff[6]*self.coeff[0], other.coeff[0]*self.coeff[7] -
        other.coeff[1]*self.coeff[6] + other.coeff[2]*self.coeff[5] +
        other.coeff[3]*self.coeff[4] + other.coeff[4]*self.coeff[3] -
        other.coeff[5]*self.coeff[2] + other.coeff[6]*self.coeff[1] +
        other.coeff[7]*self.coeff[0]]
