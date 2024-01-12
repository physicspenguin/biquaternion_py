import biquaternion_py as bq
import numpy.testing as nt


def test_pluecker_to_quat():
    coord = [1, 0, 0, 0, 1, 1]
    assert bq.BiQuaternion([0, 1, 0, 0, 0, 0, -1, -1]) == bq.pluecker_to_quat(coord)


def test_quat_to_pluecker():
    quat = bq.BiQuaternion([0, 1, 0, 0, 0, 0, 1, 1])
    fail_quat = bq.BiQuaternion([0, 1, 2, 3, 5, 4, 5, 6])

    assert bq.quat_to_pluecker(quat) == [1, 0, 0, 0, -1, -1]

    with nt.assert_raises(ValueError):
        bq.quat_to_pluecker(fail_quat)
