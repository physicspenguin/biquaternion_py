import biquaternion_py as bq
import numpy.testing as nt


def test_point_to_quat():
    assert bq.point_to_quat([1, 2, 3]) == bq.BiQuaternion([1, 0, 0, 0, 0, 1, 2, 3])


def test_quat_to_point():
    a = [1, 2, 3]
    assert bq.quat_to_point(2 * bq.point_to_quat(a)) == a
    assert bq.quat_to_point(bq.point_to_quat(a)) == a
    with nt.assert_raises(ValueError):
        bq.quat_to_point(bq.BiQuaternion([1, 2, 3, 4, 5, 6, 7, 8]))


def test_hom_point_to_quat():
    a = [5, 1, 2, 3]
    assert bq.hom_point_to_quat(a) == bq.BiQuaternion([5, 0, 0, 0, 0, 1, 2, 3])


def test_quat_to_hom_point():
    a = [5, 1, 2, 3]
    assert bq.quat_to_hom_point(2 * bq.hom_point_to_quat(a)) != a
    assert bq.quat_to_hom_point(bq.hom_point_to_quat(a)) == a
    with nt.assert_raises(ValueError):
        bq.quat_to_hom_point(bq.BiQuaternion([1, 2, 3, 4, 5, 6, 7, 8]))
