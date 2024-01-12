"""
Module containing extra generation tools and function for specialized BiQuaternions.
"""
import numpy as np
from .biquaternion import BiQuaternion, EE
from .lines import quat_to_pluecker, pluecker_to_quat, act_on_line


def point_to_quat(coord):
    """Generate BiQuaternion reperesentation of a point at specified coordinates."""
    return BiQuaternion([1, 0, 0, 0, 0, *coord])


def quat_to_point(quat):
    """Generate coordinates of a point defined by a BiQuaternion."""
    # Normalize point.
    quat *= 1 / (quat.scal)
    if quat.coeffs[1:5] == [0, 0, 0, 0]:
        return [quat.ei, quat.ej, quat.ek]
    else:
        raise ValueError("Object not a valid description of a point")


def hom_point_to_quat(coord):
    """Generate BiQuaternion reperesentation of a point at specified
    homogeneous coordinates."""
    return BiQuaternion([coord[0], 0, 0, 0, 0, *(coord[1:])])


def quat_to_hom_point(quat):
    """Generate coordinates of a homogeneous point defined by a BiQuaternion."""
    if quat.coeffs[1:5] == [0, 0, 0, 0]:
        return [quat.scal, quat.ei, quat.ej, quat.ek]
    else:
        raise ValueError("Object not a valid description of a point")


def act_on_point(quaternion, x):
    """Let a BiQuaternion act on a point."""
    return quaternion.eps_conjugate() * x * quaternion.conjugate()


def smart_act(quat, obj):
    """General purpose function for letting BiQuaternions act on object that detects
    object type.
    """
    if isinstance(obj, BiQuaternion):
        # Detect if it is a line
        if obj.scal.expand() == 0 and obj.primal() != 0:
            return act_on_line(quat, obj)
        return act_on_point(quat, obj)
    elif isinstance(obj, (list, tuple, np.ndarray)):
        if len(obj) == 3:
            return quat_to_point(act_on_point(quat, point_to_quat(obj)))
        elif len(obj) == 4:
            return quat_to_hom_point(act_on_point(quat, hom_point_to_quat(obj)))
        elif len(obj) == 6:
            return quat_to_pluecker(act_on_line(quat, pluecker_to_quat(obj)))
        else:
            raise ValueError("Object type not detected.")
    else:
        raise TypeError(
            "obj must be of type BiQuaternion, or (list, tuple, np.ndarray)."
        )


def inner(first, second):
    """Inner product of first and second, derived from the quaternion algebra."""
    return (first * second.conjugate() + second * first.conjugate()) / 2


def outer(first, second):
    """Outer product of first and second, derived from the quaternion algebra."""
    return (first * second - second * first) / 2


def fiber_project(quat):
    """Project biquaternion onto Study's quadric.

    Parameters
    ----------
    quat : BiQuaternion
        Quaternion which is to be projected

    Returns
    -------
    BiQuaternion
        Biquaternion lying on Study's quadric, corresponding to the same transformation.
    """

    primal = quat.primal()
    dual = quat.dual()
    return (
        (
            2 * primal.quadrance()
            - EE * (primal * dual.conjugate() - dual * primal.conjugate())
        )
        * primal
    ) / 2
