import numpy as np
from .biquaternion import BiQuaternion, EE


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


def pluecker_to_quat(coord):
    """Generate BiQuaternion reperesentation of a line with specified
    Pluecker coordinates."""
    return BiQuaternion([0, *coord[0:3], 0, *[-coord[i + 3] for i in range(3)]])


def quat_to_pluecker(quat):
    """Generate pluecker coordinates from a given BiQuaternion."""
    if quat.scal == 0 and quat.eps == 0:
        return [*quat.coeffs[1, 2, 3], *(-quat[5, 6, 7])]
    else:
        raise ValueError("Object not a valid description of a line")


def line_to_pluecker(direction, point):
    """Generate Pluecker coordinates for a line with given direction through a point."""
    return [*direction, -np.cross(direction, point)]


def act_on_point(quaternion, x):
    """Let a BiQuaternion act on a point."""
    return quaternion.eps_conjugate() * x * quaternion.conjugate()


def act_on_line(quaternion, lin):
    """Let a BiQuaternion act on a line."""
    return quaternion * lin * quaternion.conjugate()


def smart_act(quat, obj):
    """General purpose function for letting BiQuaternions act on object that detects
    object type.
    """
    if isinstance(obj, BiQuaternion):
        # Detect if it is a line
        if obj.scal.expand() == 0 and obj.primal() != 0:
            return act_on_line(quat, obj)
        else:
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
