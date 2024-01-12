"""Module for line creation, representation and manipulation."""
import numpy as np
from .biquaternion import BiQuaternion


def quat_to_pluecker(quat):
    """Generate pluecker coordinates from a given BiQuaternion."""
    if quat.scal == 0 and quat.eps == 0:
        return [*quat.coeffs[1:4], *[-quat.coeffs[5 + i] for i in range(3)]]
    else:
        raise ValueError("Object not a valid description of a line")


def pluecker_to_quat(coord):
    """Generate BiQuaternion reperesentation of a line with specified
    Pluecker coordinates."""
    return BiQuaternion([0, *coord[0:3], 0, *[-coord[i + 3] for i in range(3)]])


def line_to_pluecker(direction, point):
    """Generate Pluecker coordinates for a line with given direction through a point."""
    return [*direction, -np.cross(direction, point)]


def act_on_line(quaternion, lin):
    """Let a BiQuaternion act on a line."""
    return quaternion * lin * quaternion.conjugate()
