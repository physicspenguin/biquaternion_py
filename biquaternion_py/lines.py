from .biquaternion import BiQuaternion


def quaternion_to_pluecker(quat):
    """Convert BiQuaternion to Pluecker coordinates."""
    return [*quat.coeffs[1:4], -quat.coeffs[5], -quat.coeffs[6], -quat.coeffs[7]]


def pluecker_to_quaternion(line):
    """Convert Pluecker coordinates to BiQuaternion."""
    return BiQuaternion([0, *line[0:3], 0, *[-line[i + 3] for i in range(3)]])
