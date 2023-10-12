from .biquaternion import (
    define_algebra,
    BiQuaternion,
    II,
    JJ,
    KK,
    EE,
)
from .biquat_tools import (
    point_to_quat,
    quat_to_point,
    hom_point_to_quat,
    quat_to_hom_point,
    pluecker_to_quat,
    quat_to_pluecker,
    line_to_pluecker,
    act_on_point,
    act_on_line,
    smart_act,
    inner,
    outer,
    fiber_project,
)
from .polynomials import Poly, poly_div
from .lines import quaternion_to_pluecker, pluecker_to_quaternion
from .poly_tools import (
    max_real_poly_fact,
    gcd_conj_pd,
    is_poly_reduced,
    factorize_bq_poly,
    factorize_from_list,
    split_lin_factor,
    irreducible_factors,
    # is_poly_real,
)
from .random_gen import rand_bq, rand_line, rand_quat, rand_rational
