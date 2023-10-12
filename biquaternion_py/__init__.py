from .biquaternion import (
    define_algebra,
    BiQuaternion,
    II,
    JJ,
    KK,
    EE,
    act,
    act_on_line,
    fiber_project,
    inner,
    line,
    point,
    outer,
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
