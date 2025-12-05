import math
import numpy as np

GAMMA_DEFAULT = 1.4


def taylor_maccoll_rhs(theta: float,
                       y: np.ndarray,
                       gamma: float = GAMMA_DEFAULT) -> np.ndarray:
    """
    Right-hand side of Taylor–Maccoll ODE in first-order form.

    Parameters
    ----------
    theta : float
        Polar angle (rad).
    y : array_like, shape (2,)
        State vector [F, F_theta], where F_theta = dF/dtheta.
    gamma : float
        Heat capacity ratio.

    Returns
    -------
    dydtheta : ndarray, shape (2,)
        [dF/dtheta, d^2F/dtheta^2].
    """
    F, Fp = y  # F, F'

    # coefficients
    A = 0.5 * (gamma + 1.0) * Fp**2 - 0.5 * (gamma - 1.0) * (1.0 - F**2)
    cot_th = 1.0 / math.tan(theta)

    RHS = (
        (gamma - 1.0) * (1.0 - F**2) * F
        + 0.5 * (gamma - 1.0) * cot_th * (1.0 - F**2) * Fp
        - gamma * F * Fp**2
        - 0.5 * (gamma - 1.0) * cot_th * Fp**3
    )

    # dF/dtheta = F'
    dF_dtheta = Fp

    # dF'/dtheta = F'' = RHS / A
    dFp_dtheta = RHS / A

    return np.array([dF_dtheta, dFp_dtheta])
