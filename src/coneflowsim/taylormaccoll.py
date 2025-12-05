import math
import numpy as np

GAMMA_DEFAULT = 1.4


def taylor_maccoll_solution(
    M1: float,
    theta_c_deg: float,
    gamma: float = GAMMA_DEFAULT,
    n_theta: int = 200,
):
    """
    Placeholder Taylor–Maccoll solver.

    فعلاً فقط یک حل خیلی ساده و نمادین برمی‌گردانیم تا API کامل باشد.
    بعداً معادله‌ی واقعی Taylor–Maccoll را این‌جا پیاده‌سازی می‌کنیم.
    """
    theta_c = math.radians(theta_c_deg)

    # شبکه‌ی ساده از زاویه‌ها بین مخروط و محور
    theta = np.linspace(theta_c, 0.0, n_theta)

    # فعلاً F و dF/dtheta را صفر می‌گیریم (placeholder)
    F = np.zeros_like(theta)
    dF = np.zeros_like(theta)

    return theta, F, dF


def cone_field_tm(
    M1: float,
    theta_deg: float,
    L: float = 1.0,
    nx: int = 200,
    ny: int = 200,
    gamma: float = GAMMA_DEFAULT,
):
    """
    میدان «نمادین» Taylor–Maccoll روی شبکه‌ی (x, r).

    فعلاً برای این‌که شکل مقایسه‌ای داشته باشیم، از همان میدان wedge
    استفاده می‌کنیم؛ یعنی تا وقتی حل واقعی TM را پیاده نکرده‌ایم،
    دو ستون شکل‌ها تقریباً یکسان خواهند بود.
    """
    from .wedge import wedge_field

    X, Y, M, P, beta = wedge_field(
        M1, theta_deg, L=L, nx=nx, ny=ny, gamma=gamma
    )
    return X, Y, M, P
