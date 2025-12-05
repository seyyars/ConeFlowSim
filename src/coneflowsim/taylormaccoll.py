import math
import numpy as np

GAMMA_DEFAULT = 1.4


def taylor_maccoll_rhs(theta: float,
                       y: np.ndarray,
                       gamma: float = GAMMA_DEFAULT) -> np.ndarray:
    """
    Right-hand side of Taylor–Maccoll ODE in first-order form.

    فعلاً فقط یک RHS بسیار ساده (placeholder) می‌دهیم تا اسم تابع‌ها درست باشند.
    بعداً فرمول دقیق را این‌جا می‌گذاریم.
    """
    F, F_theta = y
    # Placeholder: مشتق دوم را صفر می‌گیریم
    F_theta_theta = 0.0
    return np.array([F_theta, F_theta_theta], dtype=float)


def taylor_maccoll_solution(M1: float,
                            theta_c_deg: float,
                            gamma: float = GAMMA_DEFAULT,
                            n_theta: int = 200):
    """
    Placeholder Taylor–Maccoll solution.

    خروجی فقط برای این است که API کامل باشد؛ فعلاً F(θ) و F'(θ) را صفر می‌گیریم.
    """
    theta_c = math.radians(theta_c_deg)

    # grid ساده بین سطح مخروط و محور
    theta = np.linspace(theta_c, 0.0, n_theta)

    F = np.zeros_like(theta)
    F_theta = np.zeros_like(theta)

    return theta, F, F_theta


def cone_field_tm(M1: float,
                  theta_deg: float,
                  L: float = 1.0,
                  nx: int = 200,
                  ny: int = 200,
                  gamma: float = GAMMA_DEFAULT):
    """
    Taylor–Maccoll cone field (placeholder).

    فعلاً برای این‌که شکل مقایسه‌ای داشته باشیم، از همان میدان wedge استفاده
    می‌کنیم؛ یعنی تا زمانی که حل واقعی TM را پیاده نکرده‌ایم، دو ستون شکل‌ها
    تقریباً یکسان خواهند بود.
    """
    from .wedge import wedge_field

    X, Y, M, P, beta = wedge_field(
        M1, theta_deg, L=L, nx=nx, ny=ny, gamma=gamma
    )
    return X, Y, M, P
