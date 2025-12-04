"""
Taylor–Maccoll cone-flow model for ConeFlowSim.

در حال حاضر فقط اسکلت تابع‌ها را داریم؛
در مرحله‌ی بعد معادله‌ی تیلور–مک‌کُل و رسم میدان را اینجا پیاده‌سازی می‌کنیم.
"""

import math
import numpy as np
import matplotlib.pyplot as plt

GAMMA_DEFAULT = 1.4


def taylor_maccoll_solution(M1, theta_c_deg, gamma=GAMMA_DEFAULT):
    """
    Placeholder for Taylor–Maccoll ODE solver.

    Parameters
    ----------
    M1 : float
        Free-stream Mach number.
    theta_c_deg : float
        Cone half-angle in degrees.
    gamma : float
        Ratio of specific heats.

    Returns
    -------
    dict
        For now just an empty dict. In the next step we'll return
        profiles of Mach, pressure etc. along the cone.
    """
    # TODO: implement Taylor–Maccoll solver
    return {}


def cone_field_tm(M1, theta_c_deg, L=1.0, nx=200, ny=200, gamma=GAMMA_DEFAULT):
    """
    Placeholder for 2D (r–x) field based on Taylor–Maccoll solution.
    Right now فقط یک شبکه‌ی خالی برمی‌گرداند.
    """
    x = np.linspace(0.0, L, nx)
    r = np.linspace(0.0, L, ny)
    X, R = np.meshgrid(x, r)

    Mach = np.full_like(X, np.nan, dtype=float)
    P    = np.full_like(X, np.nan, dtype=float)

    return X, R, Mach, P


def comparison_plot(M1=3.0, theta_deg=10.0, L=1.0):
    """
    Placeholder برای شکل مقایسه‌ی wedge و Taylor–Maccoll.
    فعلاً فقط یک شکل خالی می‌سازد تا import خراب نشود.
    """
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.text(0.5, 0.5,
            "Taylor–Maccoll model\n(not implemented yet)",
            ha="center", va="center")
    ax.set_axis_off()
    fig.tight_layout()
    return fig, ax
