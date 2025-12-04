import math
import numpy as np
import matplotlib.pyplot as plt

GAMMA_DEFAULT = 1.4


def oblique_shock(M1, theta, gamma=GAMMA_DEFAULT):
    """
    Compute weak oblique shock for a wedge with deflection angle theta (rad).

    Returns
    -------
    beta  : shock angle (rad)
    p2_p1 : pressure ratio
    T2_T1 : temperature ratio
    M2    : downstream Mach number
    """
    # اینجا فرمول/روت نیوتنی که در نوت‌بوک برای beta نوشته بودی را قرار بده
    # ...
    raise NotImplementedError  # این خط را بعد از کپی‌کردن کدت پاک کن


def wedge_field(M1, theta_deg, L=1.0, nx=200, ny=200, gamma=GAMMA_DEFAULT):
    """
    Build a simple Mach/pressure field between cone surface and shock
    using oblique-shock (wedge) approximation.

    بر پایه‌ی همان make_field قبلی است:
    - شبکه‌ی (x, r) می‌سازی
    - برای نقاط بین سطح مخروط و موج، M و p/p1 را می‌گذاری
    - خارج از آن NaN می‌گذاری
    """
    # کد مربوط به محاسبه‌ی Mach و p/p1 روی شبکه را از نوت‌بوک اینجا منتقل کن
    # ...
    raise NotImplementedError


def wedge_plot(M1=3.0, theta_deg=10.0, L=1.0):
    """
    Plot Mach and p/p1 fields (وضعیتی که در شکل دو تایی رنگی کشیدی).
    """
    X, Y, Mach, P, beta, p2_p1, T2_T1, M2 = wedge_field(M1, theta_deg, L=L)

    fig, axs = plt.subplots(1, 2, figsize=(10, 4), sharey=True)

    # Mach field
    m = axs[0].pcolormesh(X, Y, Mach, shading="nearest")
    axs[0].plot([0, L], [0, math.tan(math.radians(theta_deg)) * L], "k-", label="cone")
    axs[0].plot([0, L], [0, math.tan(beta) * L], "r--", label="shock")
    axs[0].set_title("Mach field")
    axs[0].set_xlabel("x (normalized)")
    axs[0].set_ylabel("r (normalized)")
    axs[0].legend()
    fig.colorbar(m, ax=axs[0], label="Mach")

    # Pressure field
    p = axs[1].pcolormesh(X, Y, P, shading="nearest")
    axs[1].plot([0, L], [0, math.tan(math.radians(theta_deg)) * L], "k-", label="cone")
    axs[1].plot([0, L], [0, math.tan(beta) * L], "r--", label="shock")
    axs[1].set_title("Pressure field p/p1")
    axs[1].set_xlabel("x (normalized)")
    axs[1].legend()
    fig.colorbar(p, ax=axs[1], label="p/p1")

    fig.suptitle(
        f"Supersonic flow over cone (wedge model)\n"
        f"M1 = {M1:.2f}, theta = {theta_deg:.1f}°"
    )
    fig.tight_layout()
    return fig, axs

