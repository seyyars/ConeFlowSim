import math

import matplotlib.pyplot as plt

from .wedge import wedge_field
from .taylormaccoll import cone_field_tm, GAMMA_DEFAULT


def comparison_plot(
    M1: float = 3.0,
    theta_deg: float = 10.0,
    gamma: float = GAMMA_DEFAULT,
    L: float = 1.0,
    nx: int = 200,
    ny: int = 200,
):
    """
    Plot wedge vs Taylor–Maccoll fields side by side.

    Left column : wedge model (Mach, p/p1)
    Right column: Taylor–Maccoll model (Mach, p/p1)

    Parameters
    ----------
    M1 : float
        Upstream Mach number.
    theta_deg : float
        Cone half-angle in degrees.
    gamma : float
        Heat capacity ratio.
    L : float
        Axial length (normalized).
    nx, ny : int
        Grid resolution in x and r.
    """
    # --- Wedge field ---------------------------------------------------------
    (
        Xw,
        Yw,
        Mw,
        Pw,
        beta_w,
        p2_p1_w,
        T2_T1_w,
        M2_w,
    ) = wedge_field(M1, theta_deg, L=L, nx=nx, ny=ny, gamma=gamma)

    # --- Taylor–Maccoll field (2D mapped) -----------------------------------
    # cone_field_tm باید یک شبکه‌ی (x, r) و میدان Mach و p/p1 بدهد.
    Xt, Yt, Mt, p_over_p1_t = cone_field_tm(
        M1, theta_deg, L=L, nx=nx, ny=ny, gamma=gamma
    )

    theta_c = math.radians(theta_deg)

    fig, axs = plt.subplots(2, 2, figsize=(10, 8), sharex=True, sharey=True)

    # -------------------------------------------------------------------------
    # Mach (wedge)
    # -------------------------------------------------------------------------
    im0 = axs[0, 0].pcolormesh(Xw, Yw, Mw, shading="nearest")
    axs[0, 0].plot([0, L], [0, math.tan(theta_c) * L], "k-", label="cone")
    axs[0, 0].plot([0, L], [0, math.tan(beta_w) * L], "r--", label="shock")
    axs[0, 0].set_title("Mach (wedge)")
    axs[0, 0].set_ylabel("r (normalized)")
    fig.colorbar(im0, ax=axs[0, 0], label="Mach")
    axs[0, 0].legend(loc="upper left")

    # -------------------------------------------------------------------------
    # Mach (Taylor–Maccoll)
    # -------------------------------------------------------------------------
    im1 = axs[0, 1].pcolormesh(Xt, Yt, Mt, shading="nearest")
    axs[0, 1].plot([0, L], [0, math.tan(theta_c) * L], "k-", label="cone")
    axs[0, 1].set_title("Mach (Taylor–Maccoll)")
    fig.colorbar(im1, ax=axs[0, 1], label="Mach")

    # -------------------------------------------------------------------------
    # p/p1 (wedge)
    # -------------------------------------------------------------------------
    # --- p/p1 (wedge) ---
    im2 = axs[1, 0].pcolormesh(Xw, Yw, Pw, shading="nearest")  
    axs[1, 0].plot([0, L], [0, math.tan(theta_c) * L], "k-", label="cone")
    axs[1, 0].plot([0, L], [0, math.tan(beta_w) * L], "r--", label="shock")
    axs[1, 0].set_title("p/p1 (wedge)")
    fig.colorbar(im2, ax=axs[1, 0], label="p/p1")

    # -------------------------------------------------------------------------
    # p/p1 (Taylor–Maccoll)
    # -------------------------------------------------------------------------
    im3 = axs[1, 1].pcolormesh(Xt, Yt, p_over_p1_t, shading="nearest")
    axs[1, 1].plot([0, L], [0, math.tan(theta_c) * L], "k-", label="cone")
    axs[1, 1].set_title("p/p1 (Taylor–Maccoll)")
    axs[1, 1].set_xlabel("x (normalized)")
    fig.colorbar(im3, ax=axs[1, 1], label="p/p1")

    # ظاهر کلی
    for ax in axs.ravel():
        ax.set_aspect("equal", adjustable="box")

    fig.suptitle(
        f"Supersonic cone flow: wedge vs Taylor–Maccoll\n"
        f"M1 = {M1:.2f}, theta = {theta_deg:.1f}°",
        fontsize=12,
    )
    fig.tight_layout(rect=[0, 0.0, 1, 0.93])

    return fig, axs
