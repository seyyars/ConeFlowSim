# همچنان داخل taylormaccoll.py یا اگر دوست داری در فایل جدا:
from .wedge import wedge_field

def comparison_plot(
    M1: float = 3.0,
    theta_deg: float = 10.0,
    gamma: float = GAMMA_DEFAULT,
    L: float = 1.0,
    nx: int = 200,
    ny: int = 200,
):
    """
    Plot wedge-model vs Taylor–Maccoll model on a single figure.

    Left column: wedge model (Mach, p/p1)
    Right column: Taylor–Maccoll model (Mach, p/p1)
    """
    # Wedge field
    Xw, Yw, Mw, Pw, beta_w = wedge_field(
        M1, theta_deg, L=L, nx=nx, ny=ny, gamma=gamma
    )

    # Taylor–Maccoll field
    Xt, Yt, Mt, Pt = cone_field_tm(
        M1, theta_deg, L=L, nx=nx, ny=ny, gamma=gamma
    )

    theta_c = math.radians(theta_deg)

    fig, axs = plt.subplots(2, 2, figsize=(10, 8), sharex=True, sharey=True)

    # --- Wedge Mach ---
    im0 = axs[0, 0].pcolormesh(Xw, Yw, Mw, shading="nearest")
    axs[0, 0].plot([0, L], [0, math.tan(theta_c) * L], "k-", label="cone")
    axs[0, 0].plot([0, L], [0, math.tan(beta_w) * L], "r--", label="shock")
    axs[0, 0].set_title("Wedge model: Mach")
    axs[0, 0].set_ylabel("r (normalized)")
    fig.colorbar(im0, ax=axs[0, 0], label="Mach")

    # --- Wedge pressure ---
    im1 = axs[1, 0].pcolormesh(Xw, Yw, Pw, shading="nearest")
    axs[1, 0].plot([0, L], [0, math.tan(theta_c) * L], "k-")
    axs[1, 0].plot([0, L], [0, math.tan(beta_w) * L], "r--")
    axs[1, 0].set_title("Wedge model: p/p1")
    axs[1, 0].set_xlabel("x (normalized)")
    axs[1, 0].set_ylabel("r (normalized)")
    fig.colorbar(im1, ax=axs[1, 0], label="p/p1")

    # --- TM Mach ---
    im2 = axs[0, 1].pcolormesh(Xt, Yt, Mt, shading="nearest")
    axs[0, 1].plot([0, L], [0, math.tan(theta_c) * L], "k-", label="cone")
    # بتای TM را از params می‌شود بیرون کشید؛ فعلاً می‌توانی همان beta_w را استفاده کنی.
    axs[0, 1].set_title("Taylor–Maccoll: Mach")
    fig.colorbar(im2, ax=axs[0, 1], label="Mach")

    # --- TM pressure ---
    im3 = axs[1, 1].pcolormesh(Xt, Yt, Pt, shading="nearest")
    axs[1, 1].plot([0, L], [0, math.tan(theta_c) * L], "k-")
    axs[1, 1].set_title("Taylor–Maccoll: p/p1")
    axs[1, 1].set_xlabel("x (normalized)")
    fig.colorbar(im3, ax=axs[1, 1], label="p/p1")

    for ax in axs.ravel():
        ax.set_aspect("equal", adjustable="box")

    fig.suptitle(
        f"Supersonic flow over cone: wedge vs Taylor–Maccoll\n"
        f"M1 = {M1:.2f}, theta = {theta_deg:.1f}°",
        fontsize=12,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.95])

    return fig, axs
