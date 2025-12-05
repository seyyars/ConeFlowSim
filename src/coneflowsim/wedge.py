import math
from typing import Tuple

import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize

GAMMA_DEFAULT: float = 1.4


def oblique_shock(
    M1: float,
    theta: float,
    gamma: float = GAMMA_DEFAULT,
) -> Tuple[float, float, float, float]:
    """
    Weak oblique shock for given upstream Mach M1 and deflection angle theta (rad).

    Parameters
    ----------
    M1 : float
        Upstream Mach number.
    theta : float
        Flow deflection angle (radians).
    gamma : float
        Heat-capacity ratio.

    Returns
    -------
    beta : float
        Shock angle (radians).
    p2_p1 : float
        Static pressure ratio across the shock.
    T2_T1 : float
        Static temperature ratio across the shock.
    M2 : float
        Downstream Mach number immediately behind the shock.
    """

    def f(beta: float) -> float:
        """Theta–beta–Mach relation: f(beta) = tan(theta) − RHS(beta)."""
        Mn1_sq = (M1 * math.sin(beta)) ** 2
        num = 2.0 * (Mn1_sq - 1.0)
        den = M1**2 * (gamma + math.cos(2.0 * beta)) + 2.0
        rhs = num / (den * math.tan(beta))  # = tan(theta_beta)
        return math.tan(theta) - rhs

    # --- robust bracketing for the weak solution in (theta, pi/2) -----------
    beta_min = theta + math.radians(0.1)          # کمی بالاتر از theta
    beta_max = 0.5 * math.pi - math.radians(0.1)  # کمی کمتر از 90 درجه

    betas = np.linspace(beta_min, beta_max, 400)
    vals = [f(b) for b in betas]

    beta_lo = None
    beta_hi = None
    for i in range(len(betas) - 1):
        v1, v2 = vals[i], vals[i + 1]
        if math.isnan(v1) or math.isnan(v2):
            continue
        if v1 == 0.0:
            beta_lo = betas[i]
            beta_hi = betas[i + 1]
            break
        if v1 * v2 < 0.0:
            beta_lo = betas[i]
            beta_hi = betas[i + 1]
            break

    if beta_lo is None or beta_hi is None:
        raise ValueError(
            f"Could not bracket oblique-shock angle for "
            f"M1={M1}, theta={theta:.4f} rad, gamma={gamma}."
        )

    beta = optimize.brentq(f, beta_lo, beta_hi)

    # --- Normal-shock relations ---------------------------------------------
    Mn1 = M1 * math.sin(beta)
    Mn1_sq = Mn1**2

    p2_p1 = 1.0 + 2.0 * gamma / (gamma + 1.0) * (Mn1_sq - 1.0)
    rho2_rho1 = (gamma + 1.0) * Mn1_sq / ((gamma - 1.0) * Mn1_sq + 2.0)
    T2_T1 = p2_p1 / rho2_rho1

    Mn2_sq = (1.0 + 0.5 * (gamma - 1.0) * Mn1_sq) / (
        gamma * Mn1_sq - 0.5 * (gamma - 1.0)
    )
    Mn2 = math.sqrt(Mn2_sq)
    M2 = Mn2 / math.sin(beta - theta)

    return beta, p2_p1, T2_T1, M2


def wedge_field(
    M1: float,
    theta_deg: float,
    L: float = 1.0,
    nx: int = 200,
    ny: int = 200,
    gamma: float = GAMMA_DEFAULT,
):
    """
    Build a simple Mach/pressure field between cone surface and shock
    using oblique-shock (wedge) approximation.

    Returns
    -------
    X, Y : 2D ndarray
        Mesh in (x, r) coordinates (normalized).
    M, P : 2D ndarray
        Mach and p/p1 fields. Points outside the cone–shock region are NaN.
    beta : float
        Shock angle (rad).
    p2_p1 : float
        Pressure ratio across the shock.
    T2_T1 : float
        Temperature ratio across the shock.
    M2 : float
        Downstream Mach number immediately behind the shock.
    """
    theta = math.radians(theta_deg)

    # Solve oblique shock (weak solution)
    beta, p2_p1, T2_T1, M2 = oblique_shock(M1, theta, gamma=gamma)

    # Grid in (x, r)
    x = np.linspace(0.0, L, nx)
    r_max = math.tan(beta) * L
    r = np.linspace(0.0, r_max, ny)
    X, Y = np.meshgrid(x, r, indexing="xy")  # (ny, nx)

    # Polar angle for each point
    with np.errstate(divide="ignore", invalid="ignore"):
        phi = np.arctan2(Y, X)

    # Initialize fields
    M = np.full_like(X, np.nan, dtype=float)
    P = np.full_like(X, np.nan, dtype=float)

    # Region between cone and shock
    inside = (X > 0.0) & (phi >= theta) & (phi <= beta)

    # Simple wedge model: uniform post-shock state between cone and shock
    M[inside] = M2
    P[inside] = p2_p1

    return X, Y, M, P, beta, p2_p1, T2_T1, M2


def wedge_plot(M1: float = 3.0, theta_deg: float = 10.0, L: float = 1.0):
    """
    Plot Mach and p/p1 fields for the wedge model.
    """
    X, Y, Mach, P, beta, p2_p1, T2_T1, M2 = wedge_field(M1, theta_deg, L=L)
    theta = math.radians(theta_deg)

    fig, axs = plt.subplots(1, 2, figsize=(10, 4), sharey=True)

    # Mach field
    im0 = axs[0].pcolormesh(X, Y, Mach, shading="nearest")
    axs[0].plot([0, L], [0, math.tan(theta) * L], "k-", label="cone")
    axs[0].plot([0, L], [0, math.tan(beta) * L], "r--", label="shock")
    axs[0].set_title("Mach (wedge model)")
    axs[0].set_xlabel("x (normalized)")
    axs[0].set_ylabel("r (normalized)")
    fig.colorbar(im0, ax=axs[0], label="Mach")
    axs[0].legend(loc="upper left")

    # Pressure field
    im1 = axs[1].pcolormesh(X, Y, P, shading="nearest")
    axs[1].plot([0, L], [0, math.tan(theta) * L], "k-", label="cone")
    axs[1].plot([0, L], [0, math.tan(beta) * L], "r--", label="shock")
    axs[1].set_title("p/p1 (wedge model)")
    axs[1].set_xlabel("x (normalized)")
    fig.colorbar(im1, ax=axs[1], label="p/p1")

    fig.suptitle(
        f"Supersonic flow over cone (wedge approximation)\n"
        f"M1 = {M1:.2f}, theta = {theta_deg:.1f}°"
    )
    fig.tight_layout(rect=[0, 0, 1, 0.90])

    return fig, axs

