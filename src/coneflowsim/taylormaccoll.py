# src/coneflowsim/taylormaccoll.py

"""
Taylor–Maccoll cone-flow model for ConeFlowSim.

این ماژول یک حل عددی برای معادله‌ی Taylor–Maccoll فراهم می‌کند و
میدان Mach / فشار را در جلوی یک مخروط تیز تولید می‌کند.
برای این‌که مدل فیزیکی کاملاً درست باشد، باید RHS معادله‌ی
Taylor–Maccoll را در تابع `taylor_maccoll_rhs` به‌درستی وارد کنید.

مراجع پیشنهادی:
- J. D. Anderson, "Modern Compressible Flow", chapter on conical flow.
- Wikipedia page: "Taylor–Maccoll equation".
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
import matplotlib.pyplot as plt

from .wedge import oblique_shock, wedge_field

GAMMA_DEFAULT = 1.4


@dataclass
class TMParams:
    gamma: float
    M1: float
    theta_c: float      # cone half-angle (rad)
    beta: float         # shock angle (rad)
    a1: float           # upstream sound speed (normalized)
    V1: float           # upstream speed  (normalized)
    V2r: float          # post-shock radial component at theta=beta
    V2t: float          # post-shock tangential component at theta=beta


# ---------- 1) معادله‌ی دیفرانسیل Taylor–Maccoll (TODO: پر کردن RHS) ----------

def taylor_maccoll_rhs(theta: float, y: np.ndarray, params: TMParams) -> np.ndarray:
    """
    Right-hand side of the Taylor–Maccoll equations in terms of
    radial and polar velocity components (V_r, V_theta).

    Parameters
    ----------
    theta : float
        Polar angle (rad) measured from the cone axis.
    y : array_like, shape (2,)
        State vector [V_r, V_theta] در دستگاه کروی (نرمال‌سازی شده مثلاً با a1).
    params : TMParams
        Flow / gas parameters.

    Returns
    -------
    dydtheta : ndarray, shape (2,)
        [dV_r/dtheta, dV_theta/dtheta]

    NOTE
    ----
    در حال حاضر این تابع فقط اسکلت است و باید با فرم دقیق معادله‌ی
    Taylor–Maccoll از مرجع دلخواهت پر شود.
    """
    Vr, Vt = y
    g = params.gamma

    # بزرگی سرعت و عدد ماخ محلی (در صورت نیاز)
    V2 = Vr * Vr + Vt * Vt

    # ------------ TODO: اینجا RHS واقعی معادله را بگذار ------------
    # جای زیر فقط یک «placeholder» است تا کد خراب نشود؛
    # آن را با معادله‌ی درست جایگزین کن.
    dVr_dtheta = Vt
    dVt_dtheta = 0.0  # <-- این باید فرم کامل Taylor–Maccoll باشد

    return np.array([dVr_dtheta, dVt_dtheta], dtype=float)


# ---------- 2) محاسبه‌ی شرایط اولیه پشت موج شوک ----------

def post_shock_state(M1: float, theta_c: float, gamma: float = GAMMA_DEFAULT) -> TMParams:
    """
    از مدل گوه‌ای برای تخمین زاویه‌ی شوک و وضعیت اولیه پشت شوک استفاده می‌کند.

    ما از تابع oblique_shock موجود در wedge.py استفاده می‌کنیم تا
    زاویه‌ی شوک β و عدد ماخ پشت شوک M2 را به‌دست آوریم، و سپس سرعت
    برداری را به مؤلفه‌های کروی Vr, Vtheta در θ = β تبدیل می‌کنیم.

    این کار تقریب است (چون مخروط، گوه‌ی دو بعدی نیست) ولی برای
    مدل آموزشی کفایت می‌کند.
    """
    theta = theta_c
    beta, p2_p1, T2_T1, M2 = oblique_shock(M1, theta, gamma=gamma)

    # فرض کنیم a1 = 1 (نرمال‌سازی): V1 = M1 * a1
    a1 = 1.0
    V1 = M1 * a1

    # جهت جریان پشت شوک بعد از انحراف تقریباً برابر θ_c است
    # (مثل مسئله‌ی گوه). در کروی، Vr در راستای شعاع از رأس، و Vt در
    # راستای افزایش θ است؛ در θ = β جهت سرعت تقریباً در زاویه θ_c است.
    flow_angle = theta_c  # نسبت به محور

    V2 = M2 * a1
    Vr2 = V2 * math.cos(flow_angle)
    Vt2 = -V2 * math.sin(flow_angle)

    return TMParams(
        gamma=gamma,
        M1=M1,
        theta_c=theta_c,
        beta=beta,
        a1=a1,
        V1=V1,
        V2r=Vr2,
        V2t=Vt2,
    )


# ---------- 3) انتگرال‌گیری عددی روی θ: حل Taylor–Maccoll ----------

def taylor_maccoll_solution(
    M1: float,
    theta_deg: float,
    gamma: float = GAMMA_DEFAULT,
    n_steps: int = 800,
) -> tuple[np.ndarray, np.ndarray, TMParams]:
    """
    Integrate the Taylor–Maccoll ODE from shock angle beta down to cone angle theta_c.

    Parameters
    ----------
    M1 : float
        Upstream Mach number.
    theta_deg : float
        Cone half-angle in degrees.
    gamma : float, optional
        Heat capacity ratio.
    n_steps : int, optional
        Number of theta grid points between beta and theta_c.

    Returns
    -------
    theta : ndarray
        Array of polar angles from cone surface (theta_c) up to shock (beta).
    sol : ndarray, shape (2, N)
        Velocities [Vr(theta), Vtheta(theta)].
    params : TMParams
        Parameters (includes beta etc.).
    """
    theta_c = math.radians(theta_deg)
    params = post_shock_state(M1, theta_c, gamma)

    beta = params.beta
    theta_grid = np.linspace(beta, theta_c, n_steps)

    Vr = np.empty_like(theta_grid)
    Vt = np.empty_like(theta_grid)

    # شرط اولیه در θ = β
    Vr[0] = params.V2r
    Vt[0] = params.V2t

    # انتگرال‌گیر ساده‌ی Runge–Kutta مرتبه‌ی 4
    for i in range(len(theta_grid) - 1):
        th = theta_grid[i]
        h = theta_grid[i + 1] - theta_grid[i]
        y = np.array([Vr[i], Vt[i]])

        k1 = taylor_maccoll_rhs(th, y, params)
        k2 = taylor_maccoll_rhs(th + 0.5 * h, y + 0.5 * h * k1, params)
        k3 = taylor_maccoll_rhs(th + 0.5 * h, y + 0.5 * h * k2, params)
        k4 = taylor_maccoll_rhs(th + h, y + h * k3, params)

        y_new = y + (h / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
        Vr[i + 1], Vt[i + 1] = y_new

    return theta_grid, np.vstack([Vr, Vt]), params


# ---------- 4) ساخت میدان Mach / فشار از روی حل TM ----------

def cone_field_tm(
    M1: float,
    theta_deg: float,
    L: float = 1.0,
    nx: int = 200,
    ny: int = 200,
    gamma: float = GAMMA_DEFAULT,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Build a 2D Mach and pressure field in the (x, r) plane using
    the Taylor–Maccoll solution.

    بر روی شبکه‌ای در مستطیل 0<=x<=L، 0<=r<=r_max میدان را می‌سازیم و
    هر نقطه را به نزدیک‌ترین θ در حل Taylor–Maccoll نگاشت می‌کنیم.
    """
    theta, sol, params = taylor_maccoll_solution(M1, theta_deg, gamma=gamma)
    Vr, Vt = sol
    Vmag = np.sqrt(Vr**2 + Vt**2)
    # فرض a1=1: Mach = V/a1
    Mach_theta = Vmag

    # فشار نسبی را با استفاده از رابطه‌ی آیزنتروپیک از Mach محلی تخمین می‌زنیم
    g = gamma
    p_over_p1_theta = (1.0 + 0.5 * (g - 1.0) * params.M1**2) / (
        1.0 + 0.5 * (g - 1.0) * Mach_theta**2
    )
    p_over_p1_theta **= g / (g - 1.0)

    theta_c = math.radians(theta_deg)
    beta = params.beta

    x = np.linspace(0.0, L, nx)
    r = np.linspace(0.0, L, ny)
    X, R = np.meshgrid(x, r)

    # زاویه‌ی هر نقطه نسبت به محور
    phi = np.arctan2(R, X + 1e-12)

    Mach = np.full_like(X, np.nan, dtype=float)
    P = np.full_like(X, np.nan, dtype=float)

    # ناحیه بین مخروط و شوک: از TM استفاده می‌کنیم
    between = (phi >= theta_c) & (phi <= beta)
    idx = between.nonzero()
    if idx[0].size:
        # برای هر نقطه نزدیک‌ترین θ را پیدا می‌کنیم
        phi_vals = phi[idx]
        theta_idx = np.searchsorted(theta[::-1], phi_vals)
        theta_idx = np.clip(theta_idx, 0, len(theta) - 1)
        theta_idx = len(theta) - 1 - theta_idx  # چون معکوس کرده‌ایم

        Mach[idx] = Mach_theta[theta_idx]
        P[idx] = p_over_p1_theta[theta_idx]

    return X, R, Mach, P
