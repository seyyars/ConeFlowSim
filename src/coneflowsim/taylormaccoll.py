import math
import numpy as np

GAMMA_DEFAULT = 1.4


def taylor_maccoll_rhs(theta: float,
                       y: np.ndarray,
                       gamma: float = GAMMA_DEFAULT) -> np.ndarray:
    """
    Right-hand side of the Taylor–Maccoll ODE in first-order form.

    Parameters
    ----------
    theta : float
        Polar angle (rad).
    y : array_like, shape (2,)
        State vector [F, F_theta], where
        F       = F(theta)
        F_theta = dF/dtheta.
    gamma : float
        Heat-capacity ratio.

    Returns
    -------
    dy_dtheta : ndarray, shape (2,)
        [dF/dtheta, dF_theta/dtheta].
    """
    F, F_theta = y

    # cot(theta) با مراقبت از تکینگی نزدیک θ = 0
    s = math.sin(theta)
    c = math.cos(theta)
    cot_theta = c / s if abs(s) > 1e-10 else 0.0

    # معادله‌ی Taylor–Maccoll
    A = 0.5 * (gamma + 1.0) * F_theta**2 - 0.5 * (gamma - 1.0) * (1.0 - F**2)
    B = ((gamma - 1.0) * (1.0 - F**2) * F
         + 0.5 * (gamma - 1.0) * cot_theta * (1.0 - F**2) * F_theta
         - gamma * F * F_theta**2
         - 0.5 * (gamma - 1.0) * cot_theta * F**3)

    dF_dtheta = F_theta

    # جلوگیری از blow-up وقتی A خیلی کوچک است
    if abs(A) < 1e-12:
        dFtheta_dtheta = 0.0
    else:
        dFtheta_dtheta = B / A

    return np.array([dF_dtheta, dFtheta_dtheta])


# ---------- حل عددی Taylor–Maccoll روی θ ----------

def _rk4_step(f, theta, y, h, gamma):
    """یک گام ساده RK4 برای y' = f(theta, y)."""
    k1 = f(theta,           y,             gamma)
    k2 = f(theta + 0.5*h,   y + 0.5*h*k1,  gamma)
    k3 = f(theta + 0.5*h,   y + 0.5*h*k2,  gamma)
    k4 = f(theta + h,       y + h*k3,      gamma)
    return y + (h/6.0) * (k1 + 2*k2 + 2*k3 + k4)


def solve_cone(M1: float,
               theta_deg: float,
               gamma: float = GAMMA_DEFAULT,
               n_theta: int = 400):
    """
    Solve the Taylor–Maccoll ODE from the axis (theta=0) to the cone surface.

    Educational 1D solver:

        F(0) = 1,  F'(0) = 0  (شرط تقارن روی محور)

    Parameters
    ----------
    M1 : float
        Free-stream Mach number.
    theta_deg : float
        Cone half-angle (degrees).
    gamma : float
        Heat-capacity ratio.
    n_theta : int
        Number of theta points.

    Returns
    -------
    theta : ndarray (n_theta,)
    F     : ndarray (n_theta,)
    F_th  : ndarray (n_theta,)
    M     : ndarray (n_theta,)
        Local Mach number along each ray.
    p_over_p1 : ndarray (n_theta,)
        Approximate p(theta)/p1 (ایزنتروپیک).
    """
    theta_c = math.radians(theta_deg)
    theta = np.linspace(0.0, theta_c, n_theta)

    # شرط اولیه روی محور
    y = np.array([1.0, 0.0])  # [F, F_theta]

    F_vals = np.zeros(n_theta)
    G_vals = np.zeros(n_theta)
    F_vals[0] = y[0]
    G_vals[0] = y[1]

    for i in range(n_theta - 1):
        h = theta[i+1] - theta[i]
        y = _rk4_step(taylor_maccoll_rhs, theta[i], y, h, gamma)
        F_vals[i+1] = y[0]
        G_vals[i+1] = y[1]

    # تبدیل F به Mach و p/p1 به صورت مدل آموزشی
    q = np.abs(F_vals)  # سرعت بی‌بعد ساده

    num = (M1**2) * q**2
    den = 1.0 + 0.5 * (gamma - 1.0) * (M1**2) * (1.0 - q**2)
    M_vals = np.sqrt(np.maximum(num / den, 1e-12))

    T_over_Tinf = 1.0 + 0.5 * (gamma - 1.0) * (M1**2) * (1.0 - q**2)
    p_over_pinf = T_over_Tinf**(gamma / (gamma - 1.0))

    # p1 را فعلاً همان p∞ می‌گیریم
    p_over_p1 = p_over_pinf

    return theta, F_vals, G_vals, M_vals, p_over_p1


# برای سازگاری با __init__.py
def taylor_maccoll_solution(M1: float,
                            theta_deg: float,
                            gamma: float = GAMMA_DEFAULT,
                            n_theta: int = 400):
    """Alias روی solve_cone برای API قدیمی‌تر."""
    return solve_cone(M1=M1, theta_deg=theta_deg,
                      gamma=gamma, n_theta=n_theta)


# ---------- ساخت میدان دوبعدی برای مقایسه با wedge ----------

def cone_field_tm(M1: float,
                  theta_deg: float,
                  L: float = 1.0,
                  nx: int = 200,
                  ny: int = 200,
                  gamma: float = GAMMA_DEFAULT):
    """
    Build a 2D Mach / p/p1 field for the Taylor–Maccoll cone solution.

    این تابع از حل ۱بعدی Taylor–Maccoll (solve_cone) استفاده می‌کند و
    آن را روی صفحه‌ی (x, r) نگاشت می‌کند.

    Parameters
    ----------
    M1 : float
        Free-stream Mach number.
    theta_deg : float
        Cone half-angle (degrees).
    L : float
        Normalized axial length.
    nx, ny : int
        Grid resolution in x (streamwise) and r (radial).
    gamma : float
        Heat-capacity ratio.

    Returns
    -------
    X, R : ndarray (ny, nx)
        Meshgrid of x and r.
    Mach : ndarray (ny, nx)
        Mach field (NaN خارج از مخروط).
    P : ndarray (ny, nx)
        p/p1 field (NaN خارج از مخروط).
    """
    theta_c = math.radians(theta_deg)

    # ۱) حل ۱بعدی روی θ
    theta_1d, _, _, M_1d, p_over_p1_1d = solve_cone(
        M1=M1, theta_deg=theta_deg, gamma=gamma, n_theta=600
    )

    # ۲) شبکه‌ی (x, r)
    x = np.linspace(0.0, L, nx)
    r = np.linspace(0.0, L * math.tan(theta_c), ny)
    X, R = np.meshgrid(x, r)

    # ۳) زاویه‌ی θ برای هر نقطه
    X_safe = np.where(X > 1e-8, X, 1e-8)
    TH = np.arctan2(R, X_safe)

    Mach = np.full_like(X, np.nan, dtype=float)
    P = np.full_like(X, np.nan, dtype=float)

    inside = TH <= theta_c + 1e-8  # داخل مخروط

    theta_min = theta_1d[0]
    theta_max = theta_1d[-1]
    t_norm = (TH[inside] - theta_min) / (theta_max - theta_min)
    t_norm = np.clip(t_norm, 0.0, 1.0)

    idx_f = t_norm * (len(theta_1d) - 1)
    i0 = np.floor(idx_f).astype(int)
    i1 = np.clip(i0 + 1, 0, len(theta_1d) - 1)
    w1 = idx_f - i0
    w0 = 1.0 - w1

    Mach_flat = w0 * M_1d[i0] + w1 * M_1d[i1]
    P_flat = w0 * p_over_p1_1d[i0] + w1 * p_over_p1_1d[i1]

    Mach[inside] = Mach_flat
    P[inside] = P_flat

    return X, R, Mach, P
