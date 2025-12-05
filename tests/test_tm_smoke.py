import pytest
import numpy as np
from coneflowsim.taylormaccoll import solve_cone, cone_field_tm
from coneflowsim import wedge_field


def test_solve_cone_runs():
    theta, F, G, M, p_over_p1 = solve_cone(M1=3.0, theta_deg=10.0)

    # طول آرایه‌ها معقول است
    assert len(theta) > 20
    assert len(M) == len(theta)

    # در ورودی، ماخ روی محور مخروط باید ~M1 باشد
    assert M[0] == pytest.approx(3.0, rel=1e-6)

    # در نزدیکی سطح مخروط، ماخ از ورودی کمتر می‌شود
    assert M[-1] < 3.0


def test_cone_field_tm_shape_matches_wedge():
    M1 = 3.0
    theta_deg = 10.0

    Xw, Yw, Mw, Pw, *_ = wedge_field(M1, theta_deg, L=1.0, nx=50, ny=100)
    Xt, Yt, Mt, Pt = cone_field_tm(M1, theta_deg, L=1.0, nx=50, ny=100)

    # شکل شبکه‌ی TM باید با wedge یکسان باشد
    assert Xt.shape == Xw.shape
    assert Yt.shape == Yw.shape
    assert Mt.shape == Mw.shape
    assert Pt.shape == Pw.shape
