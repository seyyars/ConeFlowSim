from .wedge import (
    oblique_shock,
    wedge_field,
    wedge_plot,
)

from .taylormaccoll import (
    taylor_maccoll_rhs,       # اگر دوست داری این را هم بیرون نشان بدهی
    taylor_maccoll_solution,
    cone_field_tm,
)

from .comparison_plot import (
    comparison_plot,
)

__all__ = [
    "oblique_shock",
    "wedge_field",
    "wedge_plot",
    "taylor_maccoll_rhs",      # اختیاری
    "taylor_maccoll_solution",
    "cone_field_tm",
    "comparison_plot",
]
