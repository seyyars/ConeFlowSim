If you use ConeFlowSim in your work, please cite:

S. A. Asgari, ConeFlowSim (Version 1.0.1) [Computer software], Zenodo, 2025.
https://doi.org/10.5281/zenodo.17833232


# ConeFlowSim

Educational Python tools for supersonic flow over sharp cones.

ConeFlowSim couples a simple oblique-shock (wedge) approximation with a
Taylor–Maccoll conical flow solver and plotting utilities.  
The code is intentionally small and self-contained, suitable for teaching,
early–stage design studies, and reproducible research.

## Features

- Axisymmetric cone flow with:
  - Weak oblique-shock / wedge model
  - Taylor–Maccoll conical flow (explicit RK4 integrator)
- Mach and pressure fields between cone surface and shock
- Side–by–side comparison plots (wedge vs Taylor–Maccoll)
- Pure Python, small dependency set: `numpy`, `scipy`, `matplotlib`

Repository layout:

```text
ConeFlowSim/
├ paper/              # LaTeX article (coneflowsim-paper.tex)
├ src/coneflowsim/    # Python package (wedge + Taylor–Maccoll + plots)
│   ├ wedge.py
│   ├ taylormaccoll.py
│   ├ comparison_plot.py
│   └ __init__.py
├ tests/
│   └ test_tm_smoke.py  # basic smoke test for the TM solver
├ CITATION.cff        # citation metadata
├ refs.bib            # BibTeX entry for the Zenodo record
├ LICENSE             # MIT
└ README.md
