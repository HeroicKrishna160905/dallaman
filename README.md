# dallaman
working currently on dalla man model.

Parameter Identifiability Test.
Reduction of Model Order.
Parameter Estimation.

Research repository: implementation and study of the Dalla Man meal simulation model for glucose–insulin dynamics.

currently work-in-progress

## Abstract (short)

This repository implements the canonical Dalla Man glucose–insulin meal simulation model and a compact reduced-order model (ROM). The project focuses on numerical simulation, parameter estimation, structural identifiability analysis (via GenSSI), and model-order reduction. The codebase contains Python implementations for open-loop simulations, plotting utilities, and a GenSSI-ready MATLAB symbolic model for identifiability analysis.

---

## Scientific background

The Dalla Man model (Dalla Man et al., 2007) is a physiologically-based compartmental model describing glucose absorption, distribution, and insulin dynamics following a meal. It is widely used for in silico meal experiments, controller design, and identifiability/parameter estimation studies. Key model features include:
- A multi-compartment description of glucose in plasma and tissue and gut compartments for meal absorption.
- Insulin dynamics including portal and peripheral compartments, insulin secretion and action.
- Nonlinear processes for gastric emptying and glucose appearance.

Reference:
Dalla Man, C., Rizza, R. A., & Cobelli, C. (2007). Meal simulation model of the glucose–insulin system. IEEE Trans. Biomed. Eng., 54(10), 1740–1749. https://doi.org/10.1109/TBME.2007.893506

---

## What is in this repository (file-level summary)

- `models/dallaman_openloop.py` — Primary Python implementation of the full Dalla Man model (12 states). Contains:
  - `dallaman_ode(t, x, p)` — state derivatives returning d/dt of the 12-state vector.
  - `run_simulation(t_span=(0,600), dt=0.1)` — builds parameter dictionary `p`, initial condition `x0`, integrates using `scipy.integrate.solve_ivp` (RK45 by default) and returns the solution.
  - `plot_solution(sol)` — plotting utility for plasma/tissue glucose and insulin traces.

- `dallaman_openloop.py` (root) — a near-duplicate or convenience script to run the full model from repository root with the same functions and example `if __name__ == '__main__'` execution.


- `rom_openloop.py` — Reduced-Order Model (ROM) implementation (3-state: glucose G, insulin I, insulin action X). Contains:
  - `rom_ode(t, x, p, Ra_func, u_I)` — ROM ODEs, includes endogenous insulin secretion law and input functions for glucose appearance (Ra) and insulin input.
  - `gastric_operator(t, D, ...)` — gamma-distribution shaped gastric appearance operator used to construct Ra(t).
  - `run_simulation(...)` and `plot_solution(sol)` utilities.

- `dallaman_genssi.m` — MATLAB symbolic model prepared for GenSSI structural identifiability analysis. It declares symbolic states, initials, and maps the full Dalla Man equations into the GenSSI-compatible structure. The file sets `model.sym.x`, `model.sym.x0`, and `model.sym.Nder` and exposes model state names required by GenSSI.

- Jupyter Notebooks (.ipynb) — the repository composition indicates most content is notebooks; they contain analysis/experiments.

---

## Key implementation details (for reproducibility)

State vector in the full Dalla Man implementation (12 states):
- x = [Gp, Gt, Il, Ip, Qsto1, Qsto2, Qgut, I1, Id, X, I_po, Y]
  - Gp: plasma glucose
  - Gt: tissue glucose
  - Il: insulin in liver
  - Ip: plasma insulin
  - Qsto1, Qsto2: stomach compartments for transient gastric emptying
  - Qgut: gut compartment (glucose available for absorption)
  - I1, Id: intermediary insulin compartments (delays)
  - X: insulin action (remote)
  - I_po: portal insulin
  - Y: auxiliary variable related to glucose sensor or dynamics in model formulation

Numerical solver and tolerances:
- Integration: scipy.integrate.solve_ivp with method='RK45'
- Typical tolerances used: atol=1e-6, rtol=1e-6 (set in `run_simulation`)

Example parameter set
- The model parameters are assembled into a Python dictionary `p` inside `run_simulation`. Example entries (non-exhaustive — full list is in file):
  - V_G = 1.88
  - k_1 = 0.065
  - k_2 = 0.079
  - G_b = 95.0
  - V_I = 0.05
  - m_1 = 0.19, m_2 = 0.484, m_4 = 0.194, m_5 = 0.0304, m_6 = 0.6471
  - HE_b = 0.6, I_b = 25.0
  - k_max = 0.0558, k_min = 0.008, k_abs = 0.057, k_gri = 0.0558
  - f = 0.9, b = 0.82, d = 0.01
  - BW = 78.0
  - k_p1 = 2.7, k_p2 = 0.0021, k_p3 = 0.009, k_p4 = 0.0618
  - k_i = 0.0079, U_ii = 1.0
  - V_m0 = 2.5, V_mX = 0.047, K_m0 = 225.59
  - p_2U = 0.0331, part = 0.2, K = 2.3
  - alpha = 0.05, beta = 0.11, gamma = 0.5
  - k_e1 = 5.0e-4, k_e2 = 339.0, D = 78000.0

Initial conditions example (x0) used in `run_simulation` (first entries shown):
- Gp = 178.0
- Gt = 135.0
- Il = 4.5
- Ip = 1.25
- Qsto1 = 78000.0
- ... (rest defined in file)

Reduced-order model (ROM) parameters and inputs:
- The ROM uses physiologically-motivated parameters (S_G, k_cl, k_I, k_X, alpha, beta, theta, G_b, I_b).
- Ra (rate of appearance) is generated by a gamma-distribution shaped gastric operator used in `gastric_operator(t, D, shape, scale)`.


GenSSI identifiability setup:
- `dallaman_genssi.m` contains symbolic declarations of the 12 states and initial conditions and establishes the model structure expected by GenSSI to run structural identifiability analysis (Lie-derivative depth Nder set to 4 in the file). Use MATLAB + GenSSI toolbox to run identifiability experiments

## Contribution guidelines

- Open issues for major changes or experiments you plan to add.
- For code or notebook contributions:
  - Fork the repository and open a PR with a clear description of the scientific goal.
  - Keep notebook outputs minimal in commits (or clear outputs).
  - Add unit tests for helper functions in `src/` (or `models/`) wherever possible.
- Use semantic commit messages and link PRs to issues describing experiments or results.
