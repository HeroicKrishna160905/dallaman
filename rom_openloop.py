"""
rom_openloop.py

Reduced-Order Glucose–Insulin Model (ROM)
Open-loop simulation

Run with:
    python rom_openloop.py

Dependencies:
- numpy
- scipy
- matplotlib
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.stats import gamma
import matplotlib.pyplot as plt

# ------------------------------------------------------------
# Gastric Appearance Operator (Input)
# ------------------------------------------------------------

def gastric_operator(t, D, shape=3.0, scale=20.0):
    """
    Gamma-distributed gastric glucose appearance

    Parameters
    ----------
    t : float
        Time (min)
    D : float
        Meal size (mg or scaled units)
    """
    return D * gamma.pdf(t, a=shape, scale=scale)


# ------------------------------------------------------------
# ROM ODE Definition
# ------------------------------------------------------------

def rom_ode(t, x, p, Ra_func, u_I):
    """
    States
    ------
    x[0] = G : Plasma glucose
    x[1] = I : Plasma insulin
    x[2] = X : Insulin action

    Parameters
    ----------
    p : dict
        Model parameters
    """

    # Unpack states
    G, I, X = x

    # Unpack parameters
    S_G   = p['S_G']
    k_cl  = p['k_cl']
    k_I   = p['k_I']
    k_X   = p['k_X']
    alpha = p['alpha']
    beta  = p['beta']
    theta = p['theta']
    G_b   = p['G_b']
    I_b   = p['I_b']

    # Inputs
    Ra = Ra_func(t)
    u  = u_I(t)

    # Endogenous insulin secretion
    S_endo = beta * max(G - theta, 0.0)

    # ODEs
    dG = -(S_G + k_cl) * (G - G_b) - X * (G - G_b) + Ra
    dI = -k_I * (I - I_b) + S_endo + u
    dX = -k_X * X + alpha * (I - I_b)

    return np.array([dG, dI, dX])

# ------------------------------------------------------------
# Open-loop Simulation
# ------------------------------------------------------------

def run_simulation(t_span=(0, 300), dt=0.1):
    """
    Run open-loop ROM simulation
    """

    # Nominal parameters
    p = dict(
        S_G=0.01,
        k_cl=0.01,
        k_I=0.25,
        k_X=0.05,
        alpha=0.002,
        beta=0.03,
        theta=110.0,
        G_b=95.0,
        I_b=3.6
    )

    # Initial state
    x0 = np.array([
        p['G_b'],   # G(0)
        p['I_b'],   # I(0)
        0.0         # X(0)
    ])

    # Inputs
    Ra = lambda t: gastric_operator(t, D=120.0)
    uI = lambda t: 4.0 if 60 <= t <= 65 else 0.0

    # Time grid
    t_eval = np.arange(t_span[0], t_span[1] + dt, dt)

    # Integrate
    sol = solve_ivp(
        lambda t, x: rom_ode(t, x, p, Ra, uI),
        t_span,
        x0,
        t_eval=t_eval,
        method='RK45',
        atol=1e-8,
        rtol=1e-8
    )

    return sol, p, x0

# ------------------------------------------------------------
# Plotting
# ------------------------------------------------------------

def plot_solution(sol):
    t = sol.t
    G, I, X = sol.y

    plt.figure(figsize=(10,6))

    plt.subplot(3,1,1)
    plt.plot(t, G)
    plt.ylabel("Glucose (mg/dL)")
    plt.grid(True)

    plt.subplot(3,1,2)
    plt.plot(t, I)
    plt.ylabel("Insulin (µU/mL)")
    plt.grid(True)

    plt.subplot(3,1,3)
    plt.plot(t, X)
    plt.ylabel("Insulin Action")
    plt.xlabel("Time (min)")
    plt.grid(True)

    plt.tight_layout()
    plt.show()

# ------------------------------------------------------------
# Main
# ------------------------------------------------------------

if __name__ == "__main__":
    sol, p, x0 = run_simulation((0, 300), dt=0.5)
    plot_solution(sol)
