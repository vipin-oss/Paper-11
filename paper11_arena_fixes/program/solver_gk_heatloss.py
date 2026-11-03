"""
solver_gk_heatloss.py - 1D Guyer-Krumhansl Conduction Solver with Robin Boundary Heat Loss
Part of Paper 11 (Q1-SCI Revised Research Package).

Governing Equations (dimensionless):
  tau_Delta * dT/dt + dq/dx = 0
  tau_q * dq/dt + q + tau_Delta * dT/dx - kappa^2 * d2q/dx2 = 0

Boundary Conditions:
  Front (x = 0): q(0, t) = q_pulse(t) - tau_Delta * Bi0 * T(0, t)
  Rear  (x = 1): q(1, t) = tau_Delta * BiL * T(1, t)

Spatial Discretization:
  Staggered grid: T defined at cell centers (x_{i+1/2}, i = 0 ... Nx-1)
                 q defined at cell interfaces (x_i, i = 0 ... Nx)
  Second-order boundary temperature extrapolation:
    T(0, t) = 1.5 * T_0 - 0.5 * T_1
    T(1, t) = 1.5 * T_{Nx-1} - 0.5 * T_{Nx-2}

Time Integration:
  scipy.integrate.solve_ivp with stiff BDF / Radau implicit integrator.
"""

import numpy as np
from scipy.integrate import solve_ivp

def pulse_flux(t, tau_Delta=0.04, pulse_type="raised_cosine"):
    """
    Transient pulse heat flux at front boundary.
    - raised_cosine: 1 - cos(2*pi*t/tau_Delta) for 0 <= t <= tau_Delta
    - square: 1.0 for 0 <= t <= tau_Delta
    """
    if 0.0 < t <= tau_Delta:
        if pulse_type == "raised_cosine":
            return 1.0 - np.cos(2.0 * np.pi * t / tau_Delta)
        elif pulse_type == "square":
            return 1.0
        else:
            raise ValueError(f"Unknown pulse_type: {pulse_type}")
    return 0.0

def solve_gk_heatloss(Nx=200, 
                       tau_q=0.02, 
                       kappa2=0.02, 
                       Bi0=0.0, 
                       BiL=None, 
                       tau_Delta=0.04, 
                       pulse_type="raised_cosine",
                       t_eval=None, 
                       rtol=1e-10, 
                       atol=1e-12,
                       method="BDF"):
    """
    Solves the 1D Guyer-Krumhansl system with Robin convective/radiative cooling.

    Parameters:
      Nx: Number of spatial cells (default: 200)
      tau_q: Thermal relaxation time parameter
      kappa2: Nonlocal mean-free-path parameter
      Bi0: Front-face Biot number (default: 0.0)
      BiL: Rear-face Biot number (default: Bi0 if None)
      tau_Delta: Dimensionless pulse width (default: 0.04)
      pulse_type: "raised_cosine" or "square"
      t_eval: 1D array of evaluation time points
      rtol, atol: Integrator tolerances
      method: "BDF" or "Radau"

    Returns:
      t_eval: Evaluation times
      T_rear: Rear-face temperature history T(1, t)
      T_front: Front-face temperature history T(0, t)
      full_sol: solve_ivp solution object
    """
    if BiL is None:
        BiL = Bi0
    if t_eval is None:
        t_eval = np.linspace(0.01, 1.0, 100)

    dx = 1.0 / Nx

    def rhs(t, y):
        T = y[:Nx]
        q_int = y[Nx:]  # q at interfaces 1 through Nx-1

        # Second-order boundary temperature extrapolation
        T_front = 1.5 * T[0] - 0.5 * T[1]
        T_rear = 1.5 * T[Nx - 1] - 0.5 * T[Nx - 2]

        # Robin boundary conditions for flux
        q_front = pulse_flux(t, tau_Delta, pulse_type) - tau_Delta * Bi0 * T_front
        q_rear = tau_Delta * BiL * T_rear

        # Assemble full flux vector
        q = np.empty(Nx + 1)
        q[0] = q_front
        q[1:Nx] = q_int
        q[Nx] = q_rear

        # Energy conservation: dT/dt = -1/(tau_Delta * dx) * (q_{i+1} - q_i)
        dTdt = -(q[1:] - q[:-1]) / (tau_Delta * dx)

        # Temperature gradient at interior interfaces: (T_i - T_{i-1}) / dx
        dTdx = (T[1:] - T[:-1]) / dx

        # Flux diffusion: d2q/dx2 at interior interfaces
        d2qdx2 = (q[2:] - 2.0 * q[1:Nx] + q[:-2]) / (dx**2)

        # GK constitutive rate equation
        dqdt = (-q_int - tau_Delta * dTdx + kappa2 * d2qdx2) / tau_q

        return np.concatenate([dTdt, dqdt])

    y0 = np.zeros(2 * Nx - 1, dtype=float)
    t_span = (0.0, float(t_eval[-1]))

    sol = solve_ivp(rhs, t_span, y0, method=method, t_eval=t_eval, rtol=rtol, atol=atol)

    if not sol.success:
        raise RuntimeError(f"Solver failed: {sol.message}")

    T_grid = sol.y[:Nx, :]
    T_rear = 1.5 * T_grid[Nx - 1, :] - 0.5 * T_grid[Nx - 2, :]
    T_front = 1.5 * T_grid[0, :] - 0.5 * T_grid[1, :]

    return sol.t, T_rear, T_front, sol

if __name__ == "__main__":
    t, Tr, Tf, _ = solve_gk_heatloss(Nx=100, tau_q=0.02, kappa2=0.02, Bi0=0.1)
    print(f"Test run completed successfully. Rear peak: {np.max(Tr):.4f} at t = {t[np.argmax(Tr)]:.3f}")
