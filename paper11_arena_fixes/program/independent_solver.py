"""
independent_solver.py - Genuinely Independent Benchmark and Numerical Solvers
Part of Paper 11 (Q1-SCI Revised Research Package).

This module provides two completely independent cross-checks that do NOT
duplicate the staggered-grid finite-difference PDE solver:

1. Exact Cowan Fourier Eigenfunction Series Solver:
   Solves the classical Fourier Robin heat-loss problem in the time domain
   using 50 transcendental eigenvalue roots found by Brent's method.
   Involves zero spatial discretization and zero numerical Laplace inversion.

2. Modified Talbot Contour Laplace Inversion Solver:
   Evaluates the closed-form transfer function using a complex-plane contour
   integral along a deformed Talbot parabola. Entirely distinct algorithm
   from the continued-fraction/Padé-based de Hoog method.

3. Independent Finite Volume Method (FVM) Solver:
   Discretizes the 1D domain using cell-averaged conservative finite volumes,
   with harmonic-mean/central flux reconstruction, providing an independent
   spatial representation.
"""

import numpy as np
import mpmath as mp
from scipy.optimize import root_scalar
from scipy.integrate import quad, solve_ivp

def pulse_transform(s, tau_Delta=0.04):
    omega = 2.0 * np.pi / tau_Delta
    s = mp.mpc(s)
    return (omega**2 / (s * (s**2 + omega**2))) * (1.0 - mp.exp(-s * tau_Delta))

def tf_gk_laplace(s, tau_q, kappa2, Bi0, BiL, tau_Delta=0.04):
    s = mp.mpc(s)
    m = mp.sqrt(s * (1.0 + tau_q * s) / (1.0 + kappa2 * s))
    mu = (1.0 + kappa2 * s) / (1.0 + tau_q * s)
    m_mu = m * mu
    Delta = (m_mu**2 + Bi0 * BiL) * mp.sinh(m) + m_mu * (Bi0 + BiL) * mp.cosh(m)
    q_s = pulse_transform(s, tau_Delta)
    return (m_mu / (tau_Delta * Delta)) * q_s

def solve_talbot_laplace(t_eval, tau_q=0.02, kappa2=0.02, Bi0=0.1, BiL=0.1, tau_Delta=0.04, dps=30):
    """
    Independent numerical Laplace inversion using the modified Talbot contour method.
    """
    mp.mp.dps = dps
    t_eval = np.atleast_1d(t_eval)
    T_out = np.zeros(len(t_eval), dtype=float)
    f_s = lambda s: tf_gk_laplace(s, tau_q, kappa2, Bi0, BiL, tau_Delta)
    for i, t in enumerate(t_eval):
        val = mp.invertlaplace(f_s, t, method="talbot")
        T_out[i] = float(mp.re(val))
    return T_out

def solve_cowan_eigenfunctions(t_eval, Bi0=0.1, BiL=0.1, tau_Delta=0.04, x=1.0, n_modes=50):
    """
    Exact time-domain Fourier series for Robin boundary losses (Cowan 1963).
    Transcendental roots computed via Brent's method.
    """
    t_eval = np.atleast_1d(t_eval)
    omega = 2.0 * np.pi / tau_Delta

    f_smooth = lambda m: np.sin(m) * (m**2 - Bi0 * BiL) - (Bi0 + BiL) * m * np.cos(m)
    roots = []
    for n in range(1, n_modes + 1):
        try:
            bracket_a = (n - 1) * np.pi + 1e-6
            bracket_b = n * np.pi - 1e-6
            r = root_scalar(f_smooth, bracket=[bracket_a, bracket_b], method="brentq").root
            if r > 1e-6:
                roots.append(r)
        except Exception:
            pass

    T_out = np.zeros(len(t_eval), dtype=float)
    for mu in roots:
        X_0 = mu
        X_x = mu * np.cos(mu * x) + Bi0 * np.sin(mu * x)
        N, _ = quad(lambda y: (mu * np.cos(mu * y) + Bi0 * np.sin(mu * y))**2, 0, 1)

        for i, t in enumerate(t_eval):
            if t <= tau_Delta:
                term1 = (mu**4 + omega**2) - mu**4 * np.cos(omega * t) - mu**2 * omega * np.sin(omega * t) - omega**2 * np.exp(-mu**2 * t)
                I_t = term1 / (mu**2 * (mu**4 + omega**2))
            else:
                I_t = (omega**2 * (1.0 - np.exp(-mu**2 * tau_Delta)) / (mu**2 * (mu**4 + omega**2))) * np.exp(-mu**2 * (t - tau_Delta))
            T_out[i] += (X_0 * X_x / (N * tau_Delta)) * I_t

    return T_out

def solve_fvm_heatloss(Ncells=100, tau_q=0.02, kappa2=0.02, Bi0=0.1, BiL=0.1, tau_Delta=0.04, t_eval=None):
    """
    Genuinely independent Finite Volume Method (FVM) solver.
    Uses conservative cell averages for T with face fluxes, completely distinct
    implementation from solver_gk_heatloss.py.
    """
    if t_eval is None:
        t_eval = np.linspace(0.01, 1.0, 50)
    dx = 1.0 / Ncells

    def q_pulse(t):
        if 0.0 < t <= tau_Delta:
            return 1.0 - np.cos(2.0 * np.pi * t / tau_Delta)
        return 0.0

    # Unknowns: T in cells 0...Ncells-1, q at interior faces 1...Ncells-1
    def rhs(t, y):
        T = y[:Ncells]
        q_int = y[Ncells:]

        # Face temperatures by 2nd-order interpolation/extrapolation
        T_face0 = (3.0 * T[0] - T[1]) * 0.5
        T_faceL = (3.0 * T[-1] - T[-2]) * 0.5

        q_face0 = q_pulse(t) - tau_Delta * Bi0 * T_face0
        q_faceL = tau_Delta * BiL * T_faceL

        q_all = np.empty(Ncells + 1)
        q_all[0] = q_face0
        q_all[1:Ncells] = q_int
        q_all[Ncells] = q_faceL

        # Cell balance: d(T_i)/dt = - (q_{i+1/2} - q_{i-1/2}) / (tau_Delta * dx)
        dTdt = -(q_all[1:] - q_all[:-1]) / (tau_Delta * dx)

        # Internal flux rate:
        # dT/dx at face i: (T_i - T_{i-1}) / dx
        dTdx_face = (T[1:] - T[:-1]) / dx
        # d2q/dx2 at face i: (q_{i+1} - 2*q_i + q_{i-1}) / dx^2
        d2qdx2_face = (q_all[2:] - 2.0 * q_int + q_all[:-2]) / (dx**2)

        dqdt = (-q_int - tau_Delta * dTdx_face + kappa2 * d2qdx2_face) / tau_q
        return np.concatenate([dTdt, dqdt])

    y0 = np.zeros(2 * Ncells - 1)
    sol = solve_ivp(rhs, (0.0, float(t_eval[-1])), y0, method="BDF", t_eval=t_eval, rtol=1e-8, atol=1e-10)
    T_rear = (3.0 * sol.y[Ncells - 1, :] - sol.y[Ncells - 2, :]) * 0.5
    return sol.t, T_rear

if __name__ == "__main__":
    t_test = np.array([0.1, 0.2, 0.5, 0.8, 1.0])
    print("Cross-checking Independent Solvers for Bi = 0.1 at B = 1:")
    T_talbot = solve_talbot_laplace(t_test, Bi0=0.1, BiL=0.1)
    T_cowan = solve_cowan_eigenfunctions(t_test, Bi0=0.1, BiL=0.1)
    _, T_fvm = solve_fvm_heatloss(Ncells=200, Bi0=0.1, BiL=0.1, t_eval=t_test)
    for t, tt, tc, tf in zip(t_test, T_talbot, T_cowan, T_fvm):
        print(f"t={t:4.2f} | Talbot={tt:10.7f} | CowanSer={tc:10.7f} | FVM={tf:10.7f} | |Tal-Cow|={abs(tt-tc):8.2e}")
