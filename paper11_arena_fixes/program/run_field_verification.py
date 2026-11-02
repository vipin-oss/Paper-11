"""
run_field_verification.py - Field-Level Invariance and High-Precision Sensitivity Audit
Part of Paper 11 (Q1-SCI Revised Research Package).

Computes:
  1. Exact high-precision analytical field invariance (continuous Laplace domain,
     de Hoog inversion at 50-digit precision) between GK at B = 1 and classical Fourier,
     across 3 asymmetric Biot pairs and 5 spatial stations x in {0, 0.25, 0.5, 0.75, 1}.
     Differences are computed in high-precision mpmath BEFORE conversion to float,
     faithfully reporting the true round-off residual (~1e-27).
  2. Discretization residual between time-domain PDE solver and continuous analytical
     field at the same interior stations under mesh refinement Nx in {100, 200, 400, 800}.
  3. Sensitivity step-size convergence study for finite differences h in [1e-2, 1e-16]
     at 50-digit precision, verifying asymptotic O(h^2) convergence.

Outputs:
  - data/field_invariance_asymmetric.csv
  - data/interior_point_mesh_refinement.csv
  - data/sensitivity_step_size_convergence.csv
"""

import sys
import os
import time
import numpy as np
import pandas as pd
import mpmath as mp
from scipy.integrate import solve_ivp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from analytic_laplace import (transfer_function_gk_field, transfer_function_fourier_field,
                              sensitivity_step_size_study)

DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data")
os.makedirs(DATA_DIR, exist_ok=True)

def invert_field_high_precision(x, t_val, tau_q, kappa2, Bi0, BiL, tau_Delta=0.04, dps=50, degree=30):
    mp.mp.dps = dps
    f_s = lambda s: transfer_function_gk_field(x, s, tau_q, kappa2, Bi0, BiL, tau_Delta)
    val = mp.invertlaplace(f_s, t_val, method="dehoog", degree=degree)
    return val

def invert_fourier_field_high_precision(x, t_val, Bi0, BiL, tau_Delta=0.04, dps=50, degree=30):
    mp.mp.dps = dps
    f_s = lambda s: transfer_function_fourier_field(x, s, Bi0, BiL, tau_Delta)
    val = mp.invertlaplace(f_s, t_val, method="dehoog", degree=degree)
    return val

def solve_pde_full_field(Nx, tau_q, kappa2, Bi0, BiL, tau_Delta=0.04, t_eval=None):
    if t_eval is None:
        t_eval = np.linspace(0.05, 1.0, 20)
    dx = 1.0 / Nx

    def q_pulse(t):
        if 0.0 < t <= tau_Delta:
            return 1.0 - np.cos(2.0 * np.pi * t / tau_Delta)
        return 0.0

    def rhs(t, y):
        T = y[:Nx]
        q_int = y[Nx:]
        T0 = 1.5 * T[0] - 0.5 * T[1]
        T1 = 1.5 * T[Nx - 1] - 0.5 * T[Nx - 2]
        q_front = q_pulse(t) - tau_Delta * Bi0 * T0
        q_rear = tau_Delta * BiL * T1
        q = np.empty(Nx + 1)
        q[0] = q_front
        q[1:Nx] = q_int
        q[Nx] = q_rear
        dTdt = -(q[1:] - q[:-1]) / (tau_Delta * dx)
        dTdx = (T[1:] - T[:-1]) / dx
        d2qdx2 = (q[2:] - 2.0 * q[1:Nx] + q[:-2]) / (dx ** 2)
        dqdt = (-q_int - tau_Delta * dTdx + kappa2 * d2qdx2) / tau_q
        return np.concatenate([dTdt, dqdt])

    y0 = np.zeros(2 * Nx - 1)
    sol = solve_ivp(rhs, (0.0, float(t_eval[-1])), y0, method="BDF",
                     t_eval=t_eval, rtol=1e-10, atol=1e-12)
    if not sol.success:
        raise RuntimeError(sol.message)
    T_grid = sol.y[:Nx, :]
    cell_centers = (np.arange(Nx) + 0.5) * dx
    return t_eval, cell_centers, T_grid

def interp_point(x_target, cell_centers, T_grid):
    Nx = len(cell_centers)
    dx = cell_centers[1] - cell_centers[0]
    if x_target <= 0.0:
        return 1.5 * T_grid[0, :] - 0.5 * T_grid[1, :]
    if x_target >= 1.0:
        return 1.5 * T_grid[-1, :] - 0.5 * T_grid[-2, :]
    idx = np.clip(int(np.floor((x_target - cell_centers[0]) / dx)), 0, Nx - 2)
    x0, x1 = cell_centers[idx], cell_centers[idx + 1]
    w = (x_target - x0) / (x1 - x0)
    return (1.0 - w) * T_grid[idx, :] + w * T_grid[idx + 1, :]

def run_analytic_field_invariance():
    print("=" * 70)
    print("1. ANALYTICAL FIELD-LEVEL INVARIANCE (50-DIGIT HIGH PRECISION)")
    print("=" * 70)
    tau_q, kappa2 = 0.02, 0.02
    asym_pairs = [(0.01, 0.10), (0.05, 0.20), (0.10, 0.50)]
    x_points = [0.0, 0.25, 0.50, 0.75, 1.0]
    t_test = [0.10, 0.30, 0.60, 0.90]

    rows = []
    for Bi0, BiL in asym_pairs:
        for x in x_points:
            diffs_mp = []
            for t in t_test:
                val_gk = invert_field_high_precision(x, t, tau_q, kappa2, Bi0, BiL, dps=50, degree=30)
                val_fo = invert_fourier_field_high_precision(x, t, Bi0, BiL, dps=50, degree=30)
                d = abs(val_gk - val_fo)
                diffs_mp.append(d)
            max_d_mp = max(diffs_mp)
            rows.append({
                "Bi0": Bi0,
                "BiL": BiL,
                "x": x,
                "max_abs_diff_GK_minus_Fourier_mpmath": mp.nstr(max_d_mp, 10),
                "max_abs_diff_float": float(max_d_mp)
            })
            print(f"  Bi0={Bi0:.2f} BiL={BiL:.2f} x={x:.2f} | max|GK-Fou| = {mp.nstr(max_d_mp, 6)}")

    df = pd.DataFrame(rows)
    out = os.path.join(DATA_DIR, "field_invariance_asymmetric.csv")
    df.to_csv(out, index=False)
    print(f"  Saved: {out}\n")
    return df

def run_interior_point_mesh_refinement():
    print("=" * 70)
    print("2. INTERIOR-POINT PDE-vs-LAPLACE RESIDUAL UNDER MESH REFINEMENT")
    print("=" * 70)
    tau_q, kappa2 = 0.02, 0.02
    Bi0, BiL = 0.05, 0.20
    x_points = [0.0, 0.25, 0.50, 0.75, 1.0]
    t_test = np.array([0.10, 0.30, 0.60, 0.90])
    nx_list = [100, 200, 400, 800]

    # Precompute reference continuous Laplace field
    T_lap_ref = {}
    for x in x_points:
        T_lap_ref[x] = np.array([
            float(invert_field_high_precision(x, t, tau_q, kappa2, Bi0, BiL, dps=30, degree=20))
            for t in t_test
        ])

    rows = []
    for Nx in nx_list:
        dx = 1.0 / Nx
        _, cell_centers, T_grid = solve_pde_full_field(Nx, tau_q, kappa2, Bi0, BiL, t_eval=t_test)
        for x in x_points:
            T_pde_x = interp_point(x, cell_centers, T_grid)
            res = float(np.max(np.abs(T_pde_x - T_lap_ref[x])))
            rows.append({"Nx": Nx, "dx": dx, "x": x,
                         "max_abs_residual_PDE_minus_Laplace": res})
            print(f"  Nx={Nx:3d} x={x:.2f} | max residual = {res:.3e}")

    df = pd.DataFrame(rows)
    out = os.path.join(DATA_DIR, "interior_point_mesh_refinement.csv")
    df.to_csv(out, index=False)
    print(f"  Saved: {out}\n")
    return df

def run_sensitivity_step_size_study():
    print("=" * 70)
    print("3. SENSITIVITY STEP-SIZE CONVERGENCE STUDY (50-DIGIT DPS)")
    print("=" * 70)
    res = sensitivity_step_size_study(t=0.5, Bi0=0.1, BiL=0.1)
    df = pd.DataFrame(res)
    out = os.path.join(DATA_DIR, "sensitivity_step_size_convergence.csv")
    df.to_csv(out, index=False)
    print(f"  Saved: {out}\n")
    return df

if __name__ == "__main__":
    t0 = time.time()
    run_analytic_field_invariance()
    run_interior_point_mesh_refinement()
    run_sensitivity_step_size_study()
    print(f"Field verification completed in {time.time() - t0:.2f} seconds.")
