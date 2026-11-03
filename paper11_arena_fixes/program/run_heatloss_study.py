"""
run_heatloss_study.py - Systematic Boundary Heat Loss Parameter Sweeps for Paper 11
Part of Paper 11 (Q1-SCI Revised Research Package).

Computes:
  1. Systematic Biot number sweep at Fourier resonance (B = 1.0):
     Bi in [0.0, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5]
     Evaluates:
     - J_tau_q, J_kappa^2, J_Bi
     - Collinearity correlation rho, defect |1 + rho|, residual norm R_J
     - SVD spectrum: sigma_1 (bulk diffusion), sigma_2 (convective cooling tail),
       sigma_3 (finite-difference numerical artifact floor of the theoretical structural null vector)
     - Observable condition ratio sigma_1 / sigma_2
     - Condition numbers for 2x2 and 3x3 Fisher information matrices (under i.i.d. noise W = sigma^-2 * I).
  2. Time-resolved sensitivity profiles for representative Biot numbers (Bi in {0.0, 0.05, 0.2}).

Generates:
  - data/heatloss_parameter_sweep.csv
  - data/sensitivity_profiles.csv
"""

import sys
import os
import time
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from solver_gk_heatloss import solve_gk_heatloss

DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data")
os.makedirs(DATA_DIR, exist_ok=True)

def run_resonance_heatloss_sweep():
    print("=" * 70)
    print("1. RUNNING FOURIER-RESONANCE (B = 1.0) HEAT-LOSS SWEEP")
    print("=" * 70)

    bi_list = [0.0, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5]
    t_eval = np.linspace(0.01, 1.0, 100)
    tau_q_0 = 0.02
    kappa2_0 = 0.02
    alpha_0 = 1.0
    Nx = 100
    h_rel = 1e-4

    results = []

    for Bi in bi_list:
        print(f"  Computing Bi = {Bi:5.3f} ...", end=" ", flush=True)
        t_start = time.perf_counter()

        d_tq = h_rel * tau_q_0
        d_k2 = h_rel * kappa2_0
        d_bi = h_rel * Bi if Bi > 0 else 1e-5

        # Base solution
        _, T_base, _, _ = solve_gk_heatloss(Nx=Nx, tau_q=tau_q_0, kappa2=kappa2_0, Bi0=Bi, BiL=Bi, t_eval=t_eval)

        # Perturbation in tau_q
        _, Tp_tq, _, _ = solve_gk_heatloss(Nx=Nx, tau_q=tau_q_0 + d_tq, kappa2=kappa2_0, Bi0=Bi, BiL=Bi, t_eval=t_eval)
        _, Tm_tq, _, _ = solve_gk_heatloss(Nx=Nx, tau_q=tau_q_0 - d_tq, kappa2=kappa2_0, Bi0=Bi, BiL=Bi, t_eval=t_eval)
        J_tq = (Tp_tq - Tm_tq) / (2.0 * d_tq)

        # Perturbation in kappa2
        _, Tp_k2, _, _ = solve_gk_heatloss(Nx=Nx, tau_q=tau_q_0, kappa2=kappa2_0 + d_k2, Bi0=Bi, BiL=Bi, t_eval=t_eval)
        _, Tm_k2, _, _ = solve_gk_heatloss(Nx=Nx, tau_q=tau_q_0, kappa2=kappa2_0 - d_k2, Bi0=Bi, BiL=Bi, t_eval=t_eval)
        J_k2 = (Tp_k2 - Tm_k2) / (2.0 * d_k2)

        # Perturbation in Bi
        _, Tp_bi, _, _ = solve_gk_heatloss(Nx=Nx, tau_q=tau_q_0, kappa2=kappa2_0, Bi0=Bi + d_bi, BiL=Bi + d_bi, t_eval=t_eval)
        _, Tm_bi, _, _ = solve_gk_heatloss(Nx=Nx, tau_q=tau_q_0, kappa2=kappa2_0, Bi0=max(0.0, Bi - d_bi), BiL=max(0.0, Bi - d_bi), t_eval=t_eval)
        denom_bi = (2.0 * d_bi) if Bi > 0 else d_bi
        J_bi = (Tp_bi - Tm_bi) / denom_bi

        # Collinearity metrics
        norm_tq = np.linalg.norm(J_tq)
        norm_k2 = np.linalg.norm(J_k2)
        rho_val = float(np.dot(J_tq, J_k2) / (norm_tq * norm_k2))
        one_plus_rho = float(1.0 + rho_val)
        res_norm = float(np.linalg.norm(J_tq + alpha_0 * J_k2) / norm_tq)

        # 2x2 Fisher matrix
        J2 = np.column_stack([J_tq, J_k2])
        s2 = np.linalg.svd(J2, compute_uv=False)
        cond_F2 = float((s2[0] / s2[1])**2)

        # 3x3 Fisher matrix
        J3 = np.column_stack([J_tq, J_k2, J_bi])
        s3 = np.linalg.svd(J3, compute_uv=False)
        cond_F3 = float((s3[0] / s3[2])**2)

        # Observable subspace conditioning
        ratio_sigma1_sigma2 = float(s3[0] / s3[1])

        elapsed = time.perf_counter() - t_start

        results.append({
            "Bi": Bi,
            "tau_q": tau_q_0,
            "kappa2": kappa2_0,
            "norm_J_tq": norm_tq,
            "norm_J_k2": norm_k2,
            "norm_J_bi": float(np.linalg.norm(J_bi)),
            "rho": rho_val,
            "one_plus_rho": one_plus_rho,
            "R_J": res_norm,
            "sigma1_2p": float(s2[0]),
            "sigma2_2p": float(s2[1]),
            "cond_F2": cond_F2,
            "sigma1_3p": float(s3[0]),
            "sigma2_3p": float(s3[1]),
            "sigma3_3p": float(s3[2]),
            "sigma1_over_sigma2": ratio_sigma1_sigma2,
            "cond_F3": cond_F3
        })
        print(f"Done in {elapsed:.1f}s | rho = {rho_val:.10f}, R_J = {res_norm:.2e}, sigma_3 = {s3[2]:.2e}")

    df = pd.DataFrame(results)
    out_path = os.path.join(DATA_DIR, "heatloss_parameter_sweep.csv")
    df.to_csv(out_path, index=False)
    print(f"  Saved heat loss parameter sweep to: {out_path}\n")
    return df

def run_sensitivity_profiles():
    print("=" * 70)
    print("2. GENERATING TIME-RESOLVED SENSITIVITY PROFILES")
    print("=" * 70)

    bi_samples = [0.0, 0.05, 0.2]
    t_eval = np.linspace(0.01, 1.0, 100)
    tau_q_0 = 0.02
    kappa2_0 = 0.02
    Nx = 100
    h_rel = 1e-4

    records = {"time": t_eval}

    for Bi in bi_samples:
        d_tq = h_rel * tau_q_0
        d_k2 = h_rel * kappa2_0
        d_bi = h_rel * Bi if Bi > 0 else 1e-5

        _, T_base, _, _ = solve_gk_heatloss(Nx=Nx, tau_q=tau_q_0, kappa2=kappa2_0, Bi0=Bi, BiL=Bi, t_eval=t_eval)
        _, Tp_tq, _, _ = solve_gk_heatloss(Nx=Nx, tau_q=tau_q_0 + d_tq, kappa2=kappa2_0, Bi0=Bi, BiL=Bi, t_eval=t_eval)
        _, Tm_tq, _, _ = solve_gk_heatloss(Nx=Nx, tau_q=tau_q_0 - d_tq, kappa2=kappa2_0, Bi0=Bi, BiL=Bi, t_eval=t_eval)
        J_tq = (Tp_tq - Tm_tq) / (2.0 * d_tq)

        _, Tp_k2, _, _ = solve_gk_heatloss(Nx=Nx, tau_q=tau_q_0, kappa2=kappa2_0 + d_k2, Bi0=Bi, BiL=Bi, t_eval=t_eval)
        _, Tm_k2, _, _ = solve_gk_heatloss(Nx=Nx, tau_q=tau_q_0, kappa2=kappa2_0 - d_k2, Bi0=Bi, BiL=Bi, t_eval=t_eval)
        J_k2 = (Tp_k2 - Tm_k2) / (2.0 * d_k2)

        _, Tp_bi, _, _ = solve_gk_heatloss(Nx=Nx, tau_q=tau_q_0, kappa2=kappa2_0, Bi0=Bi + d_bi, BiL=Bi + d_bi, t_eval=t_eval)
        _, Tm_bi, _, _ = solve_gk_heatloss(Nx=Nx, tau_q=tau_q_0, kappa2=kappa2_0, Bi0=max(0.0, Bi - d_bi), BiL=max(0.0, Bi - d_bi), t_eval=t_eval)
        denom_bi = (2.0 * d_bi) if Bi > 0 else d_bi
        J_bi = (Tp_bi - Tm_bi) / denom_bi

        tag = f"Bi_{Bi:04.2f}".replace(".", "p")
        records[f"T_{tag}"] = T_base
        records[f"J_tau_q_{tag}"] = J_tq
        records[f"J_kappa2_{tag}"] = J_k2
        records[f"J_Bi_{tag}"] = J_bi

    df = pd.DataFrame(records)
    out_path = os.path.join(DATA_DIR, "sensitivity_profiles.csv")
    df.to_csv(out_path, index=False)
    print(f"  Saved sensitivity profiles to: {out_path}\n")
    return df

if __name__ == "__main__":
    t0 = time.time()
    run_resonance_heatloss_sweep()
    run_sensitivity_profiles()
    print(f"Heat loss study completed in {time.time() - t0:.2f} seconds.")
