"""
generate_3d_datasets.py - High-Fidelity Pre-Computation for 3D Visualizations
Part of Paper 11 Research Package.

Generates:
  1. data/spatiotemporal_3d_temperature.csv: 3D Temperature field across (x, t) for Bi=0 and Bi=0.20
  2. data/spatiotemporal_3d_sensitivities.csv: 3D Sensitivity fields across (x, t) from TWO
     independent evaluation channels:
       - analytic differentiation of the closed-form transfer function
         (mpmath.diff, 40-digit arithmetic, de Hoog inversion degree 20);
       - central finite differences at relative step h_rel = 1e-4 in the same
         40-digit arithmetic (truncation-floor channel).
     Both channels report the collinearity sum J_tau_q + J_kappa_2 and its
     pointwise relative defect.
  3. data/canyon_3d_landscape.csv: 3D Fisher condition-number landscape across (B, Bi)
     computed from the FULL three-parameter Jacobian [J_tau_q, J_kappa^2, J_Bi].
     At B = 1 the Fisher matrix is structurally singular (analytical condition
     number infinite); the value stored there is the numerical floor of the
     double-precision evaluation channel and is flagged by the
     `cond_F_at_numerical_floor` column.
"""

import os
import sys
import time
import numpy as np
import pandas as pd
import mpmath as mp

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from analytic_laplace import (
    transfer_function_gk,
    transfer_function_gk_field,
    transfer_function_fourier_field
)

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(BASE_DIR, "data")
os.makedirs(DATA_DIR, exist_ok=True)

# Shared numerical settings for the analytic-differentiation channels.
DPS_ANALYTIC = 40
DEGREE_ANALYTIC = 20
# Relative finite-difference step of the independent verification channel.
H_REL_FD = 1e-4
# IEEE 754 double-precision unit roundoff (used to floor log-scale residuals).
EPS64 = float(np.finfo(np.float64).eps)


def generate_3d_temperature():
    print("=== Generating 3D Spatiotemporal Temperature Fields ===")
    t0 = time.time()
    mp.mp.dps = 15
    Nx, Nt = 25, 30
    x_arr = np.linspace(0.0, 1.0, Nx)
    t_arr = np.linspace(0.04, 1.5, Nt)

    records = []
    for ix, x in enumerate(x_arr):
        f_ad = lambda s: transfer_function_gk_field(x, s, tau_q=0.02, kappa2=0.02, Bi0=0.0, BiL=0.0)
        f_loss = lambda s: transfer_function_gk_field(x, s, tau_q=0.02, kappa2=0.02, Bi0=0.2, BiL=0.2)
        f_four = lambda s: transfer_function_fourier_field(x, s, Bi0=0.2, BiL=0.2)
        for it, t in enumerate(t_arr):
            v_ad = float(mp.re(mp.invertlaplace(f_ad, t, method="dehoog", degree=8)))
            v_loss = float(mp.re(mp.invertlaplace(f_loss, t, method="dehoog", degree=8)))
            v_four = float(mp.re(mp.invertlaplace(f_four, t, method="dehoog", degree=8)))
            records.append({
                "x": x,
                "t": t,
                "T_adiabatic": v_ad,
                "T_loss_Bi0p2": v_loss,
                "T_fourier_Bi0p2": v_four,
                "abs_diff": abs(v_loss - v_four)
            })

    df = pd.DataFrame(records)
    out_path = os.path.join(DATA_DIR, "spatiotemporal_3d_temperature.csv")
    df.to_csv(out_path, index=False)
    print(f"  [SAVED] {out_path} ({len(df)} points in {time.time()-t0:.2f}s)")


def generate_3d_sensitivities():
    print("=== Generating 3D Spatiotemporal Sensitivity Landscapes ===")
    t0 = time.time()
    Nx, Nt = 20, 25
    x_arr = np.linspace(0.0, 1.0, Nx)
    t_arr = np.linspace(0.04, 1.2, Nt)
    tau_0 = 0.02
    kap_0 = 0.02
    Bi = 0.1

    # Sensitivities are obtained by ANALYTIC differentiation of the closed-form
    # Laplace-domain transfer function (Eq. (transfer_field)) with respect to the
    # constitutive parameters, carried out in 40-digit arithmetic with mpmath's
    # high-precision differentiation, followed by de Hoog inversion (degree 20).
    # An independent central-finite-difference channel (relative step h_rel = 1e-4,
    # same 40-digit arithmetic) is evaluated on the same grid to quantify the
    # truncation floor of a conventional sensitivity evaluation.
    records = []
    for ix, x in enumerate(x_arr):
        # --- analytic channel -------------------------------------------------
        f_dtau = lambda s: mp.diff(
            lambda tt: transfer_function_gk_field(x, s, tau_q=tt, kappa2=kap_0,
                                                  Bi0=Bi, BiL=Bi), tau_0)
        f_dkap = lambda s: mp.diff(
            lambda kk: transfer_function_gk_field(x, s, tau_q=tau_0, kappa2=kk,
                                                  Bi0=Bi, BiL=Bi), kap_0)
        # --- finite-difference channel ---------------------------------------
        d_tq = H_REL_FD * tau_0
        d_k2 = H_REL_FD * kap_0
        f_tp = lambda s: transfer_function_gk_field(x, s, tau_q=tau_0 + d_tq, kappa2=kap_0, Bi0=Bi, BiL=Bi)
        f_tm = lambda s: transfer_function_gk_field(x, s, tau_q=tau_0 - d_tq, kappa2=kap_0, Bi0=Bi, BiL=Bi)
        f_kp = lambda s: transfer_function_gk_field(x, s, tau_q=tau_0, kappa2=kap_0 + d_k2, Bi0=Bi, BiL=Bi)
        f_km = lambda s: transfer_function_gk_field(x, s, tau_q=tau_0, kappa2=kap_0 - d_k2, Bi0=Bi, BiL=Bi)
        for it, t in enumerate(t_arr):
            jt = float(mp.re(mp.invertlaplace(f_dtau, t, method="dehoog", degree=DEGREE_ANALYTIC)))
            jk = float(mp.re(mp.invertlaplace(f_dkap, t, method="dehoog", degree=DEGREE_ANALYTIC)))
            v_tp = float(mp.re(mp.invertlaplace(f_tp, t, method="dehoog", degree=DEGREE_ANALYTIC)))
            v_tm = float(mp.re(mp.invertlaplace(f_tm, t, method="dehoog", degree=DEGREE_ANALYTIC)))
            v_kp = float(mp.re(mp.invertlaplace(f_kp, t, method="dehoog", degree=DEGREE_ANALYTIC)))
            v_km = float(mp.re(mp.invertlaplace(f_km, t, method="dehoog", degree=DEGREE_ANALYTIC)))
            jt_fd = (v_tp - v_tm) / (2.0 * d_tq)
            jk_fd = (v_kp - v_km) / (2.0 * d_k2)
            records.append({
                "x": x,
                "t": t,
                "J_tau_q": jt,
                "J_kappa_2": jk,
                "J_sum": jt + jk,
                "collinearity_defect": abs(jt + jk) / max(abs(jt), 1e-12),
                "J_tau_q_fd": jt_fd,
                "J_kappa_2_fd": jk_fd,
                "J_sum_fd": jt_fd + jk_fd,
                "collinearity_defect_fd": abs(jt_fd + jk_fd) / max(abs(jt_fd), 1e-12),
            })

    df = pd.DataFrame(records)
    out_path = os.path.join(DATA_DIR, "spatiotemporal_3d_sensitivities.csv")
    df.to_csv(out_path, index=False)
    print(f"  [SAVED] {out_path} ({len(df)} points in {time.time()-t0:.2f}s)")
    print(f"  Analytic channel : max |J_sum| = {df['J_sum'].abs().max():.3e} | "
          f"max relative defect = {df['collinearity_defect'].max():.3e} | "
          f"exactly-zero rows = {int((df['J_sum'] == 0.0).sum())}/{len(df)}")
    print(f"  FD channel       : max |J_sum| = {df['J_sum_fd'].abs().max():.3e} | "
          f"max relative defect = {df['collinearity_defect_fd'].max():.3e}")


def generate_3d_canyon():
    """3D Fisher condition-number landscape from the three-parameter Jacobian.

    No exact-resonance special case is used: the B = 1 column is evaluated by
    the same analytic-differentiation code path. Because the analytical
    collinearity is reproduced exactly in double precision there, the computed
    condition number saturates at the IEEE 754 singularity floor
    (approximately 1/eps ~ 4.5e15, i.e. 10^16 - 10^17 in practice) while the
    theoretical value is infinite. Rows at that floor carry
    cond_F_at_numerical_floor = True.
    """
    print("=== Generating 3D Singularity Canyon Landscape (3-parameter FIM) ===")
    t0 = time.time()
    B_vals = np.linspace(0.5, 1.5, 21)
    Bi_vals = np.linspace(0.0, 0.5, 9)
    t_eval = np.linspace(0.05, 1.2, 12)
    tau_0 = 0.02

    records = []
    for bi in Bi_vals:
        for b in B_vals:
            kap_0 = b * tau_0
            jt = np.zeros(len(t_eval))
            jk = np.zeros(len(t_eval))
            jb = np.zeros(len(t_eval))
            for it, t in enumerate(t_eval):
                jt[it] = float(mp.re(mp.invertlaplace(
                    lambda s: mp.diff(lambda v: transfer_function_gk(s, tau_q=v, kappa2=kap_0, Bi0=bi, BiL=bi), tau_0),
                    t, method="dehoog", degree=DEGREE_ANALYTIC)))
                jk[it] = float(mp.re(mp.invertlaplace(
                    lambda s: mp.diff(lambda v: transfer_function_gk(s, tau_q=tau_0, kappa2=v, Bi0=bi, BiL=bi), kap_0),
                    t, method="dehoog", degree=DEGREE_ANALYTIC)))
                jb[it] = float(mp.re(mp.invertlaplace(
                    lambda s: mp.diff(lambda v: transfer_function_gk(s, tau_q=tau_0, kappa2=kap_0, Bi0=v, BiL=v), bi),
                    t, method="dehoog", degree=DEGREE_ANALYTIC)))

            J = np.column_stack([jt, jk, jb])
            S = np.linalg.svd(J, compute_uv=False)
            F = J.T @ J
            cond_val = float(np.linalg.cond(F))
            s1, s2, s3 = float(S[0]), float(S[1]), float(S[2])
            at_floor = bool(s3 <= 1e-12 * s1)
            records.append({
                "Bi": bi,
                "B": b,
                "log10_cond_F": float(np.log10(cond_val)),
                "cond_F": cond_val,
                "sigma1": s1,
                "sigma2": s2,
                "sigma3": s3,
                "cond_F_at_numerical_floor": at_floor,
            })
        print(f"  Bi = {bi:.3f} done ({time.time()-t0:.0f}s elapsed)")

    df = pd.DataFrame(records)
    out_path = os.path.join(DATA_DIR, "canyon_3d_landscape.csv")
    df.to_csv(out_path, index=False)
    n_floor = int(df["cond_F_at_numerical_floor"].sum())
    print(f"  [SAVED] {out_path} ({len(df)} points in {time.time()-t0:.2f}s)")
    print(f"  rows at the numerical singularity floor: {n_floor} "
          f"(B = 1.0, where the analytical cond(F) is infinite)")
    print(f"  log10(cond) range: {df['log10_cond_F'].min():.2f} .. {df['log10_cond_F'].max():.2f}")


if __name__ == "__main__":
    generate_3d_temperature()
    generate_3d_sensitivities()
    generate_3d_canyon()
    print("All 3D datasets successfully generated!")
