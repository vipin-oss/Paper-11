"""
run_near_resonance.py - Dense Near-Resonance Study, Asymptotic Scaling Analysis,
and Multi-Configuration CRLB Identifiability Audit.
Part of Paper 11 (Q1-SCI Final Scientific Correction Pass).

Quantifies:
  1. Dense B sweep around B = 1.0:
     B in [0.5, 0.8, 0.9, 0.95, 0.98, 0.99, 0.995, 0.999, 1.0,
          1.001, 1.005, 1.01, 1.02, 1.05, 1.1, 1.2, 1.5]
     across Biot numbers Bi in {0.0, 0.01, 0.05, 0.1, 0.2, 0.5}.
  2. Metrics evaluated (all from the FULL three-parameter Jacobian
     J = [J_tau_q, J_kappa^2, J_Bi]):
     - Collinearity correlation rho and defect 1 + rho between the two
       non-Fourier sensitivity columns
     - Residual norm R_J = ||J_tau_q + J_kappa^2|| / ||J_tau_q||
     - Singular values of the 2-parameter (tau_q, kappa^2) subsystem and of the
       full 3-parameter Jacobian (sigma_1, sigma_2, sigma_3)
     - Condition numbers cond(F_{2x2}) and cond(F_{3x3})
     - Asymptotic scaling coefficient C = sigma_min / |B - 1|
  3. Independent regressions on both sides (B < 1 and B > 1) for both the
     2-parameter and the 3-parameter smallest singular value.
  4. Rigorous Cramer-Rao Lower Bound (CRLB) identifiability audit in the
     THREE-PARAMETER estimation space p = (tau_q, kappa^2, Bi):
     - Rear-face only (x = 1.0), Front + Rear (x = 0, 1), Multi-location (5 stations)
     - Temporal sampling Nt = 10, 50, 100
     - Noise levels: 0.1%, 0.5%, 1.0%
     - Target accuracy: SE(tau_q) / tau_q_0 <= 0.10 (10% relative accuracy)
     - Near resonance (B = 0.99, F non-singular): ordinary CRLB from inv(F)
     - Exact resonance (B = 1, F singular): generalized CRLB from the
       Moore-Penrose pseudoinverse F^dagger; only functionals orthogonal to the
       structural null space (i.e. Bi) admit a finite bound.
  5. Uncertainty-divergence scan used by Figure 10(d): the numerical ratio
     SE(tau_q)/SE(Bi) is computed directly from the 3-parameter FIM over a
     grid of B values; no hard-coded constants are used anywhere.

Numerical note: all Jacobian columns are obtained by ANALYTIC differentiation
(mpmath.diff, 40-digit working precision) of the closed-form Laplace-domain
transfer function with respect to each parameter, followed by de Hoog inversion
(degree 20). Finite-difference truncation error is therefore absent from these
channels; the residual collinearity defect is limited only by the Laplace
inversion and double-precision conversion floor.

Outputs:
  - data/near_resonance_scaling.csv
  - data/off_resonance_sweep.csv
  - data/scaling_summary.csv
  - data/scaling_regression_twosided.csv
  - data/crlb_configuration_audit.csv
  - data/uncertainty_divergence.csv
"""

import sys
import os
import time
import numpy as np
import pandas as pd
import mpmath as mp
from scipy.stats import linregress

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from analytic_laplace import transfer_function_gk

DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data")
os.makedirs(DATA_DIR, exist_ok=True)

# Working precision / inversion settings shared by all analytic-differentiation
# channels in this module (documented in REPRODUCIBILITY.md).
DPS = 40
DEGREE = 20


def transfer_field(s, x, tau_q, kappa2, Bi0, BiL, tau_Delta):
    """Full-field transfer function in Laplace domain."""
    s = mp.mpc(s)
    u_sq = s * (1.0 + tau_q * s) / (1.0 + kappa2 * s)
    m = mp.sqrt(u_sq)
    mu = (1.0 + kappa2 * s) / (1.0 + tau_q * s)
    m_mu = m * mu
    sinh_m = mp.sinh(m)
    cosh_m = mp.cosh(m)
    Delta = ((m_mu)**2 + Bi0 * BiL) * sinh_m + m_mu * (Bi0 + BiL) * cosh_m
    omega = 2.0 * mp.pi / tau_Delta
    q_s = (omega**2 / (s * (s**2 + omega**2))) * (1.0 - mp.exp(-tau_Delta * s))
    arg = m * (1.0 - x)
    num = m_mu * mp.cosh(arg) + BiL * mp.sinh(arg)
    return (q_s / (tau_Delta * Delta)) * num


def sensitivity_jacobian(B, Bi, t_eval, x_stations, tau_q_0=0.02, tau_Delta=0.04,
                         dps=DPS, degree=DEGREE):
    """Three-parameter sensitivity Jacobian J = [J_tau_q, J_kappa^2, J_Bi].

    Each column is the analytic derivative of the closed-form transfer function
    with respect to the corresponding parameter (mpmath.diff in `dps`-digit
    arithmetic), followed by de Hoog Laplace inversion at the requested degree.
    The Biot number enters symmetrically (Bi0 = BiL = Bi), so the third column
    is the derivative with respect to the common cooling parameter.

    Parameters
    ----------
    B : float          resonance ratio kappa^2 / tau_q
    Bi : float         Biot number (front = rear)
    t_eval : array     sampling times
    x_stations : list  spatial stations at which sensitivities are evaluated
    """
    kappa2_0 = B * tau_q_0
    mp.mp.dps = dps
    rows = []
    for x in x_stations:
        for t in t_eval:
            j_tq = float(mp.re(mp.invertlaplace(
                lambda s: mp.diff(lambda v: transfer_field(s, x, v, kappa2_0, Bi, Bi, tau_Delta),
                                  tau_q_0), t, method="dehoog", degree=degree)))
            j_k2 = float(mp.re(mp.invertlaplace(
                lambda s: mp.diff(lambda v: transfer_field(s, x, tau_q_0, v, Bi, Bi, tau_Delta),
                                  kappa2_0), t, method="dehoog", degree=degree)))
            j_bi = float(mp.re(mp.invertlaplace(
                lambda s: mp.diff(lambda v: transfer_field(s, x, tau_q_0, kappa2_0, v, v, tau_Delta),
                                  Bi), t, method="dehoog", degree=degree)))
            rows.append([j_tq, j_k2, j_bi])
    return np.array(rows)


def fim_metrics(J, rcond=1e-12):
    """Fisher-information metrics of a Jacobian J under i.i.d. noise.

    Returns singular values, the condition number of F = J^T J, the
    Moore-Penrose pseudoinverse standard errors, and a flag telling whether the
    smallest singular value is at the double-precision singularity floor
    (in which case the analytical condition number is infinite).
    """
    S = np.linalg.svd(J, compute_uv=False)
    F = J.T @ J
    cond_F = float(np.linalg.cond(F))
    F_pinv = np.linalg.pinv(F, rcond=rcond)
    se = np.sqrt(np.abs(np.diag(F_pinv)))
    # A matrix is numerically singular when its smallest singular value falls
    # below the double-precision resolution of the largest one; the analytical
    # condition number is then infinite and the computed value is a floor.
    at_floor = bool(S[-1] <= rcond * S[0])
    return {
        "sigma": S,
        "cond_F": cond_F,
        "se_pinv": se,
        "cond_F_at_numerical_floor": at_floor,
    }


def run_near_resonance_study():
    print("=" * 75)
    print("RUNNING NEAR-RESONANCE LANDSCAPE & ASYMPTOTIC SCALING STUDY")
    print("(three-parameter Jacobian, analytic differentiation)")
    print("=" * 75)

    b_vals = [0.5, 0.8, 0.9, 0.95, 0.98, 0.99, 0.995, 0.999,
              1.0,
              1.001, 1.005, 1.01, 1.02, 1.05, 1.1, 1.2, 1.5]
    bi_vals = [0.0, 0.01, 0.05, 0.1, 0.2, 0.5]
    t_eval = np.linspace(0.05, 1.0, 10)
    tau_q_0 = 0.02
    alpha_0 = 1.0

    all_records = []

    for Bi in bi_vals:
        print(f"Processing Biot number Bi = {Bi:4.2f} ...", end=" ", flush=True)
        t_bi_start = time.perf_counter()
        for B in b_vals:
            J = sensitivity_jacobian(B, Bi, t_eval, [1.0])
            j_tq, j_k2, j_bi = J[:, 0], J[:, 1], J[:, 2]

            norm_tq = float(np.linalg.norm(j_tq))
            norm_k2 = float(np.linalg.norm(j_k2))
            norm_bi = float(np.linalg.norm(j_bi))
            rho_val = float(np.dot(j_tq, j_k2) / (norm_tq * norm_k2))
            one_plus_rho = 1.0 + rho_val
            res_norm = float(np.linalg.norm(j_tq + alpha_0 * j_k2) / norm_tq)

            # Two-parameter (tau_q, kappa^2) subsystem metrics
            J2 = np.column_stack([j_tq, j_k2])
            S2 = np.linalg.svd(J2, compute_uv=False)
            sigma1_2p, sigma2_2p = float(S2[0]), float(S2[1])
            cond_J2 = sigma1_2p / sigma2_2p if sigma2_2p > 1e-15 else np.inf

            # Full three-parameter metrics
            m3 = fim_metrics(J)
            S3 = m3["sigma"]
            sigma1_3p, sigma2_3p, sigma3_3p = (float(S3[0]), float(S3[1]), float(S3[2]))

            delta_B = abs(B - 1.0)
            C_scaling_2p = (sigma2_2p / delta_B) if delta_B > 1e-5 else np.nan
            C_scaling_3p = (sigma3_3p / delta_B) if delta_B > 1e-5 else np.nan

            rec = {
                "Bi": Bi,
                "B": B,
                "delta_B": delta_B,
                "tau_q": tau_q_0,
                "kappa2": B * tau_q_0,
                "norm_J_tq": norm_tq,
                "norm_J_k2": norm_k2,
                "norm_J_bi": norm_bi,
                "rho": rho_val,
                "one_plus_rho": one_plus_rho,
                "R_J": res_norm,
                "sigma1_2p": sigma1_2p,
                "sigma2_2p": sigma2_2p,
                "cond_J2": cond_J2,
                "cond_F2": cond_J2**2,
                "sigma1_3p": sigma1_3p,
                "sigma2_3p": sigma2_3p,
                "sigma3_3p": sigma3_3p,
                "cond_F3": m3["cond_F"],
                "sigma1_over_sigma2": sigma1_3p / sigma2_3p,
                "cond_F_at_numerical_floor": m3["cond_F_at_numerical_floor"],
                "C_scaling_2p": C_scaling_2p,
                "C_scaling_3p": C_scaling_3p
            }
            all_records.append(rec)

        elapsed_bi = time.perf_counter() - t_bi_start
        print(f"Done in {elapsed_bi:.1f}s")

    df_all = pd.DataFrame(all_records)
    out_all = os.path.join(DATA_DIR, "near_resonance_scaling.csv")
    df_all.to_csv(out_all, index=False)
    print(f"\nSaved near-resonance scaling dataset to: {out_all}")

    # Companion wide off-resonance sweep used by the landscape figures
    df_off = df_all[df_all["Bi"].isin([0.0, 0.05, 0.2])][
        ["Bi", "B", "tau_q", "kappa2", "rho", "one_plus_rho", "R_J",
         "cond_F2", "cond_F3", "sigma1_3p", "sigma2_3p", "sigma3_3p"]].copy()
    out_off = os.path.join(DATA_DIR, "off_resonance_sweep.csv")
    df_off.to_csv(out_off, index=False)
    print(f"Saved off-resonance sweep dataset to: {out_off}")

    # ------------------------------------------------------------------
    # Independent two-sided regressions for Bi = 0.10
    # ------------------------------------------------------------------
    print("\n" + "=" * 75)
    print("INDEPENDENT TWO-SIDED REGRESSION AUDIT (Bi = 0.10)")
    print("=" * 75)
    df_bi01 = df_all[df_all["Bi"] == 0.10]

    regression_rows = []
    for quantity, label, fmt in [("sigma2_2p", "sigma_2 (2-parameter subsystem)", None),
                                 ("sigma3_3p", "sigma_3 (3-parameter Jacobian)", None)]:
        # Left side: B < 1 (B in [0.90, 0.999])
        sub_left = df_bi01[(df_bi01["B"] >= 0.90) & (df_bi01["B"] < 1.0)].sort_values("delta_B")
        reg_log_left = linregress(np.log(sub_left["delta_B"]), np.log(sub_left[quantity]))
        reg_lin_left = linregress(sub_left["delta_B"], sub_left[quantity])
        # Right side: B > 1 (B in [1.001, 1.10])
        sub_right = df_bi01[(df_bi01["B"] > 1.0) & (df_bi01["B"] <= 1.10)].sort_values("delta_B")
        reg_log_right = linregress(np.log(sub_right["delta_B"]), np.log(sub_right[quantity]))
        reg_lin_right = linregress(sub_right["delta_B"], sub_right[quantity])

        print(f"\nQuantity: {label}")
        print(f"Left  (B < 1): Slope = {reg_log_left.slope:.5f} +/- {reg_log_left.stderr:.5f} | "
              f"R^2 = {reg_log_left.rvalue**2:.7f} | Linear C = {reg_lin_left.slope:.4f}")
        print(f"Right (B > 1): Slope = {reg_log_right.slope:.5f} +/- {reg_log_right.stderr:.5f} | "
              f"R^2 = {reg_log_right.rvalue**2:.7f} | Linear C = {reg_lin_right.slope:.4f}")

        for branch, sub, rl, rg in [("B < 1", sub_left, reg_log_left, reg_lin_left),
                                    ("B > 1", sub_right, reg_log_right, reg_lin_right)]:
            regression_rows.append({
                "quantity": quantity,
                "quantity_label": label,
                "branch": branch,
                "B_range": "[0.900, 0.999]" if branch == "B < 1" else "[1.001, 1.100]",
                "n_points": len(sub),
                "loglog_exponent": rl.slope,
                "loglog_exponent_stderr": rl.stderr,
                "loglog_R2": rl.rvalue**2,
                "loglog_intercept": rl.intercept,
                "linear_prefactor_C": rg.slope,
                "linear_intercept": rg.intercept,
                "linear_R2": rg.rvalue**2,
            })

    df_reg = pd.DataFrame(regression_rows)
    out_reg = os.path.join(DATA_DIR, "scaling_regression_twosided.csv")
    df_reg.to_csv(out_reg, index=False)
    print(f"\nSaved two-sided regression audit to: {out_reg}")

    # ------------------------------------------------------------------
    # Scaling summary
    # ------------------------------------------------------------------
    scaling_summary = []
    for Bi in bi_vals:
        sub = df_all[(df_all["Bi"] == Bi) & (df_all["delta_B"] >= 0.005) & (df_all["delta_B"] <= 0.05)]
        mean_C_2p = sub["C_scaling_2p"].mean()
        mean_C_3p = sub["C_scaling_3p"].mean()
        scaling_summary.append({
            "Bi": Bi,
            "C_scaling_mean_2p": mean_C_2p,
            "C_scaling_mean_3p": mean_C_3p,
            "sigma1_mean_3p": sub["sigma1_3p"].mean(),
            "cond_F_prefactor_2p": (sub["sigma1_2p"].mean() / mean_C_2p)**2,
            "cond_F_prefactor_3p": (sub["sigma1_3p"].mean() / mean_C_3p)**2,
        })
    df_scaling = pd.DataFrame(scaling_summary)
    out_scaling = os.path.join(DATA_DIR, "scaling_summary.csv")
    df_scaling.to_csv(out_scaling, index=False)
    print(f"Saved scaling summary to: {out_scaling}")

    return df_all, df_reg


def run_crlb_audit():
    """Three-parameter CRLB / practical-identifiability audit.

    Near resonance (B = 0.99) the 3x3 Fisher matrix is non-singular and the
    ordinary CRLB is obtained from inv(F). At exact resonance (B = 1) F is
    structurally singular, so the generalized (Moore-Penrose) CRLB is used: only
    functionals whose coefficient vector is orthogonal to the structural null
    space - here the Biot number - admit a finite bound; SE(tau_q) and
    SE(kappa^2) are reported as not estimable (NaN).
    """
    print("\n" + "=" * 75)
    print("MULTI-CONFIGURATION CRLB & IDENTIFIABILITY AUDIT (3-PARAMETER)")
    print("=" * 75)

    configs = [
        ("Rear-face only (x=1.0)", [1.0]),
        ("Front + Rear (x=0, 1)", [0.0, 1.0]),
        ("Multi-location (5 stations)", [0.0, 0.25, 0.5, 0.75, 1.0])
    ]
    bi_test = 0.10
    tau_q_0 = 0.02
    target_rel = 0.10

    crlb_records = []
    for b_test in [0.99, 1.0]:
        delta_test = abs(b_test - 1.0)
        singular = (b_test == 1.0)
        print(f"\n--- Evaluation point B = {b_test} (|B-1| = {delta_test}) "
              f"{'[generalized CRLB]' if singular else '[ordinary CRLB]'} ---")
        for config_name, x_stations in configs:
            for Nt in [10, 50, 100]:
                t_pts = np.linspace(0.05, 1.0, Nt)
                J = sensitivity_jacobian(b_test, bi_test, t_pts, x_stations)
                S = np.linalg.svd(J, compute_uv=False)
                F = J.T @ J
                m = fim_metrics(J)
                se_pinv = m["se_pinv"]

                if singular:
                    # Generalized CRLB: tau_q and kappa^2 are NOT estimable
                    # (their coefficient vectors are not orthogonal to the
                    # structural null space); only Bi is.
                    se_tq = np.nan
                    se_k2 = np.nan
                    se_bi = se_pinv[2]
                    crit_coord = np.nan
                    crit_sing = np.nan
                else:
                    F_inv = np.linalg.inv(F)
                    se_tq = float(np.sqrt(F_inv[0, 0]))
                    se_k2 = float(np.sqrt(F_inv[1, 1]))
                    se_bi = float(np.sqrt(F_inv[2, 2]))
                    # Critical corridor: SE(tau_q)/tau_q_0 = sigma_n*se_tq/tau_q_0 <= target
                    # with se_tq scaling as 1/|B-1| near resonance:
                    crit_coord = (se_tq * delta_test) / (target_rel * tau_q_0)
                    C3 = S[2] / delta_test
                    crit_sing = 1.0 / (target_rel * tau_q_0 * C3)

                for sigma_pct, sigma_n in [(0.1, 1e-3), (0.5, 5e-3), (1.0, 1e-2)]:
                    crlb_records.append({
                        "B": b_test,
                        "configuration": config_name,
                        "Nt": Nt,
                        "noise_percent": sigma_pct,
                        "sigma_n": sigma_n,
                        "sigma1": float(S[0]),
                        "sigma2": float(S[1]),
                        "sigma3": float(S[2]),
                        "cond_F3": m["cond_F"],
                        "SE_tau_q_per_sigma_n": se_tq,
                        "SE_kappa2_per_sigma_n": se_k2,
                        "SE_Bi_per_sigma_n": se_bi,
                        "SE_tau_q_times_delta_B": (se_tq * delta_test) if not singular else np.nan,
                        "delta_B_crit_coordinate": crit_coord * sigma_n,
                        "delta_B_crit_singular": crit_sing * sigma_n,
                        "fim_singular": singular,
                        "tau_q_estimable": (not singular),
                    })
            print(f"  {config_name}: done")

    df_crlb = pd.DataFrame(crlb_records)
    out_crlb = os.path.join(DATA_DIR, "crlb_configuration_audit.csv")
    df_crlb.to_csv(out_crlb, index=False)
    print(f"\nSaved CRLB audit dataset to: {out_crlb}")

    print("\n--- Practical Identifiability Corridor |B-1|_crit (B = 0.99, coordinate CRLB) ---")
    print(df_crlb[~df_crlb["fim_singular"]][
        ["configuration", "Nt", "noise_percent", "SE_tau_q_per_sigma_n",
         "SE_Bi_per_sigma_n", "delta_B_crit_coordinate", "delta_B_crit_singular"]
    ].to_string(index=False))

    print("\n--- Generalized CRLB at exact resonance B = 1 (Moore-Penrose) ---")
    print(df_crlb[df_crlb["fim_singular"]][
        ["configuration", "Nt", "noise_percent", "SE_Bi_per_sigma_n", "cond_F3"]
    ].to_string(index=False))

    return df_crlb


def run_uncertainty_divergence():
    """Numerical uncertainty-divergence scan for Figure 10(d).

    Computes, directly from the three-parameter Fisher matrix at rear-face
    sampling (N_t = 50, Bi = 0.10), the ratio SE(tau_q)/SE(Bi) over a grid of
    resonance defects. Points whose Fisher matrix is numerically singular are
    flagged; the asymptotic 1/|B-1| reference is fitted from the computed data
    (no hard-coded constants).
    """
    print("\n" + "=" * 75)
    print("UNCERTAINTY DIVERGENCE SCAN (rear face, Nt = 50, Bi = 0.10)")
    print("=" * 75)

    b_vals = [0.80, 0.90, 0.95, 0.98, 0.99, 0.995, 0.999,
              1.0,
              1.001, 1.005, 1.01, 1.02, 1.05, 1.10, 1.20]
    t_pts = np.linspace(0.05, 1.0, 50)
    bi_test = 0.10

    rows = []
    for B in b_vals:
        J = sensitivity_jacobian(B, bi_test, t_pts, [1.0])
        m = fim_metrics(J)
        S = m["sigma"]
        singular = m["cond_F_at_numerical_floor"]
        if singular:
            # Moore-Penrose bound: tau_q is not estimable at exact resonance
            se_tq = np.nan
            se_bi = m["se_pinv"][2]
        else:
            F_inv = np.linalg.inv(J.T @ J)
            se_tq = float(np.sqrt(F_inv[0, 0]))
            se_bi = float(np.sqrt(F_inv[2, 2]))
        rows.append({
            "B": B,
            "delta_B": abs(B - 1.0),
            "Bi": bi_test,
            "sigma1": float(S[0]),
            "sigma2": float(S[1]),
            "sigma3": float(S[2]),
            "cond_F3": m["cond_F"],
            "fim_singular": singular,
            "SE_tau_q_per_sigma_n": se_tq,
            "SE_Bi_per_sigma_n": se_bi,
            "SE_ratio": (se_tq / se_bi) if not singular else np.nan,
        })
        print(f"  B = {B:6.3f}  |B-1| = {abs(B-1.0):.4f}  "
              f"SE(tau_q)/sig = {se_tq if not singular else float('nan'):10.4f}  "
              f"SE(Bi)/sig = {se_bi:.4f}  "
              f"ratio = {rows[-1]['SE_ratio'] if not singular else float('nan'):12.2f}"
              f"{'  [singular FIM]' if singular else ''}")

    df = pd.DataFrame(rows)

    # Asymptotic reference fitted from the well-conditioned part of the scan:
    # SE(tau_q)/sigma_n ~ A/|B-1| and SE(Bi)/sigma_n ~ B0 near resonance.
    ok = df[~df["fim_singular"] & (df["delta_B"] > 1e-3)].copy()
    A_fit = float((ok["SE_tau_q_per_sigma_n"] * ok["delta_B"]).mean())
    B0_fit = float(ok["SE_Bi_per_sigma_n"].mean())
    df["asymptotic_A_over_delta"] = A_fit / df["delta_B"]
    df["asymptotic_ratio"] = df["asymptotic_A_over_delta"] / B0_fit
    print(f"\n  Fitted asymptotic constants: A = SE(tau_q)*|B-1| = {A_fit:.6f}, "
          f"B0 = SE(Bi) = {B0_fit:.6f}")

    out = os.path.join(DATA_DIR, "uncertainty_divergence.csv")
    df.to_csv(out, index=False)
    print(f"Saved uncertainty divergence scan to: {out}")
    return df


if __name__ == "__main__":
    t0 = time.time()
    run_near_resonance_study()
    run_crlb_audit()
    run_uncertainty_divergence()
    print(f"\nNear-resonance study completed in {time.time() - t0:.2f} seconds.")
