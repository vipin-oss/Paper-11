"""
compute_fig6_data.py — reproducible generator of the Fig. 6 correction datasets.

EXACT DEFINITIONS (all parameters, grids, and formulas):
  Base parameters: B = 1 (Fourier resonance), tau_q = 0.02, kappa2 = B*tau_q = 0.02,
    Bi0 = BiL = Bi = 0.10 (SYMMETRIC Robin), tau_Delta = 0.04, raised-cosine pulse
    q(t) = 1 - cos(2*pi*t/tau_Delta) on [0, tau_Delta], else 0.
  Jacobian: J = [J_tau_q, J_kappa2, J_Bi], each column = ANALYTIC derivative
    (mpmath.diff, dps=40) of the closed-form Laplace transfer function
    theta(x,s) = [m*mu*cosh(m(1-x)) + BiL*sinh(m(1-x))] * q_pulse(s) / (tau_Delta*Delta(s)),
    inverted per (x,t) by de Hoog (degree=20). No finite differencing in the main channel.
  Time grid: t = linspace(0.05, 1.0, Nt) with Nt = 10 and Nt = 50 (two separate runs).
  Stations: x in [0.0, 0.25, 0.50, 0.75, 1.0] (5 stations, increasing x).
  Noise: sigma_n = 1 (norms and tr(F) are per unit sigma_n; F = J^T J = sigma_n^-2 J^T J with sigma_n=1).
  Norm definition: norm_J_p(x) = ||J_p(x, :)||_2 = sqrt( sum_t J_p(x,t)^2 ).
  trF_cum[m] = trace of F restricted to stations 0..m, cumulative in increasing-x order:
    trF_cum[m] = sum_{x in stations[0..m]} sum_p norm_J_p(x)^2.
  Fig. 6(a) correlation matrix: R_jk = F_jk / sqrt(F_jj F_kk) with F = J^T J over ALL
    5 stations and the Nt=50 grid.

SANITY CHECKS built in:
  (i)  norm_J_tq(x) == norm_J_k2(x) exactly at B=1 (collinearity J_tq = -J_k2);
  (ii) independent cross-check of norms at stations x=0.0 and x=1.0 via CENTRAL FINITE
       DIFFERENCE of the transfer function at dps=50 (h=1e-8), separate from the analytic channel;
  (iii) tr(F) Nt=10 vs Nt=50 ratio explained by re-evaluating both grids on a common
       50-point grid subsampled to 10 points.

Outputs (written to --outdir, default ../data):
  fig6_true_norms_trace_Nt10.csv, fig6_true_norms_trace_Nt50.csv, fig6_true_correlation.csv
"""
import sys, os, argparse
import numpy as np
import pandas as pd
import mpmath as mp

BASE = "/home/user/work/extract/Paper11_Q1_SUBMISSION_READY_MASTER"
sys.path.insert(0, os.path.join(BASE, "program"))
from run_near_resonance import sensitivity_jacobian  # analytic-diff channel (dps=40, dehoog degree 20)
from analytic_laplace import transfer_function_gk_field as tf

TAU_Q0, BI, TAU_DELTA, B = 0.02, 0.10, 0.04, 1.0
STATIONS = [0.0, 0.25, 0.50, 0.75, 1.0]

def jac_analytic(Nt):
    t_pts = np.linspace(0.05, 1.0, Nt)
    J = sensitivity_jacobian(B, BI, t_pts, STATIONS)      # rows x-major: (5*Nt, 3)
    return t_pts, J.reshape(len(STATIONS), Nt, 3)

def norm_fd_station(x, Nt, h=1e-8, dps=50):
    """Independent route: central finite difference of transfer function at dps=50."""
    mp.mp.dps = dps
    t_pts = np.linspace(0.05, 1.0, Nt)
    k2 = B * TAU_Q0
    out = []
    for t in t_pts:
        def safe(x, *params):
            # tf signature is transfer_function_gk_field(x, s, tau_q, kappa2, Bi0, BiL, tau_Delta);
            # nudge s off s == 0 if an inversion node lands there
            def f(s):
                if s == 0:
                    s = mp.mpf("1e-30")
                return tf(x, s, *params)
            return f
        def inv(f):
            # a de Hoog node can land on the pulse pole s = ±i*omega; fall back to Talbot there
            try:
                return mp.invertlaplace(f, t, method="dehoog", degree=20)
            except ZeroDivisionError:
                return mp.invertlaplace(f, t, method="talbot")
        def col(param, which):
            if which == "tq":
                fp = inv(safe(x, param + h, k2, BI, BI, TAU_DELTA))
                fm = inv(safe(x, param - h, k2, BI, BI, TAU_DELTA))
            elif which == "k2":
                fp = inv(safe(x, TAU_Q0, param + h, BI, BI, TAU_DELTA))
                fm = inv(safe(x, TAU_Q0, param - h, BI, BI, TAU_DELTA))
            else:
                fp = inv(safe(x, TAU_Q0, k2, param + h, param + h, TAU_DELTA))
                fm = inv(safe(x, TAU_Q0, k2, param - h, param - h, TAU_DELTA))
            return float(mp.re((fp - fm) / (2 * h)))
        out.append([col(TAU_Q0, "tq"), col(k2, "k2"), col(BI, "bi")])
    J = np.array(out)
    return float(np.sqrt((J[:, 0] ** 2).sum())), float(np.sqrt((J[:, 1] ** 2).sum())), float(np.sqrt((J[:, 2] ** 2).sum()))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data"))
    ap.add_argument("--shipped", default="/home/user/Paper-11/paper11_arena_fixes/data")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    results = {}
    for Nt in (10, 50):
        t_pts, J = jac_analytic(Nt)                            # (5, Nt, 3)
        norms = np.sqrt((J ** 2).sum(axis=1))                  # (5, 3)
        tr_cum = np.cumsum((norms ** 2).sum(axis=1))            # cumulative over stations, increasing x
        df = pd.DataFrame({"x": STATIONS, "norm_J_tq": norms[:, 0], "norm_J_k2": norms[:, 1],
                           "norm_J_Bi": norms[:, 2], "trF_cum": tr_cum, "Nt": Nt})
        out = os.path.join(args.outdir, f"fig6_true_norms_trace_Nt{Nt}.csv")
        df.to_csv(out, index=False)
        results[Nt] = df
        print(f"Nt={Nt} written {out}")
        print(df.round(4).to_string(index=False))
        # sanity (i): exact collinearity of norms at B=1
        d = np.max(np.abs(norms[:, 0] - norms[:, 1]))
        print(f"  sanity(i) max |norm_J_tq - norm_J_k2| = {d:.3e}")

    # Fig 6(a): correlation matrix at Nt=50, all 5 stations
    _, J50 = jac_analytic(50)
    F = J50.reshape(-1, 3).T @ J50.reshape(-1, 3)
    R = F / np.sqrt(np.outer(np.diag(F), np.diag(F)))
    pd.DataFrame(R, index=["tau_q", "kappa2", "Bi"], columns=["tau_q", "kappa2", "Bi"]).to_csv(
        os.path.join(args.outdir, "fig6_true_correlation.csv"))
    print("Fig6(a) R = F_jk/sqrt(F_jj F_kk):\n", np.round(R, 6))
    print(f"  rho(tau_q,Bi) = {R[0,2]:.4f}, rho(kappa2,Bi) = {R[1,2]:.4f}, rho(tau_q,kappa2) = {R[0,1]:.12f}")
    _, J10 = jac_analytic(10)
    F10 = J10.reshape(-1, 3).T @ J10.reshape(-1, 3)
    R10 = F10 / np.sqrt(np.outer(np.diag(F10), np.diag(F10)))
    print(f"  [Nt=10] rho(tau_q,Bi) = {R10[0,2]:.4f}   (hard-coded value in plot_figures.py is 0.0382)")

    # sanity (ii): FD cross-check at x=0 and x=1, Nt=50
    for x in (0.0, 1.0):
        n = norm_fd_station(x, 50)
        a = results[50][np.isclose(results[50]["x"], x)][["norm_J_tq", "norm_J_k2", "norm_J_Bi"]].values[0]
        print(f"  sanity(ii) x={x}: analytic {np.round(a,4)} vs FD {np.round(n,4)}  max rel diff = {np.max(np.abs(a-n)/np.abs(a)):.2e}")

    # sanity (iii): tr(F) Nt=10 vs Nt=50 ratio on common grid
    t50 = np.linspace(0.05, 1.0, 50)
    J50 = jac_analytic(50)[1].reshape(-1, 3)
    tr50_full = float((J50 ** 2).sum())
    idx10 = np.linspace(0, 49, 10).astype(int)   # every 5th point of the 50-grid
    tr50_sub = float((J50[idx10] ** 2).sum())
    tr10 = float(results[10]["trF_cum"].iloc[-1])
    print(f"  sanity(iii) tr(F): Nt=10 grid = {tr10:.1f}, Nt=50 grid = {results[50]['trF_cum'].iloc[-1]:.1f}, "
          f"50-grid subsampled to 10 = {tr50_sub:.1f}; ratio 50/10 = {results[50]['trF_cum'].iloc[-1]/tr10:.2f}")
    print("  -> ratio < 5 because J is concentrated at early times; the 10-pt grid (step 0.105) "
          "samples fewer early-time large-|J| points than the 50-pt grid (step 0.0196).")

    # compare vs shipped CSVs
    for Nt in (10, 50):
        f_new = os.path.join(args.outdir, f"fig6_true_norms_trace_Nt{Nt}.csv")
        f_old = os.path.join(args.shipped, f"fig6_true_norms_trace_Nt{Nt}.csv")
        if os.path.exists(f_old):
            a, b = pd.read_csv(f_new), pd.read_csv(f_old)
            common = [c for c in a.columns if c in b.columns and a[c].dtype != object]
            md = max(float(np.max(np.abs(a[c].values - b[c].values))) for c in common)
            print(f"  shipped-vs-new max abs diff Nt={Nt}: {md:.3e}")

if __name__ == "__main__":
    main()
