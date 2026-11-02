"""
compute_crlb_corridor.py — DIRECT full-FIM critical-corridor computation (Arena round 3).

WHY: Table 5 of the manuscript quotes |B-1|_crit values obtained by LOCAL LINEAR
EXTRAPOLATION from a single test point B = 0.99 (run_near_resonance.py:
crit_coord = se_tq*delta_test/(target_rel*tau_q_0)). Those are "local asymptotic
estimates", not measured thresholds. This script recomputes every Table 5
configuration by DIRECT full-FIM root-finding on a two-sided B scan.

EXACT DEFINITIONS:
  Base parameters: tau_q_0 = 0.02, kappa2 = B*tau_q_0, Bi0 = BiL = Bi = 0.10,
    tau_Delta = 0.04, raised-cosine pulse q(t) = 1 - cos(2*pi*t/tau_Delta).
  Time grid: t = linspace(0.05, 1.0, Nt)  (same as run_near_resonance.py).
  Stations / configurations (exactly the Table 5 rows):
    "Rear-face only (x=1.0)"        x = [1.0],                 Nt in {10, 50, 100}
    "Front + Rear (x=0, 1)"         x = [0.0, 1.0],             Nt = 50
    "Multi-location (5 stations)"   x = [0, .25, .5, .75, 1.0], Nt = 50
  Noise model: homoscedastic Gaussian, sigma_n in {1e-3, 5e-3, 1e-2}
    (0.1%, 0.5%, 1.0%); F = J^T J / sigma_n^2.
  Jacobian: J = [J_tau_q, J_kappa2, J_Bi], each column the ANALYTIC derivative
    (closed-form chain rule, outer mpmath dps=40) of the Laplace transfer function,
    inverted per (x, t) by de Hoog (degree=20) — same transfer function and inversion as
    compute_fig6_data.py. No finite differencing in the main channel.
  B scan: B in (0, 2), both sides of 1, deltas |B-1| in
    geomspace(1e-4, 0.1, 10) plus {0.15,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,
    0.95,0.99,0.999} on each side: B=[0.001,0.9999] and [1.0001,1.999].
    B=1 is excluded (tau_q is structurally not estimable). A finite scan
    does NOT certify absence of crossings between unsampled points.
    Scan all brackets and retain every crossing; select nearest to resonance.
    brentq xtol=2e-8, rtol=1e-12; accepted relative-error residual <2e-6.
    de Hoog degree=20 internally sets 27 working digits; outer dps=40
    is NOT 40-digit inversion accuracy. FD degree=38 internally uses 52 digits.
    Analytic chain-rule channel is checked against the original mpmath.diff.
  Criteria (solved separately, reported in separate columns):
    (param)   marginal CRLB:  SE(tau_q) = sqrt([F^-1]_11) at nonsingular B;
              root of SE(tau_q)/tau_q_0 - 0.10 = 0 via brentq, B<1 and B>1 separately.
    (singular) smallest-singular-value bound: sigma_n/sigma_min(J)/tau_q_0 - 0.10 = 0
              via brentq, both sides separately.
  If no sign change is found in the scanned range on a side, the entry is
  "not attained in scan" — NO extrapolation and no claim about all (0,2).
  Every entry is labelled method = "direct full-FIM".

INDEPENDENT CHECK: central finite-difference Jacobian (dps=50, h=1e-8) at
  B = 0.95 and B = 1.05 for two configurations (rear Nt=50; front+rear Nt=50),
  compared against the analytic channel (rel. diff of SE(tau_q)/sigma_n and
  sigma_min), and the full Jacobian against the original mp.diff channel.
  The singular criterion is worst-direction SE in the UNRESCALED parameter
  coordinates, normalized by tau_q_0; it is NOT a marginal tau_q CRLB.

OUTPUTS (relative to this script, into ../round3/data/):
  crlb_B_scan.csv            per (config, Nt, B): SE_tau_q_per_sigma_n, sigma_min, cond
  crlb_corridor_direct.csv   per (config, Nt, noise, criterion, side): |B-1|_crit, B_crit, method
  table5_comparison.csv      old Table 5 value (local asymptotic estimate) vs new direct value
  fig10a_direct.csv          corridor vs noise at Nt=50 (Fig. 10(a) replacement data)
  fd_crosscheck.csv          FD vs analytic check results
"""
import os
import argparse
from pathlib import Path
from functools import lru_cache
from concurrent.futures import ProcessPoolExecutor
import numpy as np
import pandas as pd
import mpmath as mp
from scipy.optimize import brentq

# ----------------------------------------------------------------------------
# Closed-form Laplace-domain transfer function and analytic Jacobian
# (identical to the master program/run_near_resonance.py channel).
# ----------------------------------------------------------------------------
DPS = 40
DEGREE = 20

def transfer_field(s, x, tau_q, kappa2, Bi0, BiL, tau_Delta):
    s = mp.mpc(s)
    u_sq = s * (1.0 + tau_q * s) / (1.0 + kappa2 * s)
    m = mp.sqrt(u_sq)
    mu = (1.0 + kappa2 * s) / (1.0 + tau_q * s)
    m_mu = m * mu
    Delta = ((m_mu) ** 2 + Bi0 * BiL) * mp.sinh(m) + m_mu * (Bi0 + BiL) * mp.cosh(m)
    omega = 2.0 * mp.pi / tau_Delta
    q_s = (omega ** 2 / (s * (s ** 2 + omega ** 2))) * (1.0 - mp.exp(-tau_Delta * s))
    arg = m * (1.0 - x)
    num = m_mu * mp.cosh(arg) + BiL * mp.sinh(arg)
    return (q_s / (tau_Delta * Delta)) * num

def mpdiff_jacobian(B, Bi, t_eval, x_stations, tau_q_0=0.02, tau_Delta=0.04,
                         dps=DPS, degree=DEGREE):
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

def deriv(s,x,tq,k,bi,td):
    """Closed-form chain rule; derivatives in independent (tau_q,kappa2,Bi).

    m*mu=s/m=:w, H=pulse*N/D. Differentiate H with respect to m,
    then dm/dtau_q=m*s/(2*a), dm/dkappa2=-m*s/(2*b).
    H_Bi differentiates BOTH identical Robin losses simultaneously.
    """
    a=1+tq*s; b=1+k*s
    m=mp.sqrt(s*a/b); w=s/m; z=1-x
    sh=mp.sinh(m); ch=mp.cosh(m); sz=mp.sinh(m*z); cz=mp.cosh(m*z)
    N=w*cz+bi*sz; D=(w*w+bi*bi)*sh+2*w*bi*ch
    omega=2*mp.pi/td
    pulse=omega**2/(s*(s*s+omega**2))*(-mp.expm1(-td*s))/td
    Nm=-w/m*cz+w*z*sz+bi*z*cz
    Dm=-2*w*w/m*sh+(w*w+bi*bi)*ch+2*bi*(-w/m*ch+w*sh)
    Hm=pulse*(Nm*D-N*Dm)/D**2
    Hbi=pulse*(sz*D-N*(2*bi*sh+2*w*ch))/D**2
    return Hm*m*s/(2*a),-Hm*m*s/(2*b),Hbi


def sensitivity_jacobian(B, Bi, ts, xs, degree=20):
    """Analytic chain-rule channel, sharing Laplace-node evaluations across columns."""
    mp.mp.dps=40
    rows=[]
    for x in xs:
        for t in ts:
            @lru_cache(None)
            def d(s): return deriv(s,x,.02,.02*B,Bi,.04)
            rows.append([float(mp.invertlaplace(lambda s:d(s)[i], t,method='dehoog',degree=degree)) for i in range(3)])
    return np.array(rows)

def fim_quantities(J):
    """SE(tau_q)/sigma_n = sqrt(inv(J^T J)[0,0]) and sigma_min(J)."""
    # SVD inversion of the full FIM: no truncation, no pseudoinverse.
    # This is algebraically inv(J.T @ J), avoiding squared conditioning.
    _, singular, vt = np.linalg.svd(J, full_matrices=False)
    if singular[-1] <= np.finfo(float).eps * singular[0]:
        return np.inf, float(singular[-1]), np.inf
    inv_F = (vt.T / singular**2) @ vt
    se = float(np.sqrt(inv_F[0, 0]))
    return se, float(singular[-1]), float((singular[0]/singular[-1])**2)

# ----------------------------------------------------------------------------
# FD cross-check channel (independent, dps=50)
# ----------------------------------------------------------------------------
def fd_jacobian(B, Bi, t_eval, x_stations, tau_q_0=0.02, tau_Delta=0.04, h=mp.mpf("1e-8"), dps=50):
    kappa2_0 = B * tau_q_0
    mp.mp.dps = dps
    rows = []
    for x in x_stations:
        for t in t_eval:
            def col(param, which):
                if which == "tq":
                    fp = mp.invertlaplace(lambda s: transfer_field(s, x, param + h, kappa2_0, Bi, Bi, tau_Delta), t, method="dehoog", degree=38)
                    fm = mp.invertlaplace(lambda s: transfer_field(s, x, param - h, kappa2_0, Bi, Bi, tau_Delta), t, method="dehoog", degree=38)
                elif which == "k2":
                    fp = mp.invertlaplace(lambda s: transfer_field(s, x, tau_q_0, param + h, Bi, Bi, tau_Delta), t, method="dehoog", degree=38)
                    fm = mp.invertlaplace(lambda s: transfer_field(s, x, tau_q_0, param - h, Bi, Bi, tau_Delta), t, method="dehoog", degree=38)
                else:
                    fp = mp.invertlaplace(lambda s: transfer_field(s, x, tau_q_0, kappa2_0, param + h, param + h, tau_Delta), t, method="dehoog", degree=38)
                    fm = mp.invertlaplace(lambda s: transfer_field(s, x, tau_q_0, kappa2_0, param - h, param - h, tau_Delta), t, method="dehoog", degree=38)
                return float(mp.re((fp - fm) / (2 * h)))
            rows.append([col(tau_q_0, "tq"), col(kappa2_0, "k2"), col(Bi, "bi")])
    return np.array(rows)

# ----------------------------------------------------------------------------
CONFIGS = [
    ("Rear-face only (x=1.0)", [1.0], [10, 50, 100]),
    ("Front + Rear (x=0, 1)", [0.0, 1.0], [50]),
    ("Multi-location (5 stations)", [0.0, 0.25, 0.5, 0.75, 1.0], [50]),
]
NOISE = [(0.1, 1e-3), (0.5, 5e-3), (1.0, 1e-2)]
TAU_Q0, BI, TARGET_REL = 0.02, 0.10, 0.10
DELTAS = sorted(set(np.geomspace(1e-4, .1, 10).tolist() +
                    [.15, .2, .3, .4, .5, .6, .7, .8, .9, .95, .99, .999]))
# Old Table 5 values (local asymptotic estimate from single point B=0.99), for comparison only
OLD_TABLE5 = {
    ("Rear-face only (x=1.0)", 10, 0.1): (0.246, 0.348), ("Rear-face only (x=1.0)", 10, 0.5): (1.228, 1.739),
    ("Rear-face only (x=1.0)", 10, 1.0): (2.455, 3.478),
    ("Rear-face only (x=1.0)", 50, 0.1): (0.092, 0.130), ("Rear-face only (x=1.0)", 50, 0.5): (0.459, 0.648),
    ("Rear-face only (x=1.0)", 50, 1.0): (0.917, 1.296),
    ("Rear-face only (x=1.0)", 100, 0.1): (0.065, 0.091), ("Rear-face only (x=1.0)", 100, 0.5): (0.323, 0.456),
    ("Rear-face only (x=1.0)", 100, 1.0): (0.645, 0.912),
    ("Front + Rear (x=0, 1)", 50, 0.1): (0.015, 0.021), ("Front + Rear (x=0, 1)", 50, 0.5): (0.075, 0.106),
    ("Front + Rear (x=0, 1)", 50, 1.0): (0.151, 0.212),
    ("Multi-location (5 stations)", 50, 0.1): (0.013, 0.018), ("Multi-location (5 stations)", 50, 0.5): (0.063, 0.089),
    ("Multi-location (5 stations)", 50, 1.0): (0.126, 0.178),
}

def compute_config(job):
    cfg, stations, Nt = job
    times = np.linspace(.05, 1., Nt)
    cache = {}
    def quantities(B):
        B = float(B)
        if B not in cache:
            cache[B] = fim_quantities(sensitivity_jacobian(B, BI, times, stations))
        return cache[B]
    scan_rows, roots, crossings, evaluations = [], [], [], []
    for side, sign in (("B<1", -1), ("B>1", 1)):
        # Near-to-far ordering: the nearest crossing bounds the connected
        # precision-failure interval around resonance. Also retain ALL crossings.
        Bs = [1 + sign*d for d in DELTAS]
        for B in Bs:
            se, s3, cond = quantities(B)
            for pct, noise in NOISE:
                scan_rows.append(dict(configuration=cfg, Nt=Nt, B=B, side=side,
                    noise_percent=pct, sigma_n=noise, SE_tau_q=noise*se,
                    SE_tau_q_per_sigma_n=se, relative_SE_tau_q=noise*se/TAU_Q0,
                    sigma_min=s3, singular_relative_bound=noise/(s3*TAU_Q0),
                    cond_F=cond, method="direct full-FIM"))
        print(f"[scan] {cfg} Nt={Nt} {side}: {len(Bs)} points", flush=True)
        for criterion in ("param", "singular"):
            for pct, noise in NOISE:
                def g(B):
                    se, s3, _ = quantities(B)
                    # Keep BOTH criteria as relative-error residuals, not
                    # reciprocal targets: this prevents a noise-scaling error.
                    bound = noise*se/TAU_Q0 if criterion == "param" else noise/(s3*TAU_Q0)
                    return bound - TARGET_REL
                vals = [g(b) for b in Bs]
                found = []
                for a, b, fa, fb in zip(Bs[:-1], Bs[1:], vals[:-1], vals[1:]):
                    if fa == 0:
                        root = a
                    elif fb == 0:
                        root = b
                    elif fa*fb < 0:
                        root = brentq(g, min(a,b), max(a,b), xtol=2e-8, rtol=1e-12)
                    else:
                        continue
                    if any(abs(root-r)<1e-7 for r in found):
                        continue
                    found.append(root)
                    residual = g(root)
                    assert abs(residual) < 2e-6, (cfg, Nt, root, residual)
                    crossings.append(dict(configuration=cfg, Nt=Nt, noise_percent=pct,
                        criterion=criterion, side=side, B_crit=root,
                        abs_B_minus_1_crit=abs(root-1), residual=residual,
                        bracket_B_min=min(a,b), bracket_B_max=max(a,b),
                        direction_outward="failure to pass" if fa>fb else "pass to failure",
                        method="direct full-FIM"))
                nearest = min(found, key=lambda b: abs(b-1)) if found else np.nan
                if found:
                    status = "attained"
                elif min(vals) <= 0:
                    status = "criterion satisfied at scan points; no bracketed crossing"
                else:
                    status = "not attained in scan"
                roots.append(dict(configuration=cfg, Nt=Nt, noise_percent=pct,
                    sigma_n=noise, criterion=criterion, side=side,
                    abs_B_minus_1_crit=abs(nearest-1), B_crit=nearest,
                    status=status, crossing_count=len(found),
                    residual=g(nearest) if found else np.nan,
                    scan_B_min=min(Bs), scan_B_max=max(Bs),
                    method="direct full-FIM"))
                print(f"[root] {cfg} Nt={Nt} {pct}% {criterion} {side}: "
                      f"{abs(nearest-1):.9f} ({status}; crossings={len(found)})", flush=True)
    for B, (se,s3,cond) in sorted(cache.items()):
        evaluations.append(dict(configuration=cfg,Nt=Nt,B=B,
            SE_tau_q_per_sigma_n=se,sigma_min=s3,cond_F=cond,method="direct full-FIM"))
    return scan_rows, roots, crossings, evaluations


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path,
        default=Path(__file__).resolve().parents[1]/"round3"/"data")
    parser.add_argument("--workers", type=int, default=2)
    args = parser.parse_args()
    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)
    jobs = [(cfg,stations,Nt) for cfg,stations,nts in CONFIGS for Nt in nts]
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        results = list(pool.map(compute_config, jobs))
    for i, name in enumerate(("crlb_B_scan", "crlb_corridor_direct", "all_crossings", "fim_evaluations")):
        pd.DataFrame([row for result in results for row in result[i]]).to_csv(outdir/f"{name}.csv",index=False)
    df = pd.DataFrame([row for result in results for row in result[1]])
    comp = []
    for (cfg,Nt,pct), (old_param,old_sing) in OLD_TABLE5.items():
        row = dict(configuration=cfg,Nt=Nt,noise_percent=pct,
                   old_param=old_param,old_singular=old_sing,
                   old_method="local asymptotic estimate",new_method="direct full-FIM")
        sub = df[(df.configuration==cfg)&(df.Nt==Nt)&(df.noise_percent==pct)]
        for criterion in ("param", "singular"):
            for side, suffix in (("B<1","left"),("B>1","right")):
                r=sub[(sub.criterion==criterion)&(sub.side==side)].iloc[0]
                row[f"{criterion}_{suffix}"]=r.abs_B_minus_1_crit
                row[f"{criterion}_{suffix}_status"]=r.status
        comp.append(row)
    pd.DataFrame(comp).to_csv(outdir/"table5_comparison.csv",index=False)
    df[(df.Nt==50)&(df.criterion=="param")].to_csv(outdir/"fig10a_direct.csv",index=False)

    checks=[]
    for cfg, stations in (("Rear-face only (x=1.0)",[1.]),("Front + Rear (x=0, 1)",[0.,1.])):
        times=np.linspace(.05,1.,50)
        for B in (.95,1.05):
            Ja=sensitivity_jacobian(B,BI,times,stations)
            Jf=fd_jacobian(B,BI,times,stations)
            Jm=mpdiff_jacobian(B,BI,times,stations)
            sea,s3a,_=fim_quantities(Ja)
            sef,s3f,_=fim_quantities(Jf)
            row=dict(configuration=cfg,Nt=50,B=B,FD_dps=50,FD_degree=38,FD_h=1e-8,
                SE_tau_q_per_sigma_n_analytic=sea,SE_tau_q_per_sigma_n_FD=sef,
                sigma_min_analytic=s3a,sigma_min_FD=s3f,
                rel_diff_SE=abs(sea-sef)/sea,rel_diff_sigma_min=abs(s3a-s3f)/s3a,
                rel_diff_J=np.linalg.norm(Ja-Jf)/np.linalg.norm(Ja),
                analytic_vs_mpdiff_rel_J=np.linalg.norm(Ja-Jm)/np.linalg.norm(Jm))
            assert row['rel_diff_SE']<1e-4, row
            assert row['rel_diff_sigma_min']<1e-4, row
            assert row['analytic_vs_mpdiff_rel_J']<1e-7, row
            checks.append(row)
            print('[FD/mpdiff check]',row,flush=True)
    pd.DataFrame(checks).to_csv(outdir/"fd_crosscheck.csv",index=False)
    print("ALL DONE",flush=True)

if __name__ == "__main__":
    main()
