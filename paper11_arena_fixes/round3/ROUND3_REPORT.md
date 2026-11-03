# Round 3 report — 2026-11-02

**All six requested items completed.** No pre-existing tracked file was edited;
`main.tex`, prior rounds, ZIPs, and ARENA_BRIEF are unchanged. This report and
proposed text changes do not alter the current manuscript or PDF.

| Item | Status | Finding / evidence |
|---|---|---|
| 1. Full three-parameter FIM scan | Done | Exactly 15 Table 5 configuration/noise rows; 44 B points per geometry/time grid, dense near resonance; 660 noise-expanded scan records. `data/crlb_B_scan.csv`. |
| 2. Separate branch roots | Done | 26 marginal roots and 4 marginal “not attained in scan” entries. Direct brentq evaluations, no B extrapolation. Maximum absolute relative-error residual: 8.01e-9. |
| 3. Singular-value criterion | Done | Separate columns: 24 roots and 6 “not attained in scan” entries. Criterion is sigma_n/(sigma_min(J)*tau_q0)=0.10, not marginal SE(tau_q). |
| 4. Old/new comparison and labels | Done | `data/table5_comparison.csv` has 15 rows, old coordinate/singular estimates and four new branch/criterion columns. Old = “local asymptotic estimate”; new = “direct full-FIM”. Old values assert-checked against main.tex. |
| 5. Generator, clean rerun, independent checks | Done | `../program/compute_crlb_corridor.py` is self-contained, uses relative paths, and documents all grids/parameters. Seven generator CSVs match byte-for-byte on an independent clean rerun; hashes in `rerun_verification.txt`. |
| 6. Proposed manuscript/figure-data updates | Done | `crlb_text_changes.md` contains the comparison, LaTeX-ready rows, and exact line-wise old -> new changes for the abstract, Sec. 5.6, Sec. 6, and Conclusions. `data/fig10a_direct.csv` supplies 18 replacement panel-(a) points. |

## What was wrong

The former thresholds used a single B=0.99 evaluation and assumed inverse-distance
uncertainty scaling beyond the local regime. They were not verified two-sided
CRLB boundaries. Direct thresholds are asymmetric and some are absent from the scan.
An intermediate draft of this audit also had a reciprocal-target error for the
singular criterion. It was corrected before committing; the superseded data/logs
are not included. The final validation checks the exact bound formula, its noise
scaling, and that the worst-direction bound is at least the marginal bound.

## Selected marginal-CRLB results (distances |B-1|, Nt=50)

| Configuration | Noise | Old estimate | Direct B<1 | Direct B>1 |
|---|---:|---:|---:|---:|
| Rear only | 0.1% | 0.092 | 0.0863 | 0.0991 |
| Rear only | 1.0% | 0.917 | 0.5190 | not attained in scan |
| Front + rear | 0.1% | 0.015 | 0.0150 | 0.0151 |
| Front + rear | 1.0% | 0.151 | 0.1499 | 0.1555 |

For the **singular** criterion at rear Nt=50, 1% noise, the old 1.296 is replaced
by 0.5962 on B<1 and “not attained in scan” on B>1. The old rear Nt=10, 1%-noise
singular estimate 3.478 becomes “not attained in scan” on both branches.

## Verification

- FD Jacobian: dps=50, h=1e-8, de Hoog degree=38, rear and front+rear Nt=50,
  each at B=0.95 and B=1.05. Maximum relative difference in marginal SE:
  **3.50e-11**; in sigma_min: **3.74e-11**.
- The analytic chain-rule Jacobian agrees with the retained original mp.diff
  channel to relative Jacobian norm **7.73e-13** or better on those checks.
- Twenty larger-distance roots rechecked at inversion degree 30: maximum
  relative-error residual **2.60e-9**. Explicit inversion of the full noise-scaled
  FIM agrees with the untruncated SVD inverse to relative SE **3.60e-13** or better.
- Exact-resonance check: J_tau_q + J_kappa2 = 0 in the tested channel.
- `main.tex` is byte-identical to pre-round3 commit `c4ee6cea`.

Reproduction commands (activate the environment pinned in `requirements.txt`):

```sh
python paper11_arena_fixes/program/compute_crlb_corridor.py
python paper11_arena_fixes/program/compute_crlb_corridor.py --outdir .cache/round3-rerun
diff -rq paper11_arena_fixes/round3/data .cache/round3-rerun
python paper11_arena_fixes/round3/validate_and_document.py
python paper11_arena_fixes/round3/check_numerics.py
```

## Limits / intentionally not done

Scan branches are **[0.001,0.9999] and [1.0001,1.999]**; “not attained in scan” is
not a claim over every B in (0,2), much less all physical B>=0. A finite scan cannot
exclude unresolved crossings between points. All bracketed crossings are retained;
no global monotonicity or symmetric corridor is assumed. CRLB precision is a
model-based bound, not guaranteed realized estimator accuracy.

The existing Fig. 10 image was not redrawn, nor was main.tex changed or a PDF
compiled. Requested replacement panel-(a) data and text proposals are delivered.
Panels (b,c) still use extrapolation and need regeneration or explicit labelling
as local asymptotic estimates; the proposed caption uses the labelling option.
