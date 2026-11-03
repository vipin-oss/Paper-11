# Round 3 — direct full-FIM corridor audit

No pre-existing manuscript, ZIP, plotting script, or ARENA_BRIEF is changed.
The new generator is `../program/compute_crlb_corridor.py`; all results and
proposed manuscript changes are in this directory.

## Reproduce (from the repository root)

```sh
python -m venv .venv
.venv/bin/pip install -r paper11_arena_fixes/round3/requirements.txt
.venv/bin/python paper11_arena_fixes/program/compute_crlb_corridor.py
# Independent clean run: no cached Jacobians or prior CSVs are read.
.venv/bin/python paper11_arena_fixes/program/compute_crlb_corridor.py --outdir .cache/round3-rerun
diff -rq paper11_arena_fixes/round3/data .cache/round3-rerun
.venv/bin/python paper11_arena_fixes/round3/validate_and_document.py
.venv/bin/python paper11_arena_fixes/round3/check_numerics.py
```

The generator defaults to two worker processes (`--workers 1` is supported).
It is CPU-intensive. All paths are script-relative or explicitly supplied by
`--outdir`; the original ZIP need not be extracted. No datasets are downloaded.
Use numerical tolerances rather than expecting byte equality across different
platforms/library builds; the recorded same-environment rerun uses byte equality.

## Definitions

- Base parameters: tau_q0=0.02, kappa2=B*tau_q0, Bi0=BiL=0.10, tau_Delta=0.04.
- Raised-cosine input `1-cos(2*pi*t/tau_Delta)` during the pulse, zero afterward.
- `t=linspace(0.05,1,Nt)`. Rear-only Nt=10,50,100; front+rear Nt=50;
  five stations [0,0.25,0.5,0.75,1] Nt=50. These are exactly the 15 Table 5 rows
  after including sigma_n=0.001,0.005,0.01 (absolute dimensionless noise SD).
- Full independent-parameter Jacobian: derivatives w.r.t. tau_q and kappa2
  hold the other parameter fixed. B labels the evaluation point, not a constraint
  imposed during differentiation. The Bi derivative changes both boundaries together.
- `F=J.T@J/sigma_n**2`. The full inverse is evaluated via untruncated SVD:
  `F^-1=sigma_n**2 * V diag(s^-2) V.T`. No rank truncation or pseudoinverse is used.
- Marginal relative SE: `sigma_n*sqrt(inv(J.T@J)[0,0])/tau_q0`.
- Distinct worst-direction relative bound: `sigma_n/(sigma_min(J)*tau_q0)`.
  This singular criterion is coordinate-dependent, in unrescaled parameter coordinates.
- Each criterion is solved at value 0.10 using brentq on each sign-changing bracket.
  All crossings are retained; the nearest to resonance is tabulated.
- B scan: `1 +/- d`, where d is the union of `geomspace(1e-4,0.1,10)` and
  [0.15,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,0.99,0.999].
  Hence [0.001,0.9999] and [1.0001,1.999], not the entire open interval (0,2).
  Exact B=1 is excluded because tau_q is structurally not estimable.
- An absent crossing is **not attained in scan**, never an extrapolation.
  A finite scan is not proof of no crossing between unresolved sample points.
  CRLB-based precision is a model-based local bound, not a guarantee of realized
  estimator error in an experiment.

## Analytic channel and checks

The default Jacobian uses closed-form chain-rule derivatives of the same Laplace
transfer function as the master `run_near_resonance.py`. This avoids repeated
mpmath.diff calls but makes no approximation in the derivative. With
`m=sqrt(s*(1+tau_q*s)/(1+kappa2*s))`, the surface factor `m*mu=s/m`;
differentiate the numerator and denominator w.r.t. m, then use
`dm/dtau_q=m*s/(2*(1+tau_q*s))` and
`dm/dkappa2=-m*s/(2*(1+kappa2*s))`.

The default channel uses outer mpmath dps=40 and de Hoog degree=20, matching the
shipped inversion degree. **Degree=20 internally sets 27 working digits**;
outer dps does not imply 40-digit inverse-Laplace accuracy.

The original mpmath.diff implementation is retained as `mpdiff_jacobian` and
cross-checked on rear-only and front+rear Nt=50, at B=0.95 and B=1.05.
An independent central FD Jacobian differentiates time-domain temperatures,
with dps=50, h=1e-8, and degree=38 (52 internal working digits). Checks compare
full Jacobians, marginal SE, and smallest singular values. Root residuals and
noise scaling are separately validated; matching CSVs alone do not establish
formula correctness.

`check_numerics.py` additionally recomputes the largest attained distance for
each configuration/Nt/branch/criterion at de Hoog degree 30 (20 checks).
It also compares the SVD inverse to an explicit inverse of the full noise-scaled
3x3 FIM; results are in `root_numerical_checks.csv`.

## Files

- `data/crlb_B_scan.csv`: 660 noise-expanded scan records (full-FIM marginal and singular bounds).
- `data/crlb_corridor_direct.csv`: 60 branch/criterion entries, explicit statuses and residuals.
- `data/all_crossings.csv`: all bracketed roots, bracketing intervals, and crossing directions.
- `data/fim_evaluations.csv`: all direct evaluations, including brentq refinement points.
- `data/table5_comparison.csv`: 15 rows, both old estimates and four separate new columns.
- `data/fig10a_direct.csv`: 18 replacement panel-(a) points, marginal criterion, Nt=50.
- `data/fd_crosscheck.csv`: two independent configuration checks, each on both sides.
- `crlb_text_changes.md`: comparison table, LaTeX rows, exact line-wise old -> new proposals.
- `table5_new_rows.tex`: directly usable table rows (NA means not attained in scan).
- `validate_and_document.py`: assertions and deterministic document generation; never edits main.tex.
- `root_numerical_checks.csv`: 20 degree-30 root/conventional full-inverse checks.
- `validation_summary.json`: machine-readable check results and manuscript source hash.

The new CSV does not redraw the existing Fig. 10. Panels (b,c) still depend on
local asymptotic extrapolation; regenerate them separately or label them as such.
No claim is made that the current manuscript/PDF is corrected until the proposals
are applied and the figures regenerated.
