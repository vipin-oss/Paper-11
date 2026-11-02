# Invariance of the Guyer–Krumhansl Fourier-Resonance Sensitivity Singularity to Boundary Heat Loss in Laser Flash Testing

**Lead Author:** Dr. Vipin Gupta
*Department of Mathematics, Gurugram University, Gurugram, Haryana 122018, India*
*Email:* `vipin.gupta@gurugramuniversity.ac.in`

---

## 1. What this package is

This is the complete, self-contained research package for **Paper 11**. It contains the
manuscript, the full numerical source code, every archived dataset that the manuscript
quotes, all publication figures, a supplementary-material document, and a machine-readable
reproducibility manifest.

The work is a **theoretical and computational study**. It is an independent mathematical and
numerical analysis of parameter identifiability in the linear one-dimensional
Guyer–Krumhansl (GK) model of non-Fourier heat conduction under laser flash boundary
conditions. **No physical experiments were performed, and no experimental data are
reported.** All results are either proved analytically or obtained by high-precision
numerical computation, and every quantitative claim in the manuscript is traceable to a
shipped dataset (see §7 below).

### Headline result

At exact Fourier resonance, $B \equiv \kappa^2/(\alpha_0 \tau_q) = 1$, the linear 1D GK slab
with linear Robin (convective/radiative) boundary heat loss is *mathematically identical*
to classical Fourier diffusion: the GK temperature field coincides with the Fourier solution
at every interior point and every time. Because the two constitutive parameters enter the
transfer function only through the combination $\kappa^2/\tau_q$ at $B=1$, the local
sensitivity vectors obey $\mathbf{J}_{\tau_q} \equiv -\alpha_0 \mathbf{J}_{\kappa^2}$, the
Fisher Information Matrix (FIM) for the estimation space $(\tau_q, \kappa^2, Bi)$ has an
exact structural null direction, and the non-Fourier pair $(\tau_q, \kappa^2)$ is **not
jointly identifiable** at $B=1$ — no matter how much, or how precisely, the transient is
measured. Boundary heat loss does **not** repair this singularity.

---

## 2. Key scientific contributions

1. **Analytical heat-loss invariance (Theorem 1).** For the linear 1D GK model in a finite
   slab with linear Robin conditions, $B = 1$ collapses both the bulk Helmholtz propagation
   operator $m(s)\to\sqrt{s}$ and the dynamic surface thermal-conductivity operator
   $\lambda_{\mathrm{eff}}(s)\to\lambda_0$ to their Fourier counterparts simultaneously.
   The GK field therefore equals the Fourier solution for all $x\in[0,1]$, $t\ge 0$ and all
   $Bi_0, Bi_L \ge 0$.

2. **Modal completeness in the Laplace domain.** In the open convergence half-plane
   $\operatorname{Re}(s)>0$ the physical IBVP has exactly two degrees of freedom
   ($C_3(s)\equiv 0$); there are no unconstrained homogeneous flux modes.

3. **Exact sensitivity collinearity and structural null direction.** The chain rule gives
   $\mathbf{J}_{\tau_q}(t)\equiv-\alpha_0\mathbf{J}_{\kappa^2}(t)$ point-by-point. In
   **dimensionless** variables $(\tau_q, \kappa^2, Bi)$ (i.e. after scaling both
   non-Fourier parameters by $\alpha_0$) the structural null vector of the 3-parameter FIM
   is
   $$\bm{v}_{\mathrm{null}} = (1,\,1,\,0)^{\top},$$
   with $\lambda_3\equiv 0$ and $\operatorname{cond}(\mathbf{F})=\infty$. In the
   **dimensional** parameter set $(\tau_q, \kappa^2, Bi)$ the same null direction is
   $\bm{v}_{\mathrm{null}} = (1,\,\alpha_0,\,0)^{\top}$; the two forms describe the same
   geometric fact and are used consistently throughout the manuscript, the derivations, the
   code, and the figure captions.

4. **Observable-subspace conditioning and generalized CRLB.** The observable 2D subspace is
   well conditioned ($\sigma_1/\sigma_2\in[2.34,\,4.90]$). Because $Bi$ is orthogonal to
   $\bm{v}_{\mathrm{null}}$, $Bi$ remains an estimable functional: the generalized
   (Moore–Penrose) CRLB gives $\operatorname{SE}(\widehat{Bi})/\sigma_n = 0.1243$ for
   rear-face monitoring with $N_t=50$ time samples at $B=1$, whereas
   $\operatorname{SE}(\widehat{\tau_q})$ and $\operatorname{SE}(\widehat{\kappa^2})$ are
   **not estimable** there.

5. **Two-sided near-resonance scaling.** The smallest singular value grows linearly away
   from resonance, $\sigma_{\min}\approx C(Bi)\,|B-1|^{1.0}$. For the 3-parameter Jacobian
   at $Bi=0.10$ the fitted one-sided prefactors are $C=1.481$ ($B<1$, exponent
   $1.00592\pm0.00224$, $R^2=0.999980$) and $C=1.426$ ($B>1$, exponent
   $0.99858\pm0.00024$, $R^2=1.000000$). The corresponding 2-parameter subsystem values
   ($C=1.512$ / $1.436$) are archived alongside them in
   `data/scaling_regression_twosided.csv`.

6. **Multi-configuration CRLB corridors.** A 54-row audit (27 near-resonance configurations
   plus 27 at $B=1$) quantifies the corridor $|B-1|_{\mathrm{crit}}$ needed for
   $\operatorname{SE}(\widehat{\tau}_q)/\tau_{q,0}\le 10\%$. For rear-face monitoring with
   $N_t=50$ the corridor is $|B-1|\ge 0.0917$ (0.1% noise), $0.459$ (0.5%) and $0.917$
   (1.0%); dual-surface monitoring tightens it to $0.0151$ / $0.0753$ / $0.1505$, and
   five-station monitoring to $0.0126$ / $0.0629$ / $0.1259$. Close to resonance the
   relative uncertainty of $\widehat{\tau}_q$ diverges as $|B-1|^{-1}$: the numerically
   evaluated ratio $\operatorname{SE}(\widehat{\tau}_q)/\operatorname{SE}(\widehat{Bi})$
   grows from $6.27$ at $|B-1|=0.2$ to $1482$ at $|B-1|=10^{-3}$.

7. **Independent numerical validation.** Every analytical statement is checked against
   independent numerical channels: a time-domain conservative finite-volume BDF PDE solver,
   Talbot-contour and de Hoog Laplace inversion at 30–80 digit working precision, the Cowan
   50-mode Fourier series, and a central finite-difference sensitivity channel.

---

## 3. Scope, assumptions and limitations

These are stated here and expanded in the manuscript's limitations discussion.

- **Model class.** Linear, constant-property, 1D GK conduction in a homogeneous slab with
  linear (Newton-type) Robin boundary conditions and instantaneous uniform pulse heating.
  Temperature-dependent properties, non-linear radiative loss ($T^4$), finite pulse
  duration, multi-layer or finite-geometry effects, and measurement-instrument transfer
  functions are **not** modelled.
- **Resonance is a knife-edge, not a neighbourhood.** The non-identifiability statement is
  exact at $B=1$. Away from resonance the pair $(\tau_q,\kappa^2)$ *is* identifiable in
  principle, but the required corridor $|B-1|_{\mathrm{crit}}$ may be far outside what a
  physical experiment can resolve; the CRLB corridors in §2(6) quantify exactly how far.
- **Noise model.** Additive, independent, identically distributed Gaussian noise with
  constant variance $\sigma_n$ on the sampled thermogram. Correlated, heteroscedastic,
  baseline-drift, or pulse-shape noise is not analysed.
- **Estimation setting.** The CRLB is a *local* bound for unbiased estimators; it does not
  address estimator bias, non-convexity of the objective, or practical convergence of
  non-linear regression. At $B=1$ the generalized (Moore–Penrose) CRLB is reported for the
  estimable functionals only.
- **Condition numbers at resonance.** The analytical condition number of $\mathbf{F}$ at
  $B=1$ is infinite. The values stored in `data/canyon_3d_landscape.csv` and plotted in
  Figure 8 are the finite values produced by the IEEE 754 double-precision evaluation
  channel (a numerical floor of order $10^{16}$–$10^{17}$); they are flagged by the
  `cond_F_at_numerical_floor` column and are **not** to be read as physical condition
  numbers. The same convention applies to the collinearity defect: the analytic channel
  gives exactly zero, while the finite-difference channel bottoms out at a round-off floor
  of order $10^{-5}$ (peak) rather than exactly zero.
- **Metrology context.** ASTM E1461-13(2022) prescribes flash-method diffusivity
  measurement geometry and data reduction; it does not prescribe an inverse-identification
  algorithm for non-Fourier constitutive parameters. This study is an independent
  theoretical analysis of such an identification and does not extend, modify, or rely on
  the standard.
- **Novelty statement.** To the best of the authors' literature search (September 2026,
  documented in `INTERNAL_AUDIT/literature/`), the combination of an exact
  heat-loss-invariant resonance singularity, its structural null direction, and the
  FIM/CRLB quantification of the resulting loss of practical identifiability has not been
  reported together for the GK flash problem.

---

## 4. Directory architecture

```
Paper11_Q1_SUBMISSION_READY/
├── README.md                           # This file
├── REPRODUCIBILITY.md                  # Environment, runtime, byte-level reproduction protocol
├── CHANGELOG_FINAL_REPAIR.md           # Every scientific value changed in the final repair pass
├── manuscript/                         # Canonical LaTeX manuscript
│   ├── main.tex                        # Master LaTeX source (single canonical version)
│   ├── main.pdf                        # Compiled manuscript (47 pages, 64 references)
│   ├── references.bib                  # DOI-verified BibTeX bibliography
│   ├── main.bbl                        # Pre-compiled bibliography (for clean-directory builds)
│   ├── elsarticle.cls                  # Elsevier article class (bundled)
│   ├── elsarticle-num.bst              # Elsevier numerical bibliography style
│   ├── compile.sh                      # One-touch pdflatex + bibtex build script
│   └── figures/                        # Figure files used by the manuscript
├── supplementary/
│   ├── supplementary_material.tex      # Supplementary tables S1–S6 (generated from data/)
│   ├── supplementary_material.pdf      # Compiled supplementary (6 pages)
│   └── elsarticle.cls / elsarticle-num.bst
├── program/                            # Complete numerical pipeline
│   ├── solver_gk_heatloss.py           # Conservative-FVM BDF time-domain PDE solver
│   ├── analytic_laplace.py             # Laplace transfer functions, Cowan series, precision studies
│   ├── independent_solver.py           # Talbot / Cowan / conservative-FVM independent benchmarks
│   ├── generate_3d_datasets.py         # 3D spatiotemporal datasets (analytic sensitivities)
│   ├── run_validation.py               # Benchmarks, asymptotic Fourier limit, mesh & precision studies
│   ├── run_heatloss_study.py           # Resonance heat-loss sweep and sensitivity profiles
│   ├── run_near_resonance.py           # Near-resonance SVD landscape, two-sided regression, CRLB audit
│   ├── run_field_verification.py       # 50-digit field invariance, interior refinement, FD step study
│   ├── plot_figures.py                 # All 10 publication figures (PDF + 300-DPI PNG)
│   ├── make_supplementary.py           # Generates supplementary_material.tex from the CSVs
│   ├── reproduce_all.py                # Master one-touch end-to-end runner (7 stages)
│   ├── requirements.txt                # Pinned Python environment
│   └── README.md                       # Program-level documentation
├── data/                               # 19 authoritative CSV datasets (every manuscript number)
├── figures/                            # Master figure suite (vector PDF + 300-DPI PNG)
├── overleaf/                           # Standalone Overleaf mirror (identical sources)
├── derivations/                        # Analytical proofs and derivations (LaTeX + PDF + Markdown)
├── run_manifest.json                   # Machine-readable execution log of the last full run
└── INTERNAL_AUDIT/                     # Internal quality-assurance records — NOT for reviewers
    ├── PACKAGE_MANIFEST.md
    ├── final_audit_report.md
    ├── RED_FLAG_FIX_REPORT.md
    ├── research_blueprint.md
    └── literature/                     # Search logs, novelty ledger, claim-to-reference matrices
```

`INTERNAL_AUDIT/` contains internal review and audit material only. It is **not** part of
the reviewer-facing manuscript, supplementary material, or code, and it is not cited by any
reviewer-facing file.

---

## 5. Publication figure suite (10 multi-panel figures)

All figures ship as vector PDF and 300-DPI PNG:

| Figure | Content |
| --- | --- |
| `fig1_heatloss_temperature_response` | Transient thermograms, heat-loss distortion, peak-amplitude ratio and half-rise-time distortion vs $Bi$ |
| `fig2_spatiotemporal_3d_temperature` | 3D spatiotemporal temperature fields and full-field GK/Fourier identity ($\sim10^{-8}$) |
| `fig3_validation_dual_benchmarks` | Talbot vs Cowan 50-mode series, Laplace discrepancy $\sim10^{-5}$, de Hoog Padé floor $10^{-56}$–$10^{-54}$, off-resonance linear recovery $C_{\mathrm{asymp}}=3.928$ |
| `fig4_field_verification_and_resolution` | Conservative-FVM second-order convergence ($p=2.00$), depth-wise thermograms, error distribution, FD step-size floor |
| `fig5_spatiotemporal_sensitivity_landscapes` | 3D sensitivity landscapes, mirror anti-symmetry, finite-difference precision-floor panel |
| `fig6_sensitivity_heatmap_and_partition` | FIM correlation heatmap, sensitivity $L_2$ norms, multi-sensor information partition, $\operatorname{tr}(\mathbf{F})$ scaling |
| `fig7_invariance_metrics_vs_biot` | Collinearity and conditioning invariance metrics over $Bi\in[0.001,0.50]$ |
| `fig8_3d_singularity_canyon_and_spectrum` | 3D condition-number canyon, iso-conditioning contours, 3-parameter SVD spectrum, observable-subspace ratio |
| `fig9_near_resonance_scaling` | Condition-number cross-sections, two-sided log–log regressions, $C(Bi)$ prefactor decline, power-law residual audit |
| `fig10_crlb_identifiability_maps` | CRLB identifiability phase diagrams, corridor vs noise and $N_t$, non-Fourier uncertainty divergence |

Sensitivities in Figures 2, 5, 6 and 8 are obtained by **analytic differentiation** of the
closed-form transfer function in 40-digit arithmetic; an independent central
finite-difference channel ($h_{\mathrm{rel}}=10^{-4}$) is used to quantify the truncation
floor rather than as the primary sensitivity channel.

---

## 6. Data archive (19 CSV datasets)

| Dataset | Content |
| --- | --- |
| `validation_laplace.csv` | Time-domain PDE vs continuous Laplace benchmark (8 Biot numbers) |
| `validation_cowan.csv` | Time-domain PDE vs exact Cowan Fourier series (5 Biot numbers) |
| `asymptotic_fourier_convergence.csv` | Off-resonance classical-limit convergence (7 decades) |
| `mesh_convergence.csv` | Conservative-FVM grid refinement ($p=2.00$) |
| `dehoog_convergence.csv` | Laplace-inversion precision table (dps 30–80) |
| `heatloss_parameter_sweep.csv` | Resonance heat-loss sweep, 3-parameter SVD (12 $Bi$) |
| `sensitivity_profiles.csv` | Time-resolved sensitivity trajectories (3 $Bi$) |
| `near_resonance_scaling.csv` | Near-resonance landscape, 2- and 3-parameter spectra (102 rows) |
| `off_resonance_sweep.csv` | Wide off-resonance conditioning sweep (51 rows) |
| `scaling_summary.csv` | Mean scaling prefactors and condition prefactors (6 $Bi$) |
| `scaling_regression_twosided.csv` | Two-sided regression audit, both singular-value definitions |
| `crlb_configuration_audit.csv` | 3-parameter CRLB audit, 54 rows (27 near-resonance + 27 at $B=1$) |
| `uncertainty_divergence.csv` | Numerical $\operatorname{SE}(\tau_q)/\operatorname{SE}(Bi)$ divergence scan (15 $B$ values) |
| `field_invariance_asymmetric.csv` | 50-digit field-level invariance residuals |
| `interior_point_mesh_refinement.csv` | Interior-point PDE residuals (20 evaluations) |
| `sensitivity_step_size_convergence.csv` | Finite-difference step truncation floor (8 steps) |
| `spatiotemporal_3d_temperature.csv` | 3D spatiotemporal temperature grid (750 points) |
| `spatiotemporal_3d_sensitivities.csv` | 3D sensitivity grid, analytic + FD channels (500 points) |
| `canyon_3d_landscape.csv` | 3D condition-number landscape (189 points) |

---

## 7. Reproduction instructions

### Figures only (~1 minute)
```bash
python3 program/plot_figures.py
```

### Supplementary tables only (seconds)
```bash
python3 program/make_supplementary.py
cd supplementary && pdflatex supplementary_material.tex && pdflatex supplementary_material.tex
```

### Full end-to-end pipeline (~11 minutes)
```bash
python3 program/reproduce_all.py
```
This re-runs all seven stages (validation → heat-loss sweeps → near-resonance landscape →
field verification → 3D datasets → figures → supplementary tables) and writes
`run_manifest.json` with per-stage runtimes, the software environment, and
artifact-completeness checks. The pipeline was verified to regenerate 18 of the 19 shipped
CSV datasets byte-for-byte; `data/mesh_convergence.csv` is scientifically deterministic
except for its wall-clock `runtime_sec` field. All 10 figures regenerate byte-for-byte apart
from their PDF creation timestamps. See `REPRODUCIBILITY.md` for the exact protocol and the verified
environment.

### Manuscript compilation
```bash
cd manuscript
bash compile.sh
```
**Required tools:** a TeX installation providing both `pdflatex` and `bibtex`, for example on
Debian/Ubuntu:

```bash
sudo apt-get install texlive-latex-base texlive-latex-recommended \
                     texlive-fonts-recommended texlive-latex-extra
```

`compile.sh` runs `pdflatex → bibtex → pdflatex → pdflatex` and leaves a clean 47-page
`main.pdf` with no undefined references and no overfull boxes. The script checks for both
`pdflatex` and `bibtex` before doing anything: if either is missing it prints the install
instructions and exits with status 1 **without producing or modifying any PDF** — it never
silently reuses a stale `main.pdf`. The pre-compiled `main.bbl` is also shipped so that a
`pdflatex`-only build succeeds if you prefer not to run `bibtex`.

### Overleaf
Create a new blank Overleaf project and upload the contents of `overleaf/`
(`main.tex`, `references.bib`, `main.bbl`, `elsarticle.cls`, `elsarticle-num.bst`,
`compile.sh`, `figures/`), then press **Recompile**. The Overleaf mirror is
byte-identical to `manuscript/` for every source file; there is exactly **one** canonical
manuscript in this package.

---

## 8. Verification and audit trail

- **Reference integrity.** All 64 bibliography entries were re-verified against the
  Crossref metadata API (DOI → title, journal, volume, pages, year) in September 2026;
  four incorrect DOIs and two inaccurate titles were corrected. No DOI, author, or
  citation in the package is fabricated.
- **Numerical reproducibility.** The complete pipeline re-runs end-to-end and regenerates
  every shipped artifact. Verified September 2026 with Python 3.13.14, NumPy 2.3.5,
  SciPy 1.17.1, mpmath 1.3.0, Matplotlib 3.10.9, pandas 2.2.3.
- **Claim traceability.** Every quantitative value quoted in the manuscript is either
  derived analytically (`derivations/`) or traceable to a shipped CSV dataset (`data/`).
  Where the manuscript states a rounded value, the full-precision value is in the
  corresponding CSV.
- **Final repair pass.** All changes made during the final submission-ready repair —
  including every scientific value that was revised (2-parameter → 3-parameter CRLB,
  asymptotic-prefactor residuals, condition-number conventions) — are recorded with their
  old and new values in `CHANGELOG_FINAL_REPAIR.md`.

---

## 9. Data availability

All datasets, code, figures, and documents in this package are provided with the submission
materials. Upon acceptance they will be archived in a public repository and a persistent
identifier will be added to the manuscript's data-availability statement. No public archive
link is claimed in the manuscript at this stage.

---

## 10. Citation

Gupta, V. *Invariance of the Guyer–Krumhansl Fourier-Resonance Sensitivity Singularity to
Boundary Heat Loss in Laser Flash Testing.* Manuscript prepared for submission to a
peer-reviewed journal in computational heat transfer / applied mathematics. Research
package: `Paper11_Q1_SUBMISSION_READY`.

---

## 11. Contact

Dr. Vipin Gupta — Department of Mathematics, Gurugram University, Gurugram, Haryana 122018,
India — `vipin.gupta@gurugramuniversity.ac.in`

## Arena round-2 run order (added 2026-11-02)

Run the data-computation scripts **before** `program/plot_figures.py`:

```bash
python3 program/compute_fig6_data.py --outdir data   # Fig. 6 norms/trace/correlation CSVs
python3 program/plot_figures.py                      # figures/fig1…fig10
```

`compute_fig6_data.py` regenerates `fig6_true_norms_trace_Nt{10,50}.csv` and
`fig6_true_correlation.csv`; `plot_figures.py` must read these CSVs, so it has to run after them.
