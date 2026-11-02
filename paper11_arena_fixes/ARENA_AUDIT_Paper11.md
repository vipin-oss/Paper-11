# ARENA AUDIT — Paper 11
**Paper:** "Invariance of the Guyer–Krumhansl Fourier-Resonance Sensitivity Singularity to Boundary Heat Loss in Laser Flash Testing"
**Branch:** `arena/f3708631-paper-11` · **Date:** 2026-11-02 · **Working dir:** `paper11_arena_fixes/`

## 0. Inventory / version resolution
Repo root contains 3 zips only. Diff result: `Paper11_Q1_SUBMISSION_READY_MASTER.zip` is the master
(latest CHANGELOG 2026-09-28 08:48; its `manuscript/` and `program/` are byte-identical to
`Paper11_Journal_Submission.zip` and `Paper11_Reproducibility_Package.zip`). No conflicting copies;
nothing was overwritten. All fixes live in `paper11_arena_fixes/`.

## 1. Audit table

| # | Finding | Class | Evidence (command + output) | Proposed fix | Applied |
|---|---------|-------|------------------------------|--------------|---------|
| 1 | Core claim TRUE: at B=1 GK field ≡ Fourier field for all x,t,Bi (Theorem 1) | — (verified) | `verify_paper11.py`: 40-digit Laplace field collapse, 8 random (s,x,asym Bi): max rel diff = 1.66e-41; time-domain GK==Fourier==Cowan max diff = 1.20e-14 | none needed | n/a |
| 2 | Collinearity J_τq = −J_κ² and FIM null vector (1,1,0) exact at B=1 | — (verified) | mpmath analytic diff, 27 (x,t,Bi) points: residual 0.00e+00; σ3/σ1 = 1.03e-16; R_J = 0.00e+00 | none needed | n/a |
| 3 | Cross-solver agreement PDE vs de Hoog vs Talbot vs FVM | — (verified) | max \|PDE−dehoog\| = 6.85e-06, \|Talbot−dehoog\| = 0, \|FVM−dehoog\| = 6.85e-06 (B=1, Bi=0.1); same at Bi=0 and B=0.5 | none needed | n/a |
| 4 | PDE residual / Robin BCs of the staggered solver | — (verified) | dense-output derivative residual: energy 1.07e-07, flux 2.00e-07; Robin BC algebraic residual 0.00e+00 | none needed | n/a |
| 5 | Limits: τq,κ²→0 = Fourier (exact 0); κ²→0 = MCV (2e-14); B=1 collapse under asymmetric Robin (exact 0); large Bi decays | — (verified) | verify3.py limit suite | none needed | n/a |
| 6 | Near-resonance two-sided linear law σmin ≈ C\|B−1\|: exponents 1.0059/0.9986, C = 1.481/1.426 at Bi=0.10; C(Bi) ∈ [1.3478, 1.4637] | — (verified) | independent recomputation + `scaling_regression_twosided.csv`, `scaling_summary.csv` match manuscript | none needed | n/a |
| 7 | CRLB Table 8 values (27 configs at B=0.99 + 27 at B=1); SE(Bi)/σn = 0.1243 | — (verified) | recomputed SE(Bi)/σn = 0.124343; CSV rear Nt=50 0.1%: crit 0.0917/0.1296 (table 0.092/0.130) | none needed | n/a |
| 8 | Validation/mesh numbers (L∞ 1.02e-5; mesh 1.425e-4→2.225e-6, p=2.00; C_asymp→3.928) | — (verified) | `validation_laplace.csv`, `mesh_convergence.csv`, `asymptotic_fourier_convergence.csv` | none needed | n/a |
| 9 | **Fig 1(b,c,d) hard-coded / fabricated**: residual curves were constants (1.2e-14…), peak ratio from formula 1/(1+1.24Bi+0.45Bi²) giving 0.578 at Bi=0.5, halftime shift from formula −(0.12Bi)/(1+0.8Bi) giving −4.3% | **A** | `plot_figures.py` lines 84-110; true PDE values: peak ratio 0.5908, halftime shift −17.39% at Bi=0.5 | regenerate from computed data | **yes** (fig1 rebuilt; text 0.578→0.591, −4.3%→−17.4%) |
| 10 | **Fig 4(b) interior thermograms fabricated** as scaled copies of rear-face curve | **A** | `plot_figures.py` lines 274-276 | plot true interior fields from `spatiotemporal_3d_temperature.csv` | **yes** |
| 11 | **Fig 6(b,c,d) hard-coded arrays** not matching computed data (norms 42.15/…/Bi 0.0 at front; partition 42.4/19.7/12.5-12.8%; tr(F) 28.2→105.8, 82.5→332.5) | **A** | `plot_figures.py` lines 421-456; computed (Nt=50): norms 31.9/18.6/15.4, Bi 8.95→8.15; partition 40.6/26.8/14.2/8.0/10.4%; tr(F) 656→1675 (Nt=10), 2112→5209 (Nt=50) | recompute from 3-parameter Jacobian | **yes** (fig6 rebuilt; text updated) |
| 12 | Fig 3(d) asymptote hard-coded 3.928 (value itself confirmed as asymptotic ratio 3.9284) | C | `plot_figures.py` line 235; per-ε ratios → 3.928 | fit from data in plot script | **yes** |
| 13 | Data Availability paragraph referenced the reproducibility package and `program/reproduce_all.py` (against submission rules) | **B** | `main.tex` ~line 741 | standard data-availability sentence | **yes** |
| 14 | Manuscript referenced "Supplementary Table S1/S2" and a CSV source in Table 8 caption; "archived dataset" phrasing | **B** | `main.tex` lines 482, 504, 583, 596, 632, 636, 670 | remove all supplementary/code-archive references | **yes** |
| 15 | Abstract 274 words (>250) and did not start with "This paper studies" | C | word count script: 274 → rewritten 188 words | rewrite abstract | **yes** |
| 16 | Conclusions as 6 prose paragraphs; Limitations as 7 bullets | C | style guide | summary + 6 numeric bullets + application; 3 limitation bullets | **yes** |
| 17 | Mixed "non-local"/"nonlocal" | C | grep count | standardize American "nonlocal" | **yes** |
| 18 | Missing Ethics and AI-use declarations | B | declarations audit | add default statements, flag author-to-confirm | **yes** (flagged) |
| 19 | Cowan series drops the zero eigenmode at Bi=0 (validation correctly uses Bi≥0.01) | D (observation) | verify3: Cowan(Bi=0) misses steady mode | none (not claimed for Bi=0) | no |
| 20 | TeX toolchain unavailable in sandbox → patched PDF not recompiled here | D | `which pdflatex` empty; tectonic binary download blocked | author/editor compiles; patched tex is text-only diff | no (open item) |
| 21 | References: 64/64 bib entries cited, none missing; ~10 key entries verified online (Kovács 2018 IJHMT 127:631-636; Cowan 1963 JAP 34:926; Cape & Lehman 1963 JAP 34:1909; Ding 2022 Nat Commun 13:285; Huberman 2019 Science 364:375; Blumm 2025 IJT 46:39; Sundar & Sudha 2025 IJT 46:13; Sellitto 2025 ZAMP 76:105; Ván 2017 EPL 118:50005; Tröger & Hartmann 2022 GAMM 45:e202200013) | — | web search + DOI checks | remaining 54 entries marked "unverified" | partial |

## 2. Claims verified correct (independent recomputation)
- Theorem 1 B=1 full-field collapse (Laplace, 40-digit; time-domain; asymmetric Robin) — see #1.
- dθ/dτq = −dθ/dκ², FIM null direction (1,1,0), σ3 ≡ 0 (floor 1e-16), R_J = 0 — see #2.
- Three/four-solver cross-checks, PDE residual, Robin BCs, limits — see #3-#5.
- σmin two-sided linear |B−1| (exponents 1.0059 ± …, 0.9986; C = 1.481/1.426; C(Bi) ∈ [1.348, 1.464]) — see #6.
- Table 8 CRLB values, 27+27 configuration counts, SE(Bi)/σn = 0.1243 — see #7.
- Validation L∞ = 1.02e-5, mesh convergence p = 2.00 with exact error ladder, C_asymp = 3.928 — see #8.
- Dataset counts: 19 CSVs, 11 Python modules — match Data Availability counts (now removed from manuscript per rule).

## 3. Claims found wrong (old → new → reason → where updated)
| Old | New | Reason | Updated in |
|-----|-----|--------|-----------|
| peak ratio 0.578 at Bi=0.5 | 0.591 | was value of a hard-coded curve, not computed; PDE/Laplace give 0.5908 | main.tex §Fig.1 text; fig1(c) |
| halftime shift −4.3% | −17.4% at Bi=0.5 | hard-coded formula; own-peak Parker definition gives −17.39% | main.tex §Fig.1 text; fig1(d) |
| Fig6(b) norms 42.1/22.4/21.5, Bi 0.0→11.3 | 31.9/18.6/15.4, Bi 9.0→8.2 | hard-coded; recomputed 3-param Jacobian at Nt=50 | main.tex §Fig.6 text; fig6(b) |
| Fig6(c) partition 42.4/19.7/12.5-12.8% | 40.6/26.8/8.0-14.2% | recomputed | main.tex text + caption; fig6(c) |
| Fig6(d) tr(F) 28.2→105.8 / 82.5→332.5 | 656→1675 / 2112→5209 | recomputed tr(JᵀJ) | main.tex text; fig6(d) |

## 4. Deliverables in `paper11_arena_fixes/`
- `ARENA_AUDIT_Paper11.md` (this file)
- `manuscript/` — patched `main.tex`, `references.bib`, elsarticle class/bst, figures
- `program/` — full program suite with patched `plot_figures.py` (fig1/3/4/6 data-driven)
- `figures/` — regenerated fig1, fig3, fig4, fig6 (PDF+PNG)
- `data/` — correction datasets (fig1_true_*, fig6_true_*)
- Compiled PDF: **not rebuilt in sandbox (no TeX)**; shipped `main.pdf` remains valid for unpatched text; patched tex is text-only.

## 5. Likely reviewer questions (1-line answers)
1. Why does heat loss not break the degeneracy? Because at B=1 both m(s)→√s and μ(s)→1, so the whole BVP becomes Fourier — no boundary data can distinguish them.
2. Is σ3 exactly zero or just small? Exactly zero analytically; computed values (1e-8 FD, 1e-16 analytic) are channel floors, flagged as such.
3. Why is Bi estimable but τq not? The null vector (1,1,0) has zero Bi component, so the Bi functional is orthogonal to the null space (Moore–Penrose bound).
4. Scope vs ASTM E1461? The standard prescribes Fourier data reduction; this paper analyzes identifiability within the GK model class, it does not modify the standard.
5. Nonlinear radiation? Out of scope — stated explicitly in Limitations; linear Robin is the model class of Theorem 1.

## 6. Open items for the author
- Journal choice (IJHMT vs Appl. Math. Modelling vs Int. J. Thermal Sciences).
- AI-use disclosure wording (default statement added, marked to confirm).
- Declarations confirmation (funding/competing interest/CRediT defaults present).
- 54 of 64 references not individually re-verified online in this session ("unverified").
- Whether to include the supplementary PDF in the journal zip.
- Recompile PDF in a TeX environment (pdflatex→bibtex→pdflatex→pdflatex) before submission.

## 7. Journal best-fit suggestions (quartiles: verify on the journal's current Scimago page before claiming)
- **International Journal of Heat and Mass Transfer (Elsevier, Q1)** — best fit: heat-transfer identifiability, laser flash metrology; already the \journal line.
- **International Journal of Thermophysics (Springer, Q1/Q2 — verify)** — flash-method community.
- **Applied Mathematical Modelling (Elsevier, Q1)** — analytical Laplace-domain theorem + inverse-problem analysis.
- **International Journal of Thermal Sciences (Elsevier, Q1)** — alternative from the brief.
- **Mathematics and Computers in Simulation (Elsevier, Q2 — verify)** — backup.

## 8a. Exact definitions behind the new Fig 1(c,d) numbers
- **Peak ratio (Fig 1c):** θ_max(Bi)/θ_max(0) from rear-face thermograms θ(1,t) of the time-domain PDE solver (`solver_gk_heatloss.solve_gk_heatloss`, Nx=200, BDF, rtol=1e-10), B=1 (τq=κ²=0.02), τΔ=0.04, raised-cosine pulse, **symmetric Robin Bi0=BiL=Bi**, t∈[0.005,2.0], 500 samples. At Bi=0.5: peak = 0.5908 → **0.591** (old hard-coded formula value 0.578).
- **Half-rise shift (Fig 1d):** (t_{1/2} − t_{1/2,0})/t_{1/2,0} where t_{1/2} is the time at which each curve reaches **half of its own peak** (standard Parker definition, per-curve normalization), same solver/settings as above. At Bi=0.5: t_{1/2}=0.1329 vs 0.1609 → **−17.4%** (old hard-coded formula value −4.3%). Note: with the alternative absolute-level definition (half of the adiabatic peak) the sign flips to +24.8%; the manuscript figure/text now use the own-peak definition explicitly.
- **Source command:** `/home/user/work/venv/bin/python` script computing the sweep (logged in `true_fig1.log` / `fig_corrections.log`); results archived in `paper11_arena_fixes/data/fig1_true_peak_halftime.csv` and `fig1_true_identity_residual.csv`.

## 8b. References NOT individually verified online this session (54 of 64; the other 10 verified — see #21)
alifanov1994inverse, alvarez2007memory, aster2019parameter, astm2022standard, baba2001improvement, beck1977parameter, beck1985inverse, blumm2002improvement, both2016deviation, cahill2003nanoscale, cahill2014nanoscale, carslaw1959conduction, cattaneo1958forme, chen2001ballistic, chester1963second, cimmelli2009different, clark1975radiation, degiovanni1986nouvelle, dehoog1982improved, dowding1995estimation, emery2000uncertainties, fadale1995two, guo2015phonon, guyer1966solution, guyer1966thermal, heckman1973finite, jou2010extended, kovacs2015generalized, kovacs2020numerical, lapka2020development, li2020appropriate, maillet2000thermal, maillet2019review, minnich2012determining, nishi2020effect, ordonez2010determination, ozisik1993heat, ozisik2000inverse, parker1961flash, ramos2023mathematical, ramos2023nonnegativity, raue2009structural, salazar2014extending, sellitto2016mesoscopic, talbot1979accurate, tzou1995unified, tzou1997macro, vernotte1958paradoxes, vozar2004flash, walter1996identifiability, walter1997identification, wang2001wellposedness, wang2002thermal, zhu2017nonlocal
(Verified online: kovacs2018analytic, cowan1963pulse, cape1963temperature, ding2022observation, huberman2019observation, blumm2025laser, sundar2025thermal, sellitto2025nonlinear, van2017guyer, troger2022identification.)

## 8c. Old → new values for Fig 3(d), 4(b), 6(b,c,d) and manuscript update locations
| Figure panel | Old (hard-coded in plot_figures.py) | New (computed) | Manuscript text updated |
|---|---|---|---|
| Fig 3(d) asymptote | constant 3.928 drawn | fitted from data: asymptotic ratio 3.9284 (last-4-point fit) — value unchanged, now data-driven | none needed (§Asymptotic Fourier Limit, C_asymp = 3.928 still correct) |
| Fig 4(b) interior curves | fabricated scalings 0.5(1+e^{−3t})T, 0.25(1+3e^{−5t})T of rear curve | true interior fields at x=0.25/0.5/1.0 from `spatiotemporal_3d_temperature.csv` (T_loss_Bi0p2) | none (no numeric claim tied to panel) |
| Fig 6(b) J norms, Nt=50 | tq/k2: 42.15, 28.60, 22.45, 21.80, 21.50; Bi: 0.00, 4.12, 7.85, 9.95, 11.29 | tq/k2: 31.88, 25.86, 18.56, 13.41, 15.44; Bi: 8.95, 7.75, 7.24, 7.39, 8.15 | §Parameter Correlation text: 42.1/22.4/21.5 → 31.9/18.6/15.4; Bi 0.0→11.3 → 9.0→8.2 |
| Fig 6(c) energy partition | 42.4%, 19.7%, 12.5–12.8% | 40.6%, 26.8%, 14.2%, 8.0%, 10.4% | text + caption: 42.4%→40.6%, 19.7%→26.8%, 12.5–12.8%→8.0–14.2%, rear 12.6%→10.4% |
| Fig 6(d) tr(F) | Nt=10: 28.2→105.8; Nt=50: 82.5→332.5 | Nt=10: 656→1675; Nt=50: 2112→5209 | §Parameter Correlation text updated to 656→1675 and 2112→5209 |
Correction datasets: `paper11_arena_fixes/data/fig6_true_norms_trace_Nt10.csv`, `fig6_true_norms_trace_Nt50.csv`, `fig1_true_peak_halftime.csv`, `fig1_true_identity_residual.csv`.

## 8d. Build log summary
**NOT BUILT in this sandbox.** No TeX toolchain available: `pdflatex`/`bibtex` absent, apt requires root (denied), and the Tectonic binary download from GitHub release assets is blocked by the sandbox egress proxy (only github.com/codeload/api endpoints reachable; release-asset redirects fail). The patched `main.tex` is a text-only diff (75 changed lines, braces balanced, no new \cite/\ref/\includegraphics), so build risk is minimal, but the author must run `pdflatex → bibtex → pdflatex → pdflatex` locally and check 0 errors / 0 undefined / 0 overfull.

## 9. Final summary
**Verdict: READY AFTER MINOR CHANGES.**
- **Verified by independent code:** the entire core claim chain (Theorem 1 collapse, collinearity, null vector, σ3≡0, two-sided |B−1| scaling with the exact quoted exponents/slopes, Table 8 CRLB values, validation/mesh convergence numbers, dataset counts).
- **Found wrong and fixed:** four figure panels contained hard-coded/fabricated numbers (Fig 1 b–d, Fig 4 b, Fig 6 b–d, Fig 3 d constant); manuscript numbers 0.578→0.591, −4.3%→−17.4%, Fig 6 norms/partition/trace updated; Data Availability + supplementary references removed per rules; abstract shortened to 188 words; declarations added.
- **Open:** PDF rebuild (no TeX in sandbox), author confirmation of declarations/AI wording/journal choice, 54 unverified references.
