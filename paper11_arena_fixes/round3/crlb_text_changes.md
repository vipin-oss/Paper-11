# Round 3: proposed CRLB manuscript corrections

**main.tex is unchanged.** Line numbers below refer to
`paper11_arena_fixes/manuscript/main.tex` at the recorded SHA-256.
These changes are proposals, not an assertion that the existing PDF or figure has been corrected.

## Definitions and scope

- Full three-parameter FIM, independent homoscedastic Gaussian noise; no parameter is fixed in the inverse.
- Marginal and singular-value criteria are distinct. The singular bound is coordinate-dependent,
  normalized using tau_q0=0.02, and is not SE(tau_q).
- Five geometry/time-grid combinations times three noise levels = 15 Table 5 rows.
  Each has four branch/criterion thresholds (60 entries).
- Scan: 22 B values per branch, [0.001,0.9999] and [1.0001,1.999], logarithmically dense near 1.
  Exact B=1 is structurally singular and excluded. All sign-change brackets are refined directly.
- A root is a crossing, not by itself a proof that every more distant B meets the criterion.
  `all_crossings.csv` retains every bracketed crossing; the table uses the closest to resonance.
- NA = **not attained in scan**, not infinity and not a result for all physical B>=0.
  Finite grids cannot rule out unresolved crossings between samples. No values are extrapolated.
- Old estimates (both coordinate and singular) are labelled **local asymptotic estimate**.
  New entries, including NA statuses, are labelled **direct full-FIM**.

## Table 5: comparison (distances |B-1|)

L/R mean B<1/B>1. Old columns are local asymptotic estimates, not direct thresholds.

| Configuration | Nt | Noise % | Old Param | Old Singular | Direct Param L | Direct Param R | Direct Singular L | Direct Singular R |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Rear-face only (x=1.0) | 10 | 0.1 | 0.246 | 0.348 | 0.2109 | 0.2486 | 0.2906 | 0.3490 |
| Rear-face only (x=1.0) | 10 | 0.5 | 1.228 | 1.739 | 0.5533 | NA | 0.7541 | NA |
| Rear-face only (x=1.0) | 10 | 1.0 | 2.455 | 3.478 | 0.6853 | NA | NA | NA |
| Rear-face only (x=1.0) | 50 | 0.1 | 0.092 | 0.130 | 0.0863 | 0.0991 | 0.1178 | 0.1461 |
| Rear-face only (x=1.0) | 50 | 0.5 | 0.459 | 0.648 | 0.3362 | 0.6618 | 0.4142 | NA |
| Rear-face only (x=1.0) | 50 | 1.0 | 0.917 | 1.296 | 0.5190 | NA | 0.5962 | NA |
| Rear-face only (x=1.0) | 100 | 0.1 | 0.065 | 0.091 | 0.0620 | 0.0683 | 0.0854 | 0.0994 |
| Rear-face only (x=1.0) | 100 | 0.5 | 0.323 | 0.456 | 0.2582 | 0.4195 | 0.3283 | 0.7094 |
| Rear-face only (x=1.0) | 100 | 1.0 | 0.645 | 0.912 | 0.4229 | NA | 0.5037 | NA |
| Front + Rear (x=0, 1) | 50 | 0.1 | 0.015 | 0.021 | 0.0150 | 0.0151 | 0.0211 | 0.0215 |
| Front + Rear (x=0, 1) | 50 | 0.5 | 0.075 | 0.106 | 0.0749 | 0.0762 | 0.1017 | 0.1132 |
| Front + Rear (x=0, 1) | 50 | 1.0 | 0.151 | 0.212 | 0.1499 | 0.1555 | 0.1965 | 0.2446 |
| Multi-location (5 stations) | 50 | 0.1 | 0.013 | 0.018 | 0.0126 | 0.0127 | 0.0176 | 0.0181 |
| Multi-location (5 stations) | 50 | 0.5 | 0.063 | 0.089 | 0.0615 | 0.0652 | 0.0839 | 0.0962 |
| Multi-location (5 stations) | 50 | 1.0 | 0.126 | 0.178 | 0.1199 | 0.1350 | 0.1585 | 0.2098 |

## Table 5: LaTeX-ready replacement

Use the existing seven-column tabular structure. All numeric cells below are direct full-FIM.
Define `NA` in the caption as "not attained in scan"; four decimal places are displayed.
Source: `data/table5_comparison.csv`; separately saved as `table5_new_rows.tex`.

```latex
% All entries: direct full-FIM. NA = not attained in scan.
% Entries are distances |B-1|, not the B coordinates.
Configuration & $N_t$ & Noise & Param $B<1$ & Param $B>1$ & Singular $B<1$ & Singular $B>1$ \\
\midrule
Rear-face ($x=1$) & 10 & 0.1\% & 0.2109 & 0.2486 & 0.2906 & 0.3490 \\
Rear-face ($x=1$) & 10 & 0.5\% & 0.5533 & NA & 0.7541 & NA \\
Rear-face ($x=1$) & 10 & 1.0\% & 0.6853 & NA & NA & NA \\
Rear-face ($x=1$) & 50 & 0.1\% & 0.0863 & 0.0991 & 0.1178 & 0.1461 \\
Rear-face ($x=1$) & 50 & 0.5\% & 0.3362 & 0.6618 & 0.4142 & NA \\
Rear-face ($x=1$) & 50 & 1.0\% & 0.5190 & NA & 0.5962 & NA \\
Rear-face ($x=1$) & 100 & 0.1\% & 0.0620 & 0.0683 & 0.0854 & 0.0994 \\
Rear-face ($x=1$) & 100 & 0.5\% & 0.2582 & 0.4195 & 0.3283 & 0.7094 \\
Rear-face ($x=1$) & 100 & 1.0\% & 0.4229 & NA & 0.5037 & NA \\
Front + rear ($x=0,1$) & 50 & 0.1\% & 0.0150 & 0.0151 & 0.0211 & 0.0215 \\
Front + rear ($x=0,1$) & 50 & 0.5\% & 0.0749 & 0.0762 & 0.1017 & 0.1132 \\
Front + rear ($x=0,1$) & 50 & 1.0\% & 0.1499 & 0.1555 & 0.1965 & 0.2446 \\
Five stations & 50 & 0.1\% & 0.0126 & 0.0127 & 0.0176 & 0.0181 \\
Five stations & 50 & 0.5\% & 0.0615 & 0.0652 & 0.0839 & 0.0962 \\
Five stations & 50 & 1.0\% & 0.1199 & 0.1350 & 0.1585 & 0.2098 \\

```

## Fig. 10(a)

Use `data/fig10a_direct.csv` (18 rows: three Nt=50 geometries, three noise levels,
two branches, marginal criterion only). Plot branches separately; omit NA markers rather
than substituting 0, 1, or extrapolated values. `B_crit` gives the actual B coordinate;
`abs_B_minus_1_crit` gives its distance from resonance. The existing figure has NOT been redrawn.
Panels (b,c) also use the old extrapolation. They require regeneration or explicit labelling
as local asymptotic estimates; the proposed caption below takes the labelling option.

## Exact line-wise old -> proposed new

The table body (original lines 643-661, 15 data rows plus separators) is replaced by the
LaTeX block above. Header and caption changes are listed separately below.
The abstract contains no numeric corridor widths, but its method description is clarified.

### main.tex:51

OLD
```latex
This paper studies whether boundary heat loss can remove the Guyer--Krumhansl (GK) parameter-identifiability singularity that arises in laser flash testing at the Fourier-resonance condition $B \equiv \kappa^2 / (\alpha_0 \tau_q) = 1$, where the sensitivity collinearity $\partial \theta / \partial \tau_q = -\partial \theta / \partial \kappa^2$ makes the Fisher Information Matrix (FIM) structurally singular. We derive the closed-form transfer function for the 1D GK system with asymmetric Robin boundary heat loss and prove (Theorem~1) that at $B = 1$ both the bulk propagation factor and the dynamic surface thermal conductivity collapse to their classical Fourier counterparts, so the GK temperature field coincides identically with the Fourier solution for all positions, all times, and all non-negative Biot numbers. The collinearity, the FIM null direction, and $\sigma_3 \equiv 0$ are therefore invariant under linear Robin cooling. Near resonance, the smallest singular value scales two-sided linearly as $\sigma_{\min} \propto |B-1|$, and a three-parameter Cramér--Rao bound analysis for $\bm{p} = (\tau_q, \kappa^2, Bi)^\top$ quantifies the practical precision corridor. Boundary dissipation and multi-sensor detection improve off-resonance precision but cannot resolve the exact singularity within the linear model class, with direct consequences for heat-loss-corrected flash data reduction.
```

NEW
```latex
This paper studies whether boundary heat loss can remove the Guyer--Krumhansl (GK) parameter-identifiability singularity that arises in laser flash testing at the Fourier-resonance condition $B \equiv \kappa^2 / (\alpha_0 \tau_q) = 1$, where the sensitivity collinearity $\partial \theta / \partial \tau_q = -\partial \theta / \partial \kappa^2$ makes the Fisher Information Matrix (FIM) structurally singular. We derive the closed-form transfer function for the 1D GK system with asymmetric Robin boundary heat loss and prove (Theorem~1) that at $B = 1$ both the bulk propagation factor and the dynamic surface thermal conductivity collapse to their classical Fourier counterparts, so the GK temperature field coincides identically with the Fourier solution for all positions, all times, and all non-negative Biot numbers. The collinearity, the FIM null direction, and $\sigma_3 \equiv 0$ are therefore invariant under linear Robin cooling. Near resonance, the smallest singular value scales two-sided linearly as $\sigma_{\min} \propto |B-1|$, and a three-parameter Cramér--Rao bound analysis for $\bm{p} = (\tau_q, \kappa^2, Bi)^\top$ evaluates branch-specific practical precision thresholds by direct full-FIM root-finding, separately from local near-resonance asymptotic estimates. Boundary dissipation and multi-sensor detection improve off-resonance precision but cannot resolve the exact singularity within the linear model class, with direct consequences for heat-loss-corrected flash data reduction.
```

### main.tex:168

OLD
```latex
    \item Numerical evaluation of the two-sided linear asymptotic scaling law $\sigma_{\min} \approx C(Bi) |B - 1|^{1.0}$ near resonance, alongside a systematic Cramér--Rao Lower Bound uncertainty analysis across 27 measurement configurations at $B = 0.99$ and 27 at $B = 1$.
```

NEW
```latex
    \item Numerical evaluation of the two-sided linear asymptotic scaling law $\sigma_{\min}\approx C(Bi)|B-1|$ near resonance, supplemented by direct full-FIM 10\%-precision threshold searches for the 15 configuration--noise combinations in Table~\ref{tab:crlb_analysis}, with separate $B<1$ and $B>1$ searches. The original evaluations at $B=0.99$ provide local asymptotic estimates, not global corridor boundaries.
```

### main.tex:632

OLD
```latex
To assess practical parameter estimation limits, we evaluate the Cramér--Rao lower bounds in the full three-parameter estimation space $\bm{p} = (\tau_q, \kappa^2, Bi)^\top$, using the sensitivity Jacobian $\bm{J} = [\bm{J}_{\tau_q} \; \bm{J}_{\kappa^2} \; \bm{J}_{Bi}]$ and the Fisher matrix $\bm{F} = \sigma_n^{-2} \bm{J}^\top \bm{J}$ of Eq.~\eqref{eq:fim_def}; every Jacobian column is obtained by analytic differentiation of the closed-form transfer function. Near resonance ($B = 0.99$) the Fisher matrix is non-singular and the ordinary CRLB applies, while at exact resonance ($B = 1$) it is structurally singular and the generalized (Moore--Penrose) bound of Eq.~\eqref{eq:generalized_crlb} applies. Table~\ref{tab:crlb_analysis} and Figure~\ref{fig:fig10} summarize the critical corridor half-width $|B - 1|_{\mathrm{crit}}$ needed to achieve a $10\%$ relative parameter standard error ($\operatorname{SE}(\hat{\tau}_q) / \tau_{q,0} \le 0.10$ with baseline $\tau_{q,0} = 0.02$) across measurement geometries, temporal sampling sizes $N_t$, and noise levels $\sigma_n \in \{0.1\%, 0.5\%, 1.0\%\}$. All values quoted below are listed in Table~\ref{tab:crlb_analysis}.
```

NEW
```latex
For the three independent parameters $\bm p=(\tau_q,\kappa^2,Bi)^\top$, we compute $\bm F=\sigma_n^{-2}\bm J^\top\bm J$ using analytic Laplace-domain sensitivities, with $Bi_0=Bi_L=0.10$, $\tau_{q,0}=0.02$, $\kappa^2=B\tau_{q,0}$, $\tau_\Delta=0.04$, and $t=\operatorname{linspace}(0.05,1,N_t)$. Table~\ref{tab:crlb_analysis} reports direct roots of $\sigma_n\sqrt{[(\bm J^\top\bm J)^{-1}]_{11}}/\tau_{q,0}=0.10$ and, separately, $\sigma_n/(\sigma_{\min}(\bm J)\tau_{q,0})=0.10$. The latter is a worst-direction bound in the unrescaled coordinates, not a marginal $\tau_q$ bound. We bracket each criterion on 22 points per side, $B\in[0.001,0.9999]$ and $B\in[1.0001,1.999]$, and refine crossings by Brent root-finding. The closest crossing to resonance on each branch defines its reported distance $|B-1|$. NA means not attained in scan; no boundary is extrapolated. At $B=1$, $\tau_q$ is not estimable.
```

### main.tex:636

OLD
```latex
\caption{Three-parameter Cramér--Rao Lower Bound (CRLB) estimation precision: critical near-resonance corridor half-width $|B - 1|_{\mathrm{crit}}$ required to achieve $\operatorname{SE}(\hat{\tau}_q) / \tau_{q,0} \le 0.10$ across observation configurations and noise floors, evaluated at $B = 0.99$ and $Bi = 0.10$ with the Fisher matrix $\bm{F} = \sigma_n^{-2}\bm{J}^\top\bm{J}$ of the three-parameter Jacobian $\bm{J} = [\bm{J}_{\tau_q}, \bm{J}_{\kappa^2}, \bm{J}_{Bi}]$. The "Singular" column uses the smallest singular value $\sigma_3$; the "Param" column uses the marginal $\tau_q$ coordinate bound. At exact resonance ($B = 1$) the Fisher matrix is singular and only $Bi$ admits a finite generalized bound ($\operatorname{SE}(\hat{Bi})/\sigma_n = 0.1243$ for rear-face sampling with $N_t = 50$).}
```

NEW
```latex
\caption{Direct full-FIM branch-specific threshold distances $|B-1|$ for the Table~\ref{tab:crlb_analysis} configurations. Param: marginal $\operatorname{SE}(\hat\tau_q)/\tau_{q,0}=0.10$. Singular: $\sigma_n/(\sigma_{\min}(\bm J)\tau_{q,0})=0.10$, a distinct worst-direction criterion. All four threshold columns are direct full-FIM results; $B<1$ and $B>1$ are evaluated separately. NA denotes ``not attained in scan'' on $[0.001,0.9999]$ or $[1.0001,1.999]$, not an infinite corridor. No claim is made outside the scan or between unresolved sample points. The former single-point $B=0.99$ values are local asymptotic estimates and are retained only in the comparison table.}
```

### main.tex:641

OLD
```latex
Sensor Configuration & Sampling $N_t$ & Noise Level $\sigma_n$ & $\sigma_3 / |B-1|$ & $\operatorname{SE}(\hat{\tau}_q)/\sigma_n \cdot |B-1|$ & $|B-1|_{\mathrm{crit}}$ (Param) & $|B-1|_{\mathrm{crit}}$ (Singular) \\
```

NEW
```latex
Configuration & $N_t$ & Noise & Param $B<1$ & Param $B>1$ & Singular $B<1$ & Singular $B>1$ \\
```

### main.tex:670

OLD
```latex
\caption{Cramér--Rao Lower Bound (CRLB) identifiability maps and practical precision phase diagrams: (a) Critical near-resonance corridor half-width $|B-1|_{\mathrm{crit}}$ required for $\operatorname{SE}(\hat{\tau}_q)/\tau_{q,0} \le 10\%$ versus measurement noise floor $\sigma_n \in [0.1\%, 1.0\%]$ ($N_t = 50$). (b) Corridor narrowing under increasing temporal sampling density $N_t \in [10, 100]$ at $\sigma_n = 0.5\%$. (c) Practical precision phase diagram in the $(B, \sigma_n)$ plane demarcating the 10\%-precision-failure corridor from the 10\%-precision regime for rear-face monitoring; the corridor boundary is computed from the three-parameter CRLB analysis. (d) Divergence of the non-Fourier parameter uncertainty ratio $\operatorname{SE}(\hat{\tau}_q)/\operatorname{SE}(\hat{Bi})$ as resonance is approached: markers are the numerical ratio computed directly from the three-parameter Fisher matrix (rear face, $N_t = 50$, $Bi = 0.10$), and the dashed curve is the asymptotic $1/|B-1|$ estimate fitted to the well-conditioned part of the same scan. At $B = 1$ the Fisher matrix is singular and $\hat{\tau}_q$ is not estimable.}
```

NEW
```latex
\caption{Cramér--Rao Lower Bound (CRLB) identifiability maps and practical precision phase diagrams: (a) Separate direct full-FIM threshold distances on the $B<1$ and $B>1$ branches required for $\operatorname{SE}(\hat{\tau}_q)/\tau_{q,0} \le 10\%$ versus measurement noise floor $\sigma_n \in [0.1\%, 1.0\%]$ ($N_t = 50$). (b) Local asymptotic corridor estimates (not direct thresholds) under increasing temporal sampling density $N_t \in [10, 100]$ at $\sigma_n = 0.5\%$. (c) Local asymptotic phase diagram from linear extrapolation at $B=0.99$; these plotted boundaries must not be interpreted as directly verified CRLB thresholds. (d) Divergence of the non-Fourier parameter uncertainty ratio $\operatorname{SE}(\hat{\tau}_q)/\operatorname{SE}(\hat{Bi})$ as resonance is approached: markers are the numerical ratio computed directly from the three-parameter Fisher matrix (rear face, $N_t = 50$, $Bi = 0.10$), and the dashed curve is the asymptotic $1/|B-1|$ estimate fitted to the well-conditioned part of the same scan. At $B = 1$ the Fisher matrix is singular and $\hat{\tau}_q$ is not estimable.}
```

### main.tex:676

OLD
```latex
    \item \textit{Origin of Coarse-Sampling Thresholds:} The values $|B - 1|_{\mathrm{crit}} \approx 0.35$ (0.1\% noise) and $3.48$ (1.0\% noise) stem from a sparse 10-point rear-face sampling evaluating parameter uncertainty along the principal singular mode ($\sigma_n / (0.10 C \tau_{q,0})$).
```

NEW
```latex
    \item \textit{Coarse sampling:} For rear-face $N_t=10$, direct marginal threshold distances are 0.2109 / 0.2486 at 0.1\% noise, 0.5533 / not attained in scan at 0.5\%, and 0.6853 / not attained in scan at 1.0\% (left / right branches). The distinct singular-value criterion gives 0.2906 / 0.3490 at 0.1\% and not attained in scan / not attained in scan at 1.0\%. The earlier 0.35 and 3.48 values were local asymptotic estimates, not verified thresholds.
```

### main.tex:677

OLD
```latex
    \item \textit{Influence of Sampling Density ($N_t = 50$):} With denser sampling ($N_t = 50$, Figure~\ref{fig:fig10}(a)), the sensitivity norm increases roughly as $\sqrt{N_t}$, contracting the critical corridor for parameter $\tau_q$ to $|B - 1|_{\mathrm{crit}} \approx 0.092$ under $0.1\%$ noise and $0.917$ under $1.0\%$ noise. At $1.0\%$ noise, simultaneous parameter estimation remains challenging across $B \in [0.083, 1.917]$ (Figure~\ref{fig:fig10}(c)).
```

NEW
```latex
    \item \textit{Rear-face sampling ($N_t=50$):} The marginal CRLB threshold distances are 0.0863 / 0.0991 under 0.1\% noise and 0.5190 / not attained in scan under 1.0\% (left / right). A symmetric exclusion interval $[0.083,1.917]$ is not justified by these direct calculations.
```

### main.tex:678

OLD
```latex
    \item \textit{Temporal Sampling Progression ($N_t$):} Figure~\ref{fig:fig10}(b) illustrates the narrowing of the corridor as $N_t$ rises from 10 to 100, following an asymptotic $\sim 1/\sqrt{N_t}$ contraction.
```

NEW
```latex
    \item \textit{Temporal sampling:} Table~\ref{tab:crlb_analysis} gives direct branch-specific results for rear-face $N_t=10,50,100$. The $1/\sqrt{N_t}$ scaling is only an asymptotic heuristic; it is not imposed on these roots, whose values depend on the actual temporal grid. Figure~\ref{fig:fig10}(b), unless regenerated, remains a plot of local asymptotic estimates.
```

### main.tex:679

OLD
```latex
    \item \textit{Dual-Surface Monitoring ($x = 0$ and $x = 1$):} Recording both front and rear surfaces substantially increases the sensitivity norm ($\sigma_3/|B-1| \approx 23.6$), reducing the corridor to $|B - 1|_{\mathrm{crit}} \approx 0.015$ (0.1\% noise) and $0.151$ (1.0\% noise).
```

NEW
```latex
    \item \textit{Dual-surface monitoring ($x=0,1$):} With $N_t=50$, the direct marginal threshold distances are 0.0150 / 0.0151 at 0.1\% noise and 0.1499 / 0.1555 at 1.0\% (left / right). The prefactor $\sigma_3/|B-1|\approx23.6$ evaluated at $B=0.99$ describes local scaling, not an exact threshold relation.
```

### main.tex:698

OLD
```latex
In photothermal and flash metrology, multi-sensor configurations (such as simultaneous front-face thermoreflectance and rear-face pyrometry~\cite{salazar2014extending}) are increasingly employed to decouple surface from bulk phenomena. As demonstrated by Emery et al.~\cite{emery2000uncertainties}, Fadale et al.~\cite{fadale1995two}, Dowding et al.~\cite{dowding1995estimation}, and {\L}apka et al.~\cite{lapka2020development}, optimal sensor placement maximizes the determinant of the Fisher Information Matrix. From Theorem~1, if $\theta_{\mathrm{GK}}(x, t; B=1) \equiv \theta_{\mathrm{Fourier}}(x, t)$ holds for all $x \in [0, 1]$ and all $t \ge 0$, then differentiating with respect to the constitutive parameters yields $\partial \theta / \partial \tau_q(x, t) = -\partial \theta / \partial \kappa^2(x, t)$ at every spatial point across the specimen (dimensionless parameterization). Stacking measurements from multiple sensors adds additional rows to the sensitivity matrix $\bm{J}$, but the two parameter columns remain strictly linearly dependent: $\bm{J}_{\tau_q} = -\bm{J}_{\kappa^2}$. Consequently, at exact resonance ($B = 1$) under the linear 1D model, multi-sensor monitoring cannot remove the structural rank deficiency ($\operatorname{rank}(\bm{J}) \le 2$ in the 3-parameter space $(\tau_q, \kappa^2, Bi)$). However, for near-resonant media ($B \ne 1$), our CRLB analysis (Table~\ref{tab:crlb_analysis}) demonstrates that multi-sensor configurations substantially improve near-resonance practical precision while leaving the exact-resonance structural null direction unchanged, amplifying the non-degenerate singular values and contracting the required noise corridor from $|B-1| \ge 0.917$ down to $|B-1| \ge 0.151$ under $1.0\%$ noise.
```

NEW
```latex
In photothermal and flash metrology, multi-sensor configurations (such as simultaneous front-face thermoreflectance and rear-face pyrometry~\cite{salazar2014extending}) are increasingly employed to decouple surface from bulk phenomena. As demonstrated by Emery et al.~\cite{emery2000uncertainties}, Fadale et al.~\cite{fadale1995two}, Dowding et al.~\cite{dowding1995estimation}, and {\L}apka et al.~\cite{lapka2020development}, optimal sensor placement maximizes the determinant of the Fisher Information Matrix. From Theorem~1, if $\theta_{\mathrm{GK}}(x, t; B=1) \equiv \theta_{\mathrm{Fourier}}(x, t)$ holds for all $x \in [0, 1]$ and all $t \ge 0$, then differentiating with respect to the constitutive parameters yields $\partial \theta / \partial \tau_q(x, t) = -\partial \theta / \partial \kappa^2(x, t)$ at every spatial point across the specimen (dimensionless parameterization). Stacking measurements from multiple sensors adds additional rows to the sensitivity matrix $\bm{J}$, but the two parameter columns remain strictly linearly dependent: $\bm{J}_{\tau_q} = -\bm{J}_{\kappa^2}$. Consequently, at exact resonance ($B = 1$) under the linear 1D model, multi-sensor monitoring cannot remove the structural rank deficiency ($\operatorname{rank}(\bm{J}) \le 2$ in the 3-parameter space $(\tau_q, \kappa^2, Bi)$). However, for near-resonant media ($B \ne 1$), our CRLB analysis (Table~\ref{tab:crlb_analysis}) demonstrates that multi-sensor configurations substantially improve near-resonance practical precision while leaving the exact-resonance structural null direction unchanged, improving the direct marginal CRLB thresholds under $1.0\%$ noise ($N_t=50$): the rear-face threshold distance is 0.5190 on the $B<1$ branch (not attained in scan on $B>1$), compared with 0.1499 / 0.1555 for front-plus-rear monitoring (left / right). These are branch-specific crossings, not a symmetric corridor extrapolated from resonance.
```

### main.tex:724

OLD
```latex
    \item For $10\%$ precision on $\tau_q$ under rear-face sampling ($N_t = 50$), the required corridor is $|B-1| \ge 0.092$ at $0.1\%$ noise and $|B-1| \ge 0.917$ at $1.0\%$ noise; multi-sensor detection contracts it to $|B-1| \ge 0.015$--$0.151$.
```

NEW
```latex
    \item For $10\%$ precision on $\tau_q$ under rear-face sampling ($N_t=50$), direct full-FIM threshold distances are 0.0863 / 0.0991 at 0.1\% noise and 0.5190 / not attained in scan at 1.0\% (left / right). Front-plus-rear monitoring gives 0.0150 / 0.0151 at 0.1\% and 0.1499 / 0.1555 at 1.0\%. ``Not attained in scan'' refers only to the stated sampled branch; no extrapolation is used.
```


Source SHA-256: `34f37143d9a1a5bd70d36b2573ac8aef67770b8d5b070bba6684f45a9ce93a56`.

Other near-resonance scaling statements (e.g. lines 613, 621, 626, 723) describe local asymptotics rather than Table 5 thresholds and are not changed by this audit.
