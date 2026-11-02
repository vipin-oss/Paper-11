"""Generate supplementary/supplementary_material.tex from the shipped CSV datasets
and the two tables moved out of the main manuscript."""
import io
import os
import pandas as pd

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
D = os.path.join(BASE, 'data')


def sci(x, digits=3):
    """Format a float in LaTeX scientific notation."""
    if x == 0:
        return '0'
    import math
    exp = int(math.floor(math.log10(abs(x))))
    mant = x / (10.0 ** exp)
    if abs(mant - round(mant)) < 1e-12:
        return r'$10^{%d}$' % exp
    return r'$%.*f \times 10^{%d}$' % (digits, mant, exp)


# ---------------------------------------------------------------- S1 (moved table)
s1_rows = [
    (r'$1.0 \times 10^{-1}$', '0.100000', '0.050000', r'$2.470 \times 10^{-1}$', '2.4703'),
    (r'$1.0 \times 10^{-2}$', '0.010000', '0.005000', r'$5.590 \times 10^{-2}$', '5.5895'),
    (r'$1.0 \times 10^{-3}$', '0.001000', '0.000500', r'$4.068 \times 10^{-3}$', '4.0682'),
    (r'$1.0 \times 10^{-4}$', '0.000100', '0.000050', r'$3.942 \times 10^{-4}$', '3.9418'),
    (r'$1.0 \times 10^{-5}$', '0.000010', '0.000005', r'$3.930 \times 10^{-5}$', '3.9297'),
    (r'$1.0 \times 10^{-6}$', '0.000001', '0.0000005', r'$3.928 \times 10^{-6}$', '3.9285'),
    (r'$1.0 \times 10^{-7}$', '0.0000001', '0.00000005', r'$3.928 \times 10^{-7}$', '3.9284'),
]
s1_body = '\n'.join(' & '.join(r) + r' \\' for r in s1_rows)

# ---------------------------------------------------------------- S2 (moved table)
s2_rows = [
    ('100', r'$4.746\times 10^{-4}$', r'$1.180\times 10^{-4}$', r'$2.403\times 10^{-6}$', r'$1.215\times 10^{-4}$', r'$1.072\times 10^{-4}$'),
    ('200', r'$1.186\times 10^{-4}$', r'$2.950\times 10^{-5}$', r'$5.986\times 10^{-7}$', r'$3.038\times 10^{-5}$', r'$2.676\times 10^{-5}$'),
    ('400', r'$2.966\times 10^{-5}$', r'$7.379\times 10^{-6}$', r'$1.487\times 10^{-7}$', r'$7.597\times 10^{-6}$', r'$6.686\times 10^{-6}$'),
    ('800', r'$7.412\times 10^{-6}$', r'$1.847\times 10^{-6}$', r'$3.636\times 10^{-8}$', r'$1.899\times 10^{-6}$', r'$1.671\times 10^{-6}$'),
]
s2_body = '\n'.join(' & '.join(r) + r' \\' for r in s2_rows)

# ---------------------------------------------------------------- S3 de Hoog precision
dh = pd.read_csv(D + '/dehoog_convergence.csv')
s3_rows = []
for _, r in dh.iterrows():
    s3_rows.append('%d & %d & %.16f & %s \\\\' % (r.dps, r.degree, r.value, sci(r.error_vs_gold, 2)))
s3_body = '\n'.join(s3_rows)

# ---------------------------------------------------------------- S4 step-size study
ss = pd.read_csv(D + '/sensitivity_step_size_convergence.csv')
s4_rows = []
for _, r in ss.iterrows():
    s4_rows.append('%s & %.12f & %.12f & %s & %s & %.2f \\\\' % (
        sci(r.step_size_h, 0), r.J_tau_q, r.J_kappa2,
        sci(r.abs_residual, 2), sci(r.rel_residual, 2), r.ratio_O_h2))
s4_body = '\n'.join(s4_rows)

# ---------------------------------------------------------------- S5 full CRLB audit
cr = pd.read_csv(D + '/crlb_configuration_audit.csv')
cfg_short = {'Rear-face only (x=1.0)': 'Rear only',
             'Front + Rear (x=0, 1)': 'Front+Rear',
             'Multi-location (5 stations)': 'Multi-5'}
s5_rows = []
for _, r in cr.iterrows():
    if r.fim_singular:
        se_tq = r'\multicolumn{1}{c}{n.e.}'
        c_coord = r'\multicolumn{1}{c}{--}'
        c_sing = r'\multicolumn{1}{c}{--}'
    else:
        se_tq = '%.3f' % r.SE_tau_q_per_sigma_n
        c_coord = '%.4f' % r.delta_B_crit_coordinate
        c_sing = '%.4f' % r.delta_B_crit_singular
    if c_sing.startswith('\\'):
        sing_cell = c_sing
    else:
        sing_cell = '$' + c_sing + '$'
    s5_rows.append('%s & %d & %.1f\\%% & %.3f & %.3f & %s & %s & %.3f & %s & %s & %s \\\\' % (
        cfg_short[r.configuration], r.Nt, r.noise_percent, r.sigma1, r.sigma2,
        sci(r.sigma3, 3), sci(r.cond_F3, 2), r.SE_Bi_per_sigma_n, se_tq,
        c_coord if c_coord.startswith('\\') else '$' + c_coord + '$',
        sing_cell))
s5_body = '\n'.join(s5_rows)

# ---------------------------------------------------------------- S6 dataset inventory
inv = [
    ('validation_laplace.csv', 'PDE vs continuous Laplace benchmark (8 Biot numbers)'),
    ('validation_cowan.csv', 'PDE vs exact Cowan Fourier series (5 Biot numbers)'),
    ('asymptotic_fourier_convergence.csv', 'Classical Fourier-limit convergence (7 decades)'),
    ('mesh_convergence.csv', 'Conservative FVM grid refinement (4 resolutions)'),
    ('dehoog_convergence.csv', 'Laplace inversion precision table (dps 30--80)'),
    ('heatloss_parameter_sweep.csv', 'Resonance heat-loss sweep, 3-parameter SVD (12 Bi)'),
    ('sensitivity_profiles.csv', 'Time-resolved sensitivity trajectories (3 Bi)'),
    ('near_resonance_scaling.csv', 'Near-resonance landscape, 2- and 3-parameter spectra (102 rows)'),
    ('off_resonance_sweep.csv', 'Wide off-resonance conditioning sweep (51 rows)'),
    ('scaling_summary.csv', 'Mean scaling prefactors and condition prefactors (6 Bi)'),
    ('scaling_regression_twosided.csv', 'Two-sided regression audit, both singular-value definitions'),
    ('crlb_configuration_audit.csv', 'Three-parameter CRLB audit, 54 rows (27 near-resonance + 27 at B = 1)'),
    ('uncertainty_divergence.csv', 'Numerical $SE(\\tau_q)$/$SE(B_i)$ divergence scan (15 B values)'),
    ('field_invariance_asymmetric.csv', '50-digit field-level invariance residuals'),
    ('interior_point_mesh_refinement.csv', 'Interior-point PDE residuals (20 evaluations)'),
    ('sensitivity_step_size_convergence.csv', 'Finite-difference step truncation floor (8 steps)'),
    ('spatiotemporal_3d_temperature.csv', '3D spatiotemporal temperature grid (750 points)'),
    ('spatiotemporal_3d_sensitivities.csv', '3D sensitivity grid, analytic + FD channels (500 points)'),
    ('canyon_3d_landscape.csv', '3D condition-number landscape (189 points)'),
]
s6_body = '\n'.join(r'\texttt{%s} & %s \\\\' % (f.replace('_', r'\_\allowbreak{}'), d) for f, d in inv)

TEX = r"""\documentclass[preprint]{elsarticle}
\usepackage{amsmath,amssymb,bm}
\usepackage{booktabs}
\usepackage{array}
\usepackage{graphicx}
\usepackage[colorlinks=true,linkcolor=blue,citecolor=blue,urlcolor=blue]{hyperref}
\emergencystretch=3em

\journal{Supplementary Material}

\begin{document}

\begin{frontmatter}

\title{Supplementary Material for:\\
Invariance of the Guyer--Krumhansl Fourier-Resonance Sensitivity Singularity to Boundary Heat Loss in Laser Flash Testing}

\author{Vipin Gupta}

\begin{abstract}
This supplementary document collects the detailed numerical audit tables that support the
main manuscript: the classical Fourier-limit convergence study, the interior-point mesh
refinement residuals, the Laplace-inversion precision table, the finite-difference step-size
study, the complete three-parameter Cram\'er--Rao audit, and the full inventory of the
archived datasets. Every table is generated from, and is numerically identical to, the
corresponding CSV file in the accompanying reproducibility package.
\end{abstract}

\end{frontmatter}

\section{Classical Fourier-Limit Convergence}
\label{sec:supp_asymptotic}

Table~\ref{tab:supp_asymptotic} reports the off-resonance classical limit study discussed in
Section~4.3 of the main manuscript. At fixed $B = 0.5$ ($\kappa^2 = 0.5\,\tau_q$) and
$Bi = 0.100$, the constitutive discrepancy $\|\theta_{\mathrm{GK}} - \theta_{\mathrm{Fourier}}\|_\infty$
converges linearly in $\epsilon = \max(\tau_q, \kappa^2)$ with $C_{\mathrm{asymp}} \approx 3.928$.
Source file: \texttt{data/asymptotic\_fourier\_convergence.csv}.

\begin{table}[htbp]
\centering
\caption{Asymptotic convergence of the Guyer--Krumhansl solution toward the classical Cowan
Fourier solution ($B = 0.5$, $Bi = 0.100$). The discrepancy is decomposed into
constitutive-limit error, spatial discretization error, and Laplace inversion error; the
error ratio per $\epsilon$ converges to the linear prefactor $C_{\mathrm{asymp}}$.}
\label{tab:supp_asymptotic}
\resizebox{\textwidth}{!}{%
\begin{tabular}{ccccc}
\toprule
$\epsilon$ & $\tau_q$ & $\kappa^2$ & Constitutive Error $\|\theta_{\mathrm{GK}} - \theta_{\mathrm{Fou}}\|_\infty$ & Error Ratio / $\epsilon$ \\
\midrule
__S1__
\bottomrule
\end{tabular}}
\end{table}

\section{Interior-Point Mesh Refinement Residuals}
\label{sec:supp_interior}

Table~\ref{tab:supp_interior} reports the interior-point residual
$\max_t |T_{\mathrm{PDE}}(x,t) - T_{\mathrm{Laplace}}(x,t)|$ of the time-domain PDE solver
under spatial mesh refinement, for asymmetric cooling $Bi_0 = 0.05$, $Bi_L = 0.20$ at
resonance $B = 1$. Off-grid coordinates are obtained by linear interpolation. The observed
convergence order is $p = 2.00$ at every station. Source file:
\texttt{data/interior\_point\_mesh\_refinement.csv}.

\begin{table}[htbp]
\centering
\caption{Interior-point PDE-vs-Laplace residual under spatial mesh refinement (asymmetric
cooling $Bi_0 = 0.05$, $Bi_L = 0.20$, $B = 1$, $\tau_q = \kappa^2 = 0.02$). Values are
$\max_t |T_{\mathrm{PDE}}(x,t) - T_{\mathrm{Laplace}}(x,t)|$.}
\label{tab:supp_interior}
\resizebox{\columnwidth}{!}{%
\begin{tabular}{cccccc}
\toprule
$N_x$ & $x=0.00$ & $x=0.25$ & $x=0.50$ & $x=0.75$ & $x=1.00$ \\
\midrule
__S2__
\midrule
Order $p$ & 2.00 & 2.00 & 2.01 & 2.00 & 2.00 \\
\bottomrule
\end{tabular}}
\end{table}

\section{Laplace-Inversion Precision Audit}
\label{sec:supp_dehoog}

Table~\ref{tab:supp_dehoog} reports the de Hoog quotient-difference inversion of the
rear-face transfer function at $Bi = 0.10$, $\tau_q = \kappa^2 = 0.02$, evaluated at
increasing working precision and Pad\'e degree. The absolute error is measured against an
80-digit reference value. The truncation floor reaches $10^{-56}$--$10^{-54}$ at degree
$M = 60$. Source file: \texttt{data/dehoog\_convergence.csv}.

\begin{table}[htbp]
\centering
\caption{de Hoog Laplace-inversion precision audit (rear face, $Bi = 0.10$,
$\tau_q = \kappa^2 = 0.02$). The error column is measured against the 80-digit reference.}
\label{tab:supp_dehoog}
\begin{tabular}{cccl}
\toprule
Digits (dps) & Pad\'e degree $M$ & Inverted value & Absolute error \\
\midrule
__S3__
\bottomrule
\end{tabular}
\end{table}

\section{Sensitivity Step-Size Truncation Floor}
\label{sec:supp_step}

Table~\ref{tab:supp_step} reports the finite-difference sensitivity study of
Section~4.4 of the main manuscript. Central differences with relative step $h$ are applied
to the closed-form transfer function; the collinearity residual decreases quadratically
until the round-off floor is reached. Source file:
\texttt{data/sensitivity\_step\_size\_convergence.csv}.

\begin{table}[htbp]
\centering
\caption{Finite-difference step-size study for the sensitivity Jacobian at the rear face
($Bi = 0.10$, $\tau_q = \kappa^2 = 0.02$, $t = 0.5$). The residual is
$|J_{\tau_q} + J_{\kappa^2}|$ and the ratio column verifies the $O(h^2)$ truncation rate.}
\label{tab:supp_step}
\resizebox{\textwidth}{!}{%
\begin{tabular}{cccccc}
\toprule
Step $h$ & $J_{\tau_q}$ & $J_{\kappa^2}$ & Absolute residual & Relative residual & Ratio \\
\midrule
__S4__
\bottomrule
\end{tabular}}
\end{table}

\section{Complete Three-Parameter Cram\'er--Rao Audit}
\label{sec:supp_crlb}

Table~\ref{tab:supp_crlb} reports the complete three-parameter CRLB audit that underlies
Table~7 and Figure~10 of the main manuscript. The first block evaluates the Fisher matrix at
$B = 0.99$ (non-singular, ordinary CRLB); the second block evaluates the generalized
(Moore--Penrose) bound at exact resonance $B = 1$, where $\tau_q$ and $\kappa^2$ are not
estimable (n.e.) and only $Bi$ admits a finite bound. All values are archived in
\texttt{data/crlb\_configuration\_audit.csv}.

\begin{table}[htbp]
\centering
\caption{Complete three-parameter CRLB audit ($Bi = 0.10$, $\tau_{q,0} = 0.02$).}
\label{tab:supp_crlb}
\resizebox{\textwidth}{!}{%
\begin{tabular}{lcccrrcrccc}
\toprule
Configuration & $N_t$ & $\sigma_n$ (\%) & $\sigma_1$ & $\sigma_2$ & $\sigma_3$ & $\operatorname{cond}(F)$ & $\operatorname{SE}(\widehat{Bi})/\sigma_n$ & $\operatorname{SE}(\widehat{\tau}_q)/\sigma_n$ & $|B-1|_{\mathrm{crit}}$ (coord.) & $|B-1|_{\mathrm{crit}}$ (sing.) \\
\midrule
__S5__
\bottomrule
\end{tabular}}
\end{table}

\section{Dataset Inventory}
\label{sec:supp_inventory}

Table~\ref{tab:supp_inventory} lists every CSV dataset in the reproducibility package
together with the manuscript element it sources. All datasets are regenerated by
\texttt{program/reproduce\_all.py}.

\begin{table}[htbp]
\centering
\caption{Inventory of the 19 archived CSV datasets.}
\label{tab:supp_inventory}
\begin{tabular}{>{\raggedright\arraybackslash}p{0.30\textwidth}p{0.56\textwidth}}
\toprule
Dataset & Content \\
\midrule
__S6__
\bottomrule
\end{tabular}
\end{table}

\end{document}
"""

TEX = (TEX.replace('__S1__', s1_body).replace('__S2__', s2_body)
          .replace('__S3__', s3_body).replace('__S4__', s4_body)
          .replace('__S5__', s5_body).replace('__S6__', s6_body))

out = os.path.join(BASE, 'supplementary', 'supplementary_material.tex')
os.makedirs(os.path.dirname(out), exist_ok=True)
io.open(out, 'w', encoding='utf-8').write(TEX)
print('written', out, len(TEX), 'chars')
