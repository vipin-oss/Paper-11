"""Validate round-3 CSVs and regenerate proposed text/LaTeX; never edits main.tex.
Run after compute_crlb_corridor.py and its independent clean rerun.
"""
from pathlib import Path
import hashlib
import json
import sys
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
DATA = HERE / 'data'
SOURCE = HERE.parent / 'manuscript' / 'main.tex'


def main():
    sys.path.insert(0,str(HERE.parent/'program'))
    import compute_crlb_corridor as c
    # Synthetic full inverse and exact-resonance null-direction checks.
    synthetic=np.array([[2.,.2,.1],[.3,1.4,.7],[1.,.1,2.],[.2,.5,.8]])
    se,_,_=c.fim_quantities(synthetic)
    expected=np.sqrt(np.linalg.inv(synthetic.T@synthetic)[0,0])
    assert np.isclose(se,expected,rtol=1e-12)
    J1=c.sensitivity_jacobian(1.,.1,np.array([.05,.5,1.]),[0.,1.])
    resonance_null_relative_error=np.linalg.norm(J1[:,0]+J1[:,1])/np.linalg.norm(J1[:,0])
    assert resonance_null_relative_error<1e-10
    roots = pd.read_csv(DATA/'crlb_corridor_direct.csv')
    comp = pd.read_csv(DATA/'table5_comparison.csv')
    scan = pd.read_csv(DATA/'crlb_B_scan.csv')
    fd = pd.read_csv(DATA/'fd_crosscheck.csv')
    crossings = pd.read_csv(DATA/'all_crossings.csv')
    assert len(comp)==15 and len(roots)==60 and len(scan)==660
    assert len(roots.drop_duplicates(['configuration','Nt','noise_percent','criterion','side']))==60
    assert (roots.method=='direct full-FIM').all()
    assert (comp.old_method=='local asymptotic estimate').all()
    assert (scan.B.between(.0009999,1.9990001) & (scan.B!=1)).all()
    assert np.allclose(scan.relative_SE_tau_q,scan.sigma_n*scan.SE_tau_q_per_sigma_n/.02)
    assert np.allclose(scan.singular_relative_bound,scan.sigma_n/(scan.sigma_min*.02))
    assert (scan.singular_relative_bound >= scan.relative_SE_tau_q*(1-1e-12)).all()
    attained=roots[roots.status=='attained']
    missing=roots[roots.status=='not attained in scan']
    assert len(attained)+len(missing)==60
    assert attained.residual.abs().max()<2e-6
    assert missing.B_crit.isna().all() and missing.abs_B_minus_1_crit.isna().all()
    assert np.allclose(attained.abs_B_minus_1_crit,abs(attained.B_crit-1))
    assert crossings.residual.abs().max()<2e-6
    for (cfg,nt,noise,side), group in roots.groupby(['configuration','Nt','noise_percent','side']):
        p=group[group.criterion=='param'].iloc[0]
        s=group[group.criterion=='singular'].iloc[0]
        if s.status=='attained':
            assert p.status=='attained'
            assert s.abs_B_minus_1_crit >= p.abs_B_minus_1_crit-1e-6
    for key, group in roots.groupby(['configuration','Nt','side','criterion']):
        ordered=group.sort_values('noise_percent')
        ds=ordered.abs_B_minus_1_crit.fillna(np.inf).to_numpy()
        assert np.all(np.diff(ds[np.isfinite(ds)])>=-1e-6),key
    assert fd.rel_diff_SE.max()<1e-4 and fd.rel_diff_sigma_min.max()<1e-4
    assert fd.analytic_vs_mpdiff_rel_J.max()<1e-7
    assert len(fd)==4 and fd.configuration.nunique()==2
    assert len(pd.read_csv(DATA/'fig10a_direct.csv'))==18

    tex=SOURCE.read_text()
    lines=tex.splitlines()
    # Assert old table values were transcribed from manuscript, not inferred.
    old_rows=[line for line in lines if (' & ' in line and '\\%' in line and
               ('Rear-face only' in line or 'Front + Rear' in line or 'Multi-location' in line))]
    assert len(old_rows)==15
    for raw,(_,row) in zip(old_rows,comp.iterrows()):
        fields=[v.strip().replace('\\\\','').strip() for v in raw.split('&')]
        assert int(fields[1])==row.Nt
        assert float(fields[2].replace('\\%',''))==row.noise_percent
        assert float(fields[-2])==row.old_param
        assert float(fields[-1])==row.old_singular

    def value(cfg,nt,pct,criterion,side):
        r=roots[(roots.configuration==cfg)&(roots.Nt==nt)&(roots.noise_percent==pct)&
                (roots.criterion==criterion)&(roots.side==side)].iloc[0]
        return 'not attained in scan' if r.status!='attained' else f'{r.abs_B_minus_1_crit:.4f}'
    rear='Rear-face only (x=1.0)'
    front='Front + Rear (x=0, 1)'
    def pair(cfg,nt,pct,criterion='param'):
        return value(cfg,nt,pct,criterion,'B<1')+' / '+value(cfg,nt,pct,criterion,'B>1')
    labels={rear:r'Rear-face ($x=1$)',front:r'Front + rear ($x=0,1$)',
            'Multi-location (5 stations)':'Five stations'}
    rows=[]
    for _,r in comp.iterrows():
        values=[]
        for criterion in ('param','singular'):
            for side in ('B<1','B>1'):
                v=value(r.configuration,r.Nt,r.noise_percent,criterion,side)
                values.append('NA' if v=='not attained in scan' else v)
        rows.append(labels[r.configuration]+' & '+str(r.Nt)+' & '+f'{r.noise_percent:.1f}'+
                    r'\% & '+' & '.join(values)+r' \\')
    header=r'Configuration & $N_t$ & Noise & Param $B<1$ & Param $B>1$ & Singular $B<1$ & Singular $B>1$ \\'
    (HERE/'table5_new_rows.tex').write_text('% All entries: direct full-FIM. NA = not attained in scan.\n'+
        '% Entries are distances |B-1|, not the B coordinates.\n'+header+'\n\\midrule\n'+'\n'.join(rows)+'\n')
    markdown=['| Configuration | Nt | Noise % | Old Param | Old Singular | Direct Param L | Direct Param R | Direct Singular L | Direct Singular R |',
              '|---|---:|---:|---:|---:|---:|---:|---:|---:|']
    for _,r in comp.iterrows():
        vals=[value(r.configuration,r.Nt,r.noise_percent,c,s).replace('not attained in scan','NA')
              for c in ('param','singular') for s in ('B<1','B>1')]
        markdown.append('| '+' | '.join([r.configuration,str(r.Nt),str(r.noise_percent),
            f'{r.old_param:.3f}',f'{r.old_singular:.3f}']+vals)+' |')

    new={}
    new[51]=lines[50].replace('quantifies the practical precision corridor.',
        'evaluates branch-specific practical precision thresholds by direct full-FIM root-finding, separately from local near-resonance asymptotic estimates.')
    new[168]=(r'    \item Numerical evaluation of the two-sided linear asymptotic scaling law '
        r'$\sigma_{\min}\approx C(Bi)|B-1|$ near resonance, supplemented by direct full-FIM '
        r'10\%-precision threshold searches for the 15 configuration--noise combinations in Table~\ref{tab:crlb_analysis}, '
        r'with separate $B<1$ and $B>1$ searches. The original evaluations at $B=0.99$ provide local asymptotic estimates, not global corridor boundaries.')
    new[632]=(r'For the three independent parameters $\bm p=(\tau_q,\kappa^2,Bi)^\top$, we compute '
        r'$\bm F=\sigma_n^{-2}\bm J^\top\bm J$ using analytic Laplace-domain sensitivities, '
        r'with $Bi_0=Bi_L=0.10$, $\tau_{q,0}=0.02$, $\kappa^2=B\tau_{q,0}$, '
        r'$\tau_\Delta=0.04$, and $t=\operatorname{linspace}(0.05,1,N_t)$. '
        r'Table~\ref{tab:crlb_analysis} reports direct roots of '
        r'$\sigma_n\sqrt{[(\bm J^\top\bm J)^{-1}]_{11}}/\tau_{q,0}=0.10$ '
        r'and, separately, $\sigma_n/(\sigma_{\min}(\bm J)\tau_{q,0})=0.10$. '
        r'The latter is a worst-direction bound in the unrescaled coordinates, not a marginal $\tau_q$ bound. '
        r'We bracket each criterion on 22 points per side, $B\in[0.001,0.9999]$ and '
        r'$B\in[1.0001,1.999]$, and refine crossings by Brent root-finding. '
        r'The closest crossing to resonance on each branch defines its reported distance $|B-1|$. '
        r'NA means not attained in scan; no boundary is extrapolated. At $B=1$, $\tau_q$ is not estimable.')
    new[636]=(r'\caption{Direct full-FIM branch-specific threshold distances $|B-1|$ for '
        r'the Table~\ref{tab:crlb_analysis} configurations. Param: marginal '
        r'$\operatorname{SE}(\hat\tau_q)/\tau_{q,0}=0.10$. Singular: '
        r'$\sigma_n/(\sigma_{\min}(\bm J)\tau_{q,0})=0.10$, a distinct worst-direction criterion. '
        r'All four threshold columns are direct full-FIM results; $B<1$ and $B>1$ are evaluated separately. '
        r"NA denotes ``not attained in scan'' on $[0.001,0.9999]$ or $[1.0001,1.999]$, "
        r'not an infinite corridor. No claim is made outside the scan or between unresolved sample points. '
        r'The former single-point $B=0.99$ values are local asymptotic estimates and are retained only in the comparison table.}')
    new[641]=header
    # Explicitly distinguish corrected panel (a) from still-asymptotic panels (b,c).
    new[670]=lines[669].replace(
        r'(a) Critical near-resonance corridor half-width $|B-1|_{\mathrm{crit}}$ required for',
        r'(a) Separate direct full-FIM threshold distances on the $B<1$ and $B>1$ branches required for').replace(
        r'(b) Corridor narrowing under increasing temporal sampling density',
        r'(b) Local asymptotic corridor estimates (not direct thresholds) under increasing temporal sampling density').replace(
        r'(c) Practical precision phase diagram in the $(B, \sigma_n)$ plane demarcating the 10\%-precision-failure corridor from the 10\%-precision regime for rear-face monitoring; the corridor boundary is computed from the three-parameter CRLB analysis.',
        r'(c) Local asymptotic phase diagram from linear extrapolation at $B=0.99$; these plotted boundaries must not be interpreted as directly verified CRLB thresholds.')
    new[676]=(r'    \item \textit{Coarse sampling:} For rear-face $N_t=10$, direct marginal threshold distances '
        f'are {pair(rear,10,.1)} at 0.1\\% noise, {pair(rear,10,.5)} at 0.5\\%, and '
        f'{pair(rear,10,1.)} at 1.0\\% (left / right branches). The distinct singular-value '
        f'criterion gives {pair(rear,10,.1,"singular")} at 0.1\\% and '
        f'{pair(rear,10,1.,"singular")} at 1.0\\%. The earlier 0.35 and 3.48 values were local asymptotic estimates, not verified thresholds.')
    new[677]=(r'    \item \textit{Rear-face sampling ($N_t=50$):} The marginal CRLB threshold distances are '
        f'{pair(rear,50,.1)} under 0.1\\% noise and {pair(rear,50,1.)} under 1.0\\% '
        r'(left / right). A symmetric exclusion interval $[0.083,1.917]$ is not justified by these direct calculations.')
    new[678]=(r'    \item \textit{Temporal sampling:} Table~\ref{tab:crlb_analysis} gives direct branch-specific '
        r'results for rear-face $N_t=10,50,100$. The $1/\sqrt{N_t}$ scaling is only an asymptotic heuristic; '
        r'it is not imposed on these roots, whose values depend on the actual temporal grid. '
        r'Figure~\ref{fig:fig10}(b), unless regenerated, remains a plot of local asymptotic estimates.')
    new[679]=(r'    \item \textit{Dual-surface monitoring ($x=0,1$):} With $N_t=50$, '
        f'the direct marginal threshold distances are {pair(front,50,.1)} at 0.1\\% noise '
        f'and {pair(front,50,1.)} at 1.0\\% (left / right). '
        r'The prefactor $\sigma_3/|B-1|\approx23.6$ evaluated at $B=0.99$ describes local scaling, not an exact threshold relation.')
    old_tail=(r'amplifying the non-degenerate singular values and contracting the required noise corridor from '
        r'$|B-1| \ge 0.917$ down to $|B-1| \ge 0.151$ under $1.0\%$ noise.')
    new_tail=(r'improving the direct marginal CRLB thresholds under $1.0\%$ noise ($N_t=50$): '
        f'the rear-face threshold distance is {value(rear,50,1.,"param","B<1")} on the $B<1$ branch '
        f'(not attained in scan on $B>1$), compared with {pair(front,50,1.)} '
        r'for front-plus-rear monitoring (left / right). These are branch-specific crossings, not a symmetric corridor extrapolated from resonance.')
    assert old_tail in lines[697]
    new[698]=lines[697].replace(old_tail,new_tail)
    new[724]=(r'    \item For $10\%$ precision on $\tau_q$ under rear-face sampling ($N_t=50$), '
        f'direct full-FIM threshold distances are {pair(rear,50,.1)} at 0.1\\% noise '
        f'and {pair(rear,50,1.)} at 1.0\\% (left / right). Front-plus-rear monitoring gives '
        f'{pair(front,50,.1)} at 0.1\\% and {pair(front,50,1.)} at 1.0\\%. '
        r"``Not attained in scan'' refers only to the stated sampled branch; no extrapolation is used.")
    for line,replacement in new.items():
        assert lines[line-1]!=replacement, line

    doc='''# Round 3: proposed CRLB manuscript corrections

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

'''+ '\n'.join(markdown)+'''

## Table 5: LaTeX-ready replacement

Use the existing seven-column tabular structure. All numeric cells below are direct full-FIM.
Define `NA` in the caption as "not attained in scan"; four decimal places are displayed.
Source: `data/table5_comparison.csv`; separately saved as `table5_new_rows.tex`.

```latex
'''+(HERE/'table5_new_rows.tex').read_text()+'''
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

'''
    for line,replacement in new.items():
        doc+=f'### main.tex:{line}\n\nOLD\n```latex\n{lines[line-1]}\n```\n\nNEW\n```latex\n{replacement}\n```\n\n'
    digest=hashlib.sha256(SOURCE.read_bytes()).hexdigest()
    doc+=f'\nSource SHA-256: `{digest}`.\n\nOther near-resonance scaling statements (e.g. lines 613, 621, 626, 723) describe local asymptotics rather than Table 5 thresholds and are not changed by this audit.\n'
    (HERE/'crlb_text_changes.md').write_text(doc)
    result=dict(table5_rows=len(comp),threshold_entries=len(roots),scan_rows=len(scan),
        attained=len(attained),not_attained=len(missing),bracketed_crossings=len(crossings),
        max_root_relative_error_residual=float(attained.residual.abs().max()),
        max_FD_relative_SE_difference=float(fd.rel_diff_SE.max()),
        max_FD_relative_sigma_min_difference=float(fd.rel_diff_sigma_min.max()),
        max_analytic_vs_mpdiff_relative_J_difference=float(fd.analytic_vs_mpdiff_rel_J.max()),
        main_tex_sha256=digest,
        exact_resonance_null_relative_error=float(resonance_null_relative_error))
    (HERE/'validation_summary.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
    print('PASS: coverage, formulas, root residuals, singular >= marginal, FD, source-table transcription, and document generation')

if __name__=='__main__':
    main()
