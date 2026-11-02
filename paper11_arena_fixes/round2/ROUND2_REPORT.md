# Arena round-2 report (2026-11-02)

All outputs in `paper11_arena_fixes/round2/`; no pre-existing file was modified.

1. Fig. 6 re-derivation — `../program/compute_fig6_data.py` run clean (EXIT=0), log in
   `fig6_rerun_log.txt`, CSVs in `data/`. Clean rerun matches shipped CSVs to max abs diff
   2.3e-13 (Nt=10) / 1.4e-12 (Nt=50). Sanity: (i) norm_J_tq == norm_J_k2 exactly at B=1;
   (ii) independent FD (dps=50, h=1e-8) at x=0 and x=1 agrees to rel diff 2.0e-08 / 5.3e-10;
   (iii) tr(F) Nt50/Nt10 = 3.11 because J is concentrated at early times (10-pt grid step 0.105
   vs 50-pt step 0.0196). Computed Fig6(a) rho(tau_q,Bi) = 0.0700 (Nt=50) / 0.0798 (Nt=10)
   vs the 0.0382 hard-coded in `plot_figures.py` — discrepancy flagged, not silently changed.
   The 34.8 vs 17.89 norm difference is explained by the time grid: user FD used t in
   [0.01,1]; this script uses linspace(0.05,1.0,Nt) as in the master `run_near_resonance.py`.
2. nishi2020effect — second author confirmed as Azuma, Naoyoshi
   (degruyterbrill.com/document/doi/10.1515/htmp-2020-0079 + author's Ibaraki Univ. list);
   see `nishi2020effect_author_confirmation.md`. No change needed.
3. Random 12-key re-check (seed 20261102, from the 58 'verified' rows of
   `references_verification.csv`) — all 12 OK field-by-field; see `recheck_12keys.csv`.
   zhu2017nonlocal published title is "Nonlocal Effects and Slip Heat Flow in Nanolayers"
   (nature.com/articles/s41598-017-10416-9); title fixed in `round2/manuscript/references.bib`.
4. `highlights.txt` — 5 Elsevier-style highlights, each <= 85 chars, from verified numbers.
5. `README.md` / `REPRODUCIBILITY.md` — copies of the master docs with run-order update:
   run `compute_fig6_data.py` (and any compute_*_data script) BEFORE `plot_figures.py`.
