# ARENA BRIEF — Paper 11 (kaam Arena par, verification Claude ke saath)

Input: `Paper11_Q1_SUBMISSION_READY_MASTER.zip` (isi folder me; manuscript/, program/, data/, figures/, INTERNAL_AUDIT/).
Paper: "Invariance of the Guyer–Krumhansl Fourier-Resonance Sensitivity Singularity to Boundary Heat Loss in Laser Flash Testing" (elsarticle, ~12.3k words tex, 53 cited keys, 64 bib entries).
Target journals (abhi): IJHMT, Appl. Math. Modelling, Int. J. Thermal Sciences (final choice author karega).

## Core claim (inhi ko code se verify karna hai)
At B = kappa^2/(alpha0*tau_q) = 1, GK field == classical Fourier field for ALL x, t, Bi >= 0 (Theorem 1), so dtheta/dtau_q = -dtheta/dkappa^2, FIM null direction, sigma_3 == 0 for all Bi. Near B=1: sigma_min ~ two-sided linear |B-1|. 3-parameter p=(tau_q, kappa^2, Bi) CRLB (Table 8).

## Hard rules
1. Koi data/result/reference invent nahi. Number badla to CHANGES me old->new->reason.
2. Citations sirf references.bib se (BibTeX). Manuscript me supplementary/code-archive/reproducibility package ka zikr nahi (main.tex line ~741 Data Availability check karo).
3. Document class / line spacing mat badlo. Style sample: FEM_1/2/3 (final hain, unhe mat chhuo).
4. Issues ko A fatal / B major / C minor / D optional classify karo. Sirf A, B fix. D mat karo.
5. Build: pdflatex -> bibtex -> pdflatex -> pdflatex, 0 errors / 0 undefined / 0 overfull.

## Kaam (priority order)
A. SOLVER/CODE VERIFY (sabse zaroori; byte-identical re-run kaafi nahi):
   - closed-form roots / Laplace analytic (analytic_laplace.py) vs time-domain PDE solver (solver_gk_heatloss.py) vs independent_solver.py
   - PDE residual: modal vector ko equation me daalo
   - limits: kappa^2->0 (MCV), tau_q->0 (Fourier), B=1 collapse, Bi->0 adiabatic, Bi->inf
   - convergence: grid Nx, time step, Laplace inversion (Talbot/Stehfest) params
   - sab numeric claims (counts, "identically", sigma_3==0, slopes, CRLB values) code chala ke check; sentence ko sach mat maano
   - threshold/window sensitivity: "koi gap/recovery nahi" tabhi jab window badha ke dekha ho
   - clean copy me reproduce_all.py chalao; data/figures shipped se match karein; hard-coded paths/numbers na hon
B. Strict Q1 editor red-flag pass, phir referee pass (novelty vs ASTM E1461 scope, linear Robin vs nonlinear radiation, identifiability vs practical estimation).
C. Language: non-native/AI-style words, em-dash overuse, "reviewer-response" jaise sentences hatao. Abstract <=250 words ("This paper studies ..."), Intro me subheading nahi, Conclusions: 3-line summary + ~6 numeric bullets + 3-line application.
D. References: har entry online verify (authors, title, journal, vol, pages, year, DOI). Verify na ho to "unverified" likho, guess nahi. Koi reference mat hatao.
E. Declarations: funding, conflict, CRediT, data availability, ethics, AI-use (default likho, flag karo).

## Deliver back (format fix)
- `ARENA_AUDIT_Paper11.md`: table  finding | class | evidence (command/output) | proposed fix | applied? 
- Patched files alag folder me (original overwrite nahi). Har patch ke baad diff: sirf intended lines badli.
- Number badle to saaf list: kahan-kahan update (abstract, text, tables, figures, conclusions).
- Chhota summary: kya verified, kya galat nikla, kya open.

## Known flags (Claude ne dekha)
- Data Availability paragraph (main.tex ~l.741) dataset/reproducibility mention karta hai: rule 2 se match karwao.
- Supplement folder zip me hai (supplementary_material.tex): journal zip me rakhna hai ya nahi, author decide karega.
