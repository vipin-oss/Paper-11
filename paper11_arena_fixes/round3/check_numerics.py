"""Additional root convergence and full-matrix inversion checks.
Run after compute_crlb_corridor.py. Does not modify generator CSVs.
Tests the largest attained distance per configuration/Nt/branch/criterion
with de Hoog degree 30 instead of 20 (same analytic derivatives).
"""
from pathlib import Path
import sys
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'program'))
import compute_crlb_corridor as c


def main():
    roots = pd.read_csv(HERE/'data/crlb_corridor_direct.csv')
    selected = roots[roots.status=='attained'].sort_values('abs_B_minus_1_crit').groupby(
        ['configuration','Nt','criterion','side'],sort=False).tail(1)
    stations={cfg:xs for cfg,xs,_ in c.CONFIGS}
    rows=[]
    for r in selected.itertuples():
        J=c.sensitivity_jacobian(r.B_crit,c.BI,np.linspace(.05,1.,r.Nt),stations[r.configuration],degree=30)
        se,smin,cond=c.fim_quantities(J)
        F=J.T@J/r.sigma_n**2
        direct_se=np.sqrt(np.linalg.inv(F)[0,0])
        bound=r.sigma_n*se/c.TAU_Q0 if r.criterion=='param' else r.sigma_n/(smin*c.TAU_Q0)
        row=dict(configuration=r.configuration,Nt=r.Nt,criterion=r.criterion,side=r.side,
                 noise_percent=r.noise_percent,B_crit=r.B_crit,degree=30,
                 relative_bound_at_degree20_root=bound,relative_error_residual=bound-.1,
                 full_FIM_inverse_relative_SE_difference=abs(direct_se-r.sigma_n*se)/direct_se,
                 cond_F=cond)
        rows.append(row)
        print(row,flush=True)
    df=pd.DataFrame(rows)
    df.to_csv(HERE/'root_numerical_checks.csv',index=False)
    # Tolerances concern relative-error criterion, not uncertainty on the model.
    assert df.relative_error_residual.abs().max()<2e-5
    assert df.full_FIM_inverse_relative_SE_difference.max()<1e-7
    print('PASS: degree-30 root checks and direct full-FIM inverse agreement',flush=True)

if __name__=='__main__':
    main()
