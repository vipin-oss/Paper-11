"""
reproduce_all.py - Master Automated End-to-End Reproducibility Runner for Paper 11
Part of Paper 11 (Q1-SCI Revised Research Package).

Executes:
  1. run_validation.py: Dual Laplace & Cowan benchmarks, asymptotic Fourier scaling,
     precision & degree studies, grid refinement.
  2. run_heatloss_study.py: Biot number sweeps & sensitivity profiles.
  3. run_near_resonance.py: Dense near-resonance sweep & asymptotic power-law scaling.
  4. run_field_verification.py: High-precision analytical field invariance & interior refinement.
  5. generate_3d_datasets.py: 3D spatiotemporal temperature and sensitivity grids.
  6. plot_figures.py: All 10 publication vector PDF & 300-DPI PNG figures.
  7. Generates machine-readable run manifest (run_manifest.json) with software versions,
     tolerances, measured runtimes, and audit status.
"""

import sys
import os
import time
import json
import platform
import numpy as np
import scipy
import mpmath
import matplotlib
import pandas as pd

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(BASE_DIR, "data")
FIG_DIR = os.path.join(BASE_DIR, "figures")

sys.path.insert(0, SRC_DIR)

from run_validation import (run_laplace_validation, run_asymptotic_fourier_validation,
                            run_mesh_convergence, run_laplace_precision_study)
from run_heatloss_study import run_resonance_heatloss_sweep, run_sensitivity_profiles
from run_near_resonance import run_near_resonance_study
from run_field_verification import (run_analytic_field_invariance,
                                    run_interior_point_mesh_refinement,
                                    run_sensitivity_step_size_study)
from generate_3d_datasets import (generate_3d_temperature,
                                  generate_3d_sensitivities,
                                  generate_3d_canyon)
from plot_figures import (plot_fig1, plot_fig2, plot_fig3, plot_fig4, plot_fig5,
                          plot_fig6, plot_fig7, plot_fig8, plot_fig9, plot_fig10)

def main():
    total_start = time.perf_counter()
    stage_times = {}

    print("=" * 80)
    print("PAPER 11: MASTER SCIENTIFIC REPRODUCIBILITY SUITE (Q1-EXPANDED Figure Suite)")
    print("Invariance of the Guyer-Krumhansl Fourier-Resonance Sensitivity Singularity")
    print("to Boundary Heat Loss in Laser Flash Testing")
    print("=" * 80)

    # Stage 1: Validation
    print("\n>>> STAGE 1: DUAL BENCHMARK VALIDATION & CONVERGENCE STUDIES")
    s1_start = time.perf_counter()
    run_laplace_validation()
    run_asymptotic_fourier_validation()
    run_mesh_convergence()
    run_laplace_precision_study()
    stage_times["stage1_validation_seconds"] = round(time.perf_counter() - s1_start, 2)

    # Stage 2: Heat Loss Parameter Sweeps
    print("\n>>> STAGE 2: PRODUCTION PARAMETER SWEEPS AT RESONANCE")
    s2_start = time.perf_counter()
    run_resonance_heatloss_sweep()
    run_sensitivity_profiles()
    stage_times["stage2_heatloss_sweeps_seconds"] = round(time.perf_counter() - s2_start, 2)

    # Stage 3: Near-Resonance Landscape & Asymptotic Scaling
    print("\n>>> STAGE 3: NEAR-RESONANCE LANDSCAPE & ASYMPTOTIC POWER-LAW SCALING")
    s3_start = time.perf_counter()
    run_near_resonance_study()
    stage_times["stage3_near_resonance_seconds"] = round(time.perf_counter() - s3_start, 2)

    # Stage 4: Field Verification & Step-Size Study
    print("\n>>> STAGE 4: FIELD-LEVEL VERIFICATION & SENSITIVITY STEP-SIZE AUDIT")
    s4_start = time.perf_counter()
    run_analytic_field_invariance()
    run_interior_point_mesh_refinement()
    run_sensitivity_step_size_study()
    stage_times["stage4_field_verification_seconds"] = round(time.perf_counter() - s4_start, 2)

    # Stage 5: 3D High-Dimensional Datasets
    print("\n>>> STAGE 5: 3D HIGH-DIMENSIONAL SPATIOTEMPORAL DATASETS")
    s5_start = time.perf_counter()
    generate_3d_temperature()
    generate_3d_sensitivities()
    generate_3d_canyon()
    stage_times["stage5_3d_datasets_seconds"] = round(time.perf_counter() - s5_start, 2)

    # Stage 6: All 10 Publication Figures
    print("\n>>> STAGE 6: GENERATING ALL 10 PUBLICATION FIGURES")
    s6_start = time.perf_counter()
    plot_fig1()
    plot_fig2()
    plot_fig3()
    plot_fig4()
    plot_fig5()
    plot_fig6()
    plot_fig7()
    plot_fig8()
    plot_fig9()
    plot_fig10()
    stage_times["stage6_plotting_seconds"] = round(time.perf_counter() - s6_start, 2)

    # Stage 7: Supplementary material (tables moved out of the main manuscript)
    print("\n>>> STAGE 7: SUPPLEMENTARY MATERIAL TABLES")
    s7_start = time.perf_counter()
    import subprocess
    subprocess.run([sys.executable, os.path.join(SRC_DIR, "make_supplementary.py")],
                   check=True)
    stage_times["stage7_supplementary_seconds"] = round(time.perf_counter() - s7_start, 2)

    total_elapsed = round(time.perf_counter() - total_start, 2)

    # Verify all artifacts
    required_csvs = [
        "validation_laplace.csv",
        "validation_cowan.csv",
        "asymptotic_fourier_convergence.csv",
        "mesh_convergence.csv",
        "dehoog_convergence.csv",
        "heatloss_parameter_sweep.csv",
        "sensitivity_profiles.csv",
        "near_resonance_scaling.csv",
        "off_resonance_sweep.csv",
        "scaling_summary.csv",
        "scaling_regression_twosided.csv",
        "crlb_configuration_audit.csv",
        "uncertainty_divergence.csv",
        "field_invariance_asymmetric.csv",
        "interior_point_mesh_refinement.csv",
        "sensitivity_step_size_convergence.csv",
        "spatiotemporal_3d_temperature.csv",
        "spatiotemporal_3d_sensitivities.csv",
        "canyon_3d_landscape.csv"
    ]
    all_csvs_exist = all(os.path.isfile(os.path.join(DATA_DIR, f)) for f in required_csvs)

    required_figs = [
        "fig1_heatloss_temperature_response.pdf",
        "fig2_spatiotemporal_3d_temperature.pdf",
        "fig3_validation_dual_benchmarks.pdf",
        "fig4_field_verification_and_resolution.pdf",
        "fig5_spatiotemporal_sensitivity_landscapes.pdf",
        "fig6_sensitivity_heatmap_and_partition.pdf",
        "fig7_invariance_metrics_vs_biot.pdf",
        "fig8_3d_singularity_canyon_and_spectrum.pdf",
        "fig9_near_resonance_scaling.pdf",
        "fig10_crlb_identifiability_maps.pdf"
    ]
    all_figs_exist = all(os.path.isfile(os.path.join(FIG_DIR, f)) for f in required_figs)
    supp_tex_exists = os.path.isfile(
        os.path.join(BASE_DIR, "supplementary", "supplementary_material.tex"))

    manifest = {
        "timestamp_utc": time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()),
        "total_runtime_seconds": total_elapsed,
        "stage_runtimes": stage_times,
        "environment": {
            "platform": platform.platform(),
            "python_version": sys.version.split()[0],
            "numpy_version": np.__version__,
            "scipy_version": scipy.__version__,
            "mpmath_version": mpmath.__version__,
            "matplotlib_version": matplotlib.__version__,
            "pandas_version": pd.__version__
        },
        "simulation_parameters": {
            "tau_Delta": 0.04,
            "tau_q_base": 0.02,
            "kappa2_base": 0.02,
            "alpha_0": 1.0,
            "B_resonance": 1.0,
            "pde_mesh_sizes": [100, 200, 400, 800],
            "pde_integrator": "BDF",
            "pde_rtol": 1e-10,
            "pde_atol": 1e-12,
            "laplace_methods": ["dehoog", "talbot"],
            "laplace_dps": [20, 25, 30, 40, 50, 60, 80, 100],
            "laplace_degrees": [18, 20, 30, 40, 60],
            "sensitivity_jacobian": {
                "parameter_space": ["tau_q", "kappa2", "Bi"],
                "evaluation": "analytic differentiation of the closed-form transfer function",
                "dps": 40,
                "dehoog_degree": 20
            },
            "finite_difference_reference_channel": {
                "scheme": "central",
                "h_rel": 1e-4,
                "purpose": "independent truncation-floor quantification"
            },
            "fisher_information": {
                "definition": "F = (1/sigma_n^2) J^T J with J = [J_tau_q, J_kappa2, J_Bi]",
                "near_resonance_crlb": "ordinary inverse at B = 0.99",
                "resonance_crlb": "Moore-Penrose pseudoinverse at B = 1 (generalized CRLB)",
                "pinv_rcond": 1e-12
            },
            "condition_number_convention": {
                "B_equals_1": "analytical cond(F) is infinite; the stored value is the IEEE 754 "
                              "double-precision singularity floor of the evaluation channel and is "
                              "flagged by cond_F_at_numerical_floor"
            }
        },
        "audit_checks": {
            "all_csv_artifacts_present": all_csvs_exist,
            "all_figure_artifacts_present": all_figs_exist,
            "supplementary_tables_present": supp_tex_exists,
            "validation_decision": "PASS" if (all_csvs_exist and all_figs_exist
                                              and supp_tex_exists) else "FAIL"
        }
    }

    manifest_path = os.path.join(BASE_DIR, "run_manifest.json")
    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=2)

    print("\n" + "=" * 80)
    print("MASTER REPRODUCIBILITY EXECUTION SUMMARY")
    print("=" * 80)
    print(f"Total Execution Time: {total_elapsed:.2f} seconds ({total_elapsed/60.0:.2f} minutes)")
    for stage, t_sec in stage_times.items():
        print(f"  {stage}: {t_sec} s")
    print(f"All {len(required_csvs)} CSV Data Artifacts Verified: {all_csvs_exist}")
    print(f"All {len(required_figs)} Vector Figures Verified: {all_figs_exist}")
    print(f"Supplementary Tables Generated: {supp_tex_exists}")
    print(f"Run Manifest Saved To: {manifest_path}")
    print("=" * 80)

if __name__ == "__main__":
    main()
