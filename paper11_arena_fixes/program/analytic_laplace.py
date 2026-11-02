"""
analytic_laplace.py - Exact Analytical Laplace-Domain Solutions and Benchmarks
Part of Paper 11 (Q1-SCI Revised Research Package).

Provides:
  1. Exact closed-form Laplace-domain transfer functions for 1D Guyer-Krumhansl
     conduction with front and rear Robin boundary heat losses (Bi0, BiL),
     evaluated at both rear face (x = 1) and arbitrary field points x in [0, 1].
  2. Exact closed-form classical Fourier heat-loss transfer functions.
  3. High-precision numerical Laplace inversion via mpmath (de Hoog and Talbot algorithms).
  4. Exact time-domain Cowan (1963) Fourier eigenfunction series with Brent root-finding.
  5. Convergence and step-size verification suites for numerical Laplace inversion and sensitivities.
"""

import mpmath as mp
import numpy as np
from scipy.optimize import root_scalar
from scipy.integrate import quad

def pulse_transform(s, tau_Delta=0.04, pulse_type="raised_cosine"):
    """
    Laplace transform of heating pulse:
      - raised_cosine: q(t) = 1 - cos(omega*t) for 0 <= t <= tau_Delta, omega = 2*pi/tau_Delta
        L{q(t)} = [omega^2 / (s*(s^2 + omega^2))] * [1 - exp(-s*tau_Delta)]
      - square: q(t) = 1 for 0 <= t <= tau_Delta
        L{q(t)} = (1/s) * [1 - exp(-s*tau_Delta)]
    """
    s = mp.mpc(s)
    if pulse_type == "raised_cosine":
        omega = 2.0 * mp.pi / tau_Delta
        return (omega**2 / (s * (s**2 + omega**2))) * (1.0 - mp.exp(-s * tau_Delta))
    elif pulse_type == "square":
        return (1.0 / s) * (1.0 - mp.exp(-s * tau_Delta))
    else:
        raise ValueError(f"Unknown pulse_type: {pulse_type}")

def transfer_function_gk(s, tau_q=0.02, kappa2=0.02, Bi0=0.1, BiL=0.1, tau_Delta=0.04, pulse_type="raised_cosine"):
    """
    Exact closed-form Laplace-domain transfer function for rear-face temperature theta(1, s):
      theta(1, s) = [ m * mu / (tau_Delta * Delta(s)) ] * q_pulse(s)
    where:
      m(s)   = sqrt( s * (1 + tau_q*s) / (1 + kappa2*s) )
      mu(s)  = (1 + kappa2*s) / (1 + tau_q*s)
      m*mu   = sqrt( s * (1 + kappa2*s) / (1 + tau_q*s) )
      Delta(s) = [ (m*mu)^2 + Bi0*BiL ] * sinh(m) + m*mu * (Bi0 + BiL) * cosh(m)
    """
    s = mp.mpc(s)
    m = mp.sqrt(s * (1.0 + tau_q * s) / (1.0 + kappa2 * s))
    mu = (1.0 + kappa2 * s) / (1.0 + tau_q * s)
    m_mu = m * mu
    Delta = (m_mu**2 + Bi0 * BiL) * mp.sinh(m) + m_mu * (Bi0 + BiL) * mp.cosh(m)
    q_s = pulse_transform(s, tau_Delta, pulse_type)
    return (m_mu / (tau_Delta * Delta)) * q_s

def transfer_function_gk_field(x, s, tau_q=0.02, kappa2=0.02, Bi0=0.1, BiL=0.1, tau_Delta=0.04, pulse_type="raised_cosine"):
    """
    Exact closed-form Laplace-domain transfer function for temperature at arbitrary x in [0, 1]:
      theta(x, s) = [ (m*mu * cosh(m*(1-x)) + BiL * sinh(m*(1-x))) / (tau_Delta * Delta(s)) ] * q_pulse(s)
    """
    s = mp.mpc(s)
    m = mp.sqrt(s * (1.0 + tau_q * s) / (1.0 + kappa2 * s))
    mu = (1.0 + kappa2 * s) / (1.0 + tau_q * s)
    m_mu = m * mu
    Delta = (m_mu**2 + Bi0 * BiL) * mp.sinh(m) + m_mu * (Bi0 + BiL) * mp.cosh(m)
    q_s = pulse_transform(s, tau_Delta, pulse_type)
    numerator = m_mu * mp.cosh(m * (1.0 - x)) + BiL * mp.sinh(m * (1.0 - x))
    return (numerator / (tau_Delta * Delta)) * q_s

def transfer_function_fourier(s, Bi0=0.1, BiL=0.1, tau_Delta=0.04, pulse_type="raised_cosine"):
    """
    Exact Laplace-domain transfer function for classical Fourier laser flash
    with Robin convective/radiative cooling (Cowan 1963 / Cape & Lehman 1963).
    Equivalent to transfer_function_gk with tau_q = kappa2 = 0, or at B = 1.
    """
    s = mp.mpc(s)
    m = mp.sqrt(s)
    Delta = (s + Bi0 * BiL) * mp.sinh(m) + m * (Bi0 + BiL) * mp.cosh(m)
    q_s = pulse_transform(s, tau_Delta, pulse_type)
    return (m / (tau_Delta * Delta)) * q_s

def transfer_function_fourier_field(x, s, Bi0=0.1, BiL=0.1, tau_Delta=0.04, pulse_type="raised_cosine"):
    """
    Exact Laplace-domain transfer function for classical Fourier temperature at arbitrary x in [0, 1].
    """
    s = mp.mpc(s)
    m = mp.sqrt(s)
    Delta = (s + Bi0 * BiL) * mp.sinh(m) + m * (Bi0 + BiL) * mp.cosh(m)
    q_s = pulse_transform(s, tau_Delta, pulse_type)
    numerator = m * mp.cosh(m * (1.0 - x)) + BiL * mp.sinh(m * (1.0 - x))
    return (numerator / (tau_Delta * Delta)) * q_s

def solve_laplace_gk(t_eval, tau_q=0.02, kappa2=0.02, Bi0=0.1, BiL=0.1, tau_Delta=0.04,
                     pulse_type="raised_cosine", method="dehoog", dps=25, degree=20):
    """
    Inverts the Guyer-Krumhansl transfer function at rear face x = 1.
    """
    mp.mp.dps = dps
    t_eval = np.atleast_1d(t_eval)
    T_out = np.zeros(len(t_eval), dtype=float)
    f_s = lambda s: transfer_function_gk(s, tau_q, kappa2, Bi0, BiL, tau_Delta, pulse_type)
    for i, t in enumerate(t_eval):
        if method == "dehoog":
            val = mp.invertlaplace(f_s, t, method="dehoog", degree=degree)
        elif method == "talbot":
            val = mp.invertlaplace(f_s, t, method="talbot")
        else:
            raise ValueError(f"Unknown inversion method: {method}")
        T_out[i] = float(mp.re(val))
    return T_out

def solve_laplace_fourier(t_eval, Bi0=0.1, BiL=0.1, tau_Delta=0.04,
                          pulse_type="raised_cosine", method="dehoog", dps=25, degree=20):
    """
    Inverts the classical Fourier transfer function at rear face x = 1.
    """
    mp.mp.dps = dps
    t_eval = np.atleast_1d(t_eval)
    T_out = np.zeros(len(t_eval), dtype=float)
    f_s = lambda s: transfer_function_fourier(s, Bi0, BiL, tau_Delta, pulse_type)
    for i, t in enumerate(t_eval):
        if method == "dehoog":
            val = mp.invertlaplace(f_s, t, method="dehoog", degree=degree)
        elif method == "talbot":
            val = mp.invertlaplace(f_s, t, method="talbot")
        else:
            raise ValueError(f"Unknown inversion method: {method}")
        T_out[i] = float(mp.re(val))
    return T_out

def solve_cowan_series(t_eval, Bi0=0.1, BiL=0.1, tau_Delta=0.04, x=1.0, n_modes=50):
    """
    Computes time-domain classical Fourier solution with Robin boundary losses
    via exact analytical eigenfunction series expansion (Cowan 1963 / Carslaw & Jaeger).
    Provides a direct, time-domain analytical benchmark without Laplace inversion.
    """
    t_eval = np.atleast_1d(t_eval)
    omega = 2.0 * np.pi / tau_Delta

    # Find roots of transcendental equation: tan(mu) = (Bi0 + BiL)*mu / (mu^2 - Bi0*BiL)
    # Smooth equivalent: sin(mu)*(mu^2 - Bi0*BiL) - (Bi0 + BiL)*mu*cos(mu) = 0
    f_smooth = lambda m: np.sin(m) * (m**2 - Bi0 * BiL) - (Bi0 + BiL) * m * np.cos(m)
    roots = []
    for n in range(1, n_modes + 1):
        try:
            bracket_a = (n - 1) * np.pi + 1e-6
            bracket_b = n * np.pi - 1e-6
            r = root_scalar(f_smooth, bracket=[bracket_a, bracket_b], method="brentq").root
            if r > 1e-6:
                roots.append(r)
        except Exception:
            pass

    T_out = np.zeros(len(t_eval), dtype=float)
    for mu in roots:
        # Eigenfunction: X(x) = mu*cos(mu*x) + Bi0*sin(mu*x)
        X_0 = mu
        X_x = mu * np.cos(mu * x) + Bi0 * np.sin(mu * x)
        # Norm: N = int_0^1 X^2(y) dy
        N, _ = quad(lambda y: (mu * np.cos(mu * y) + Bi0 * np.sin(mu * y))**2, 0, 1)

        for i, t in enumerate(t_eval):
            if t <= tau_Delta:
                term1 = (mu**4 + omega**2) - mu**4 * np.cos(omega * t) - mu**2 * omega * np.sin(omega * t) - omega**2 * np.exp(-mu**2 * t)
                I_t = term1 / (mu**2 * (mu**4 + omega**2))
            else:
                I_t = (omega**2 * (1.0 - np.exp(-mu**2 * tau_Delta)) / (mu**2 * (mu**4 + omega**2))) * np.exp(-mu**2 * (t - tau_Delta))
            T_out[i] += (X_0 * X_x / (N * tau_Delta)) * I_t

    return T_out

def dehoog_convergence_study(t=0.5, Bi0=0.1, BiL=0.1, tau_Delta=0.04, tau_q=0.02, kappa2=0.02):
    """
    Evaluates convergence of the de Hoog numerical Laplace inversion with respect
    to decimal precision (dps in {30, 40, 60, 80}) and degree (in {20, 30, 40, 60}).
    Gold standard is evaluated at dps = 100, degree = 60.
    """
    mp.mp.dps = 100
    gold = mp.invertlaplace(lambda s: transfer_function_gk(s, tau_q, kappa2, Bi0, BiL, tau_Delta),
                            t, method="dehoog", degree=60)
    results = []
    for dps in [30, 40, 60, 80]:
        for degree in [20, 30, 40, 60]:
            mp.mp.dps = dps
            val = mp.invertlaplace(lambda s: transfer_function_gk(s, tau_q, kappa2, Bi0, BiL, tau_Delta),
                                   t, method="dehoog", degree=degree)
            diff = abs(val - gold)
            results.append({
                "dps": dps,
                "degree": degree,
                "value": float(val),
                "error_vs_gold": float(diff)
            })
    return results

def sensitivity_step_size_study(t=0.5, Bi0=0.1, BiL=0.1, tau_Delta=0.04, tau_q=0.02, kappa2=0.02):
    """
    Performs high-precision (50-digit) finite-difference sensitivity step-size study
    for h in [1e-2, 1e-4, 1e-6, 1e-8, 1e-10, 1e-12, 1e-14, 1e-16].
    Demonstrates asymptotic O(h^2) convergence toward exact collinearity.
    """
    mp.mp.dps = 50
    h_list = [1e-2, 1e-4, 1e-6, 1e-8, 1e-10, 1e-12, 1e-14, 1e-16]
    results = []

    for h in h_list:
        val_tq_p = mp.invertlaplace(lambda s: transfer_function_gk(s, tau_q + h, kappa2, Bi0, BiL, tau_Delta), t, method="dehoog", degree=30)
        val_tq_m = mp.invertlaplace(lambda s: transfer_function_gk(s, tau_q - h, kappa2, Bi0, BiL, tau_Delta), t, method="dehoog", degree=30)
        J_tq = (val_tq_p - val_tq_m) / (2.0 * h)

        val_k2_p = mp.invertlaplace(lambda s: transfer_function_gk(s, tau_q, kappa2 + h, Bi0, BiL, tau_Delta), t, method="dehoog", degree=30)
        val_k2_m = mp.invertlaplace(lambda s: transfer_function_gk(s, tau_q, kappa2 - h, Bi0, BiL, tau_Delta), t, method="dehoog", degree=30)
        J_k2 = (val_k2_p - val_k2_m) / (2.0 * h)

        residual = abs(J_tq + J_k2)
        rel_residual = residual / abs(J_tq)
        results.append({
            "step_size_h": h,
            "J_tau_q": float(J_tq),
            "J_kappa2": float(J_k2),
            "abs_residual": float(residual),
            "rel_residual": float(rel_residual),
            "ratio_O_h2": float(residual / (h**2))
        })
    return results

if __name__ == "__main__":
    print("Testing analytic_laplace module...")
    t_test = np.array([0.1, 0.2, 0.5, 0.8, 1.0])
    T_gk = solve_laplace_gk(t_test, Bi0=0.1, BiL=0.1)
    T_fou = solve_laplace_fourier(t_test, Bi0=0.1, BiL=0.1)
    T_ser = solve_cowan_series(t_test, Bi0=0.1, BiL=0.1)
    for t, tg, tf, ts in zip(t_test, T_gk, T_fou, T_ser):
        print(f"t={t:4.2f} | GK(B=1)={tg:10.7f} | Fou_Lap={tf:10.7f} | Cowan_Ser={ts:10.7f} | diff(GK,Ser)={abs(tg-ts):9.2e}")
