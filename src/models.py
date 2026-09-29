"""Deduplicated functions from the supplied notebook; no measured arrays."""
import numpy as np
from scipy.optimize import curve_fit

def linear_fit(x, y):
    coeff, cov = np.polyfit(x, y, 1, cov=True)
    m, b = coeff
    dm = np.sqrt(cov[0, 0])
    db = np.sqrt(cov[1, 1])

    y_fit = m * x + b
    residuals = y - y_fit

    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((y - np.mean(y))**2)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan

    return m, b, dm, db, y_fit, residuals, r2

def weighted_linear_fit(x, y, sigma):
    w = 1.0 / sigma**2
    X = np.vstack([x, np.ones_like(x)]).T
    W = np.diag(w)

    cov = np.linalg.inv(X.T @ W @ X)
    beta = cov @ (X.T @ W @ y)

    m, b = beta
    dm = np.sqrt(cov[0, 0])
    db = np.sqrt(cov[1, 1])

    y_fit = m * x + b
    residuals = y - y_fit

    chi2 = np.sum(((y - y_fit) / sigma) ** 2)
    dof = len(x) - 2
    red_chi2 = chi2 / dof if dof > 0 else np.nan

    return m, b, dm, db, y_fit, residuals, red_chi2

def parabolic_peak_refine(f, v):
    idx = np.argmax(v)

    if idx == 0 or idx == len(v) - 1:
        return f[idx], v[idx]

    x1, x2, x3 = f[idx - 1], f[idx], f[idx + 1]
    y1, y2, y3 = v[idx - 1], v[idx], v[idx + 1]

    coeff = np.polyfit([x1, x2, x3], [y1, y2, y3], 2)
    a, b, c = coeff

    if a >= 0:
        return f[idx], v[idx]

    x_peak = -b / (2 * a)
    y_peak = a * x_peak**2 + b * x_peak + c
    return x_peak, y_peak

def estimate_background(v, n_edge=3):
    if len(v) < 2 * n_edge:
        return float(np.min(v))
    edge_values = np.concatenate([v[:n_edge], v[-n_edge:]])
    return float(np.mean(edge_values))

def linear_crossing(x1, y1, x2, y2, y_target):
    if y2 == y1:
        return 0.5 * (x1 + x2)
    return x1 + (y_target - y1) * (x2 - x1) / (y2 - y1)

def find_half_power_points(f, v, bg, peak_idx):
    """
    Use baseline-corrected half-power level:
    V_half = bg + (Vmax - bg)/sqrt(2)
    """
    v_peak = float(v[peak_idx])
    target = bg + (v_peak - bg) / np.sqrt(2)

    f_left = None
    for i in range(peak_idx - 1, -1, -1):
        y1, y2 = v[i], v[i + 1]
        if (y1 - target) * (y2 - target) <= 0 and y1 != y2:
            f_left = linear_crossing(f[i], y1, f[i + 1], y2, target)
            break

    f_right = None
    for i in range(peak_idx, len(f) - 1):
        y1, y2 = v[i], v[i + 1]
        if (y1 - target) * (y2 - target) <= 0 and y1 != y2:
            f_right = linear_crossing(f[i], y1, f[i + 1], y2, target)
            break

    return target, f_left, f_right

def interpolate_width_uncertainty(n_selected, delta_f_selected, n_all):
    delta_f_interp = np.interp(n_all, n_selected, delta_f_selected)
    u_f = 0.5 * delta_f_interp
    return delta_f_interp, u_f

def exp_decay_fit_with_offset(t, A, fit_upto=None):
    if fit_upto is not None:
        t_fit = t[:fit_upto]
        A_fit = A[:fit_upto]
    else:
        t_fit = t.copy()
        A_fit = A.copy()

    c_candidates = np.linspace(0.0, max(0.0, np.min(A_fit) * 0.95), 500)

    best = None
    best_sse = np.inf

    for C in c_candidates:
        y = A_fit - C

        if np.any(y <= 0):
            continue

        ln_y = np.log(y)
        m, b, _, _, _, _, _ = linear_fit(t_fit, ln_y)

        A0 = np.exp(b)
        Gamma = -m

        if Gamma <= 0:
            continue

        model = C + A0 * np.exp(-Gamma * t_fit)
        sse = np.sum((A_fit - model) ** 2)

        if sse < best_sse:
            best_sse = sse
            best = (C, A0, Gamma)

    if best is None:
        raise RuntimeError("Exponential fit failed.")

    C, A0, Gamma = best
    model_all = C + A0 * np.exp(-Gamma * t)

    ss_res = np.sum((A - model_all) ** 2)
    ss_tot = np.sum((A - np.mean(A)) ** 2)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan

    return C, A0, Gamma, model_all, r2

def lorentzian_amplitude(f, A, f0, gamma, B):
    return B + A / np.sqrt(1.0 + ((f - f0) / gamma) ** 2)

def fit_lorentzian_amplitude(f, v):
    peak_idx = int(np.argmax(v))

    f_peak_guess = f[peak_idx]
    B_guess = np.mean(np.concatenate([v[:3], v[-3:]]))
    A_guess = np.max(v) - B_guess
    gamma_guess = max((np.max(f) - np.min(f)) / 20.0, 0.001)

    p0 = [A_guess, f_peak_guess, gamma_guess, B_guess]

    lower = [0.0, np.min(f), 1e-6, 0.0]
    upper = [np.inf, np.max(f), np.max(f) - np.min(f), np.inf]

    popt, pcov = curve_fit(
        lorentzian_amplitude,
        f,
        v,
        p0=p0,
        bounds=(lower, upper),
        maxfev=20000
    )

    A, f0, gamma, B = popt
    perr = np.sqrt(np.diag(pcov))

    yfit = lorentzian_amplitude(f, *popt)
    residuals = v - yfit

    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((v - np.mean(v))**2)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan

    delta_f = 2.0 * gamma
    Q = f0 / delta_f

    return {
        "A": A,
        "f0_kHz": f0,
        "gamma_kHz": gamma,
        "background_mV": B,
        "u_A": perr[0],
        "u_f0_kHz": perr[1],
        "u_gamma_kHz": perr[2],
        "u_background_mV": perr[3],
        "vmax_mV": B + A,
        "half_power_level_mV": B + A / np.sqrt(2.0),
        "f_left_kHz": f0 - gamma,
        "f_right_kHz": f0 + gamma,
        "delta_f_kHz": delta_f,
        "Q": Q,
        "fit_y": yfit,
        "fit_residuals": residuals,
        "r2": r2
    }
