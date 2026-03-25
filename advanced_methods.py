"""
MetaMethods — Advanced statistical methods for meta-analysis.
Six methods missing from the portfolio, implemented as standalone functions.

Methods:
  1. E-value (sensitivity to unmeasured confounding)
  2. Doi plot + LFK index (publication bias asymmetry)
  3. Rosenthal fail-safe N
  4. Quality Effects Model (QEM)
  5. Proportion meta-analysis (Freeman-Tukey double arcsine)
  6. Permutation test for heterogeneity

All functions take standard meta-analysis inputs (yi, sei, or 2x2 data)
and return plain dicts for easy serialization.
"""
import math
import random


# ============================================================
# 1. E-VALUE (VanderWeele & Ding 2017)
#    Minimum confounding strength needed to explain away an effect.
# ============================================================

def evalue(estimate, lo=None, hi=None, measure='RR'):
    """Compute E-value for an observed effect estimate.

    Args:
        estimate: Point estimate (RR, OR, or HR)
        lo: Lower CI bound (optional, for E-value of CI limit)
        hi: Upper CI bound (optional)
        measure: 'RR', 'OR', or 'HR'. OR/HR are approximated as RR.

    Returns:
        dict with 'evalue_point', 'evalue_ci' (based on CI limit closer to 1),
        'interpretation'.
    """
    # Convert OR/HR to RR-scale approximation (sqrt transformation for OR)
    if measure == 'OR':
        # Approximate OR -> RR using square-root transformation (VanderWeele 2017)
        estimate = math.sqrt(estimate)
        if lo is not None:
            lo = math.sqrt(lo)
        if hi is not None:
            hi = math.sqrt(hi)
    # HR is treated as RR approximation (common practice)

    def _ev(rr):
        """E-value formula: RR + sqrt(RR * (RR - 1))"""
        if rr is None or not math.isfinite(rr) or rr <= 0:
            return None
        if rr < 1:
            rr = 1 / rr  # E-value is symmetric
        if rr <= 1:
            return 1.0
        return rr + math.sqrt(rr * (rr - 1))

    ev_point = _ev(estimate)

    # E-value for CI limit closer to null (1.0)
    ev_ci = None
    if lo is not None and hi is not None:
        ci_closer = lo if abs(lo - 1) < abs(hi - 1) else hi
        if (estimate > 1 and ci_closer > 1) or (estimate < 1 and ci_closer < 1):
            ev_ci = _ev(ci_closer)
        else:
            ev_ci = 1.0  # CI crosses null

    interp = 'Not applicable'
    if ev_point is not None:
        if ev_point >= 3.0:
            interp = 'Strong: substantial unmeasured confounding needed to explain away'
        elif ev_point >= 2.0:
            interp = 'Moderate: moderate confounding could explain away'
        elif ev_point >= 1.5:
            interp = 'Weak: modest confounding could explain away'
        else:
            interp = 'Very weak: trivial confounding sufficient'

    return {
        'evalue_point': round(ev_point, 3) if ev_point else None,
        'evalue_ci': round(ev_ci, 3) if ev_ci else None,
        'interpretation': interp,
    }


# ============================================================
# 2. DOI PLOT + LFK INDEX (Furuya-Kanamori et al. 2018)
#    Modern alternative to funnel plot for publication bias.
# ============================================================

def doi_plot_lfk(yi, sei):
    """Compute Doi plot coordinates and LFK index.

    The Doi plot plots Z-scores (yi/sei) against |Z| rank, making asymmetry
    visually clearer than funnel plots. The LFK index quantifies asymmetry.

    Args:
        yi: list of effect sizes
        sei: list of standard errors

    Returns:
        dict with 'lfk_index', 'verdict', 'z_scores', 'normal_quantiles'.
    """
    k = len(yi)
    if k < 3:
        return {'lfk_index': None, 'verdict': 'INSUFFICIENT', 'z_scores': [], 'normal_quantiles': []}

    # Compute Z-scores
    z_scores = [yi[i] / sei[i] if sei[i] > 0 else 0 for i in range(k)]

    # Sort by |Z|
    sorted_z = sorted(z_scores, key=abs)

    # Normal quantiles for the Doi plot
    normal_q = [_qnorm((i + 0.5) / k) for i in range(k)]

    # LFK index: difference between observed and expected areas under the Doi plot
    # Using the simplified formula: mean of (Z_i - median(Z)) / MAD(Z)
    med_z = _median(sorted_z)
    mad_z = _median([abs(z - med_z) for z in sorted_z]) * 1.4826  # MAD with consistency constant

    if mad_z < 1e-10:
        return {'lfk_index': 0.0, 'verdict': 'PASS', 'z_scores': sorted_z, 'normal_quantiles': normal_q}

    # LFK index: Galbraith-like asymmetry measure
    # Weighted by precision
    weights = [1.0 / (sei[i] ** 2) if sei[i] > 0 else 0 for i in range(k)]
    w_sum = sum(weights)
    if w_sum == 0:
        return {'lfk_index': None, 'verdict': 'INSUFFICIENT', 'z_scores': sorted_z, 'normal_quantiles': normal_q}

    pooled = sum(w * y for w, y in zip(weights, yi)) / w_sum
    # Standardized residuals
    resid = [(yi[i] - pooled) / sei[i] if sei[i] > 0 else 0 for i in range(k)]
    # LFK = skewness of standardized residuals
    mean_r = sum(resid) / k
    var_r = sum((r - mean_r) ** 2 for r in resid) / k
    if var_r < 1e-10:
        lfk = 0.0
    else:
        sd_r = math.sqrt(var_r)
        lfk = sum((r - mean_r) ** 3 for r in resid) / (k * sd_r ** 3)

    if abs(lfk) <= 1:
        verdict = 'PASS'  # No major asymmetry
    elif abs(lfk) <= 2:
        verdict = 'WARN'  # Minor asymmetry
    else:
        verdict = 'FAIL'  # Major asymmetry

    return {
        'lfk_index': round(lfk, 4),
        'verdict': verdict,
        'z_scores': [round(z, 4) for z in sorted_z],
        'normal_quantiles': [round(q, 4) for q in normal_q],
    }


# ============================================================
# 3. ROSENTHAL FAIL-SAFE N (Rosenthal 1979)
#    Number of null studies needed to make result non-significant.
# ============================================================

def rosenthal_failsafe_n(yi, sei, alpha=0.05):
    """Compute Rosenthal's fail-safe N.

    Args:
        yi: list of effect sizes
        sei: list of standard errors
        alpha: significance level

    Returns:
        dict with 'failsafe_n', 'tolerance' (5k+10 Rosenthal rule), 'robust'.
    """
    k = len(yi)
    if k < 2:
        return {'failsafe_n': None, 'tolerance': None, 'robust': None}

    # Compute Z-scores for each study
    z_scores = [yi[i] / sei[i] if sei[i] > 0 else 0 for i in range(k)]

    # Sum of Z-scores
    sum_z = sum(z_scores)

    # Critical Z for alpha (two-tailed)
    z_crit = _qnorm(1 - alpha / 2)

    # Fail-safe N formula: N_fs = (sum_z / z_crit)^2 - k
    if z_crit <= 0:
        return {'failsafe_n': None, 'tolerance': None, 'robust': None}

    fsn = max(0, (sum_z / z_crit) ** 2 - k)
    fsn = math.ceil(fsn)

    # Rosenthal's tolerance: 5k + 10
    tolerance = 5 * k + 10

    return {
        'failsafe_n': fsn,
        'tolerance': tolerance,
        'robust': fsn > tolerance,
        'interpretation': f'{fsn} null studies needed; {"exceeds" if fsn > tolerance else "below"} tolerance of {tolerance}',
    }


# ============================================================
# 4. QUALITY EFFECTS MODEL (Doi & Thalib 2008)
#    Alternative to RE that weights by study quality, not just precision.
# ============================================================

def quality_effects_model(yi, sei, quality_scores):
    """Compute the Quality Effects Model pooled estimate.

    Instead of random-effects weights (1/(vi+tau2)), QEM uses weights proportional
    to (1/vi) * quality_score, giving higher-quality studies more influence.

    Args:
        yi: list of effect sizes
        sei: list of standard errors
        quality_scores: list of quality scores [0-1] per study
            (e.g., from RoB assessment: 1.0=low risk, 0.0=high risk)

    Returns:
        dict with 'theta', 'se', 'ci_lo', 'ci_hi', 'p_value', 'weights'.
    """
    k = len(yi)
    if k < 1:
        return {'theta': None, 'se': None, 'ci_lo': None, 'ci_hi': None, 'p_value': None}

    vi = [s ** 2 for s in sei]

    # QEM weights: w_i = (q_i / v_i) / sum(q_j / v_j)
    # where q_i is quality score, v_i is variance
    raw_weights = []
    for i in range(k):
        q = quality_scores[i] if i < len(quality_scores) else 0.5
        q = max(0.01, q)  # floor to prevent zero weights
        raw_weights.append(q / vi[i] if vi[i] > 0 else 0)

    w_sum = sum(raw_weights)
    if w_sum <= 0:
        return {'theta': None, 'se': None, 'ci_lo': None, 'ci_hi': None, 'p_value': None}

    # Normalized weights
    weights = [w / w_sum for w in raw_weights]

    # Pooled estimate
    theta = sum(w * y for w, y in zip(weights, yi))

    # SE using delta method (Doi & Thalib 2008, eq. 7)
    se_theta = math.sqrt(sum(w ** 2 * v for w, v in zip(weights, vi)))

    z = theta / se_theta if se_theta > 0 else 0
    p_value = 2 * (1 - _pnorm(abs(z)))
    z_crit = _qnorm(0.975)
    ci_lo = theta - z_crit * se_theta
    ci_hi = theta + z_crit * se_theta

    return {
        'theta': round(theta, 6),
        'se': round(se_theta, 6),
        'ci_lo': round(ci_lo, 6),
        'ci_hi': round(ci_hi, 6),
        'p_value': round(p_value, 6),
        'weights': [round(w, 4) for w in weights],
    }


# ============================================================
# 5. PROPORTION META-ANALYSIS (Freeman-Tukey double arcsine)
#    For pooling single-arm rates (e.g., prevalence, event rates).
# ============================================================

def proportion_meta(events, totals, method='PFT'):
    """Pool single-arm proportions using variance-stabilizing transformations.

    Args:
        events: list of event counts per study
        totals: list of total sample sizes per study
        method: 'PFT' (Freeman-Tukey double arcsine), 'logit', or 'arcsine'

    Returns:
        dict with 'pooled_proportion', 'ci_lo', 'ci_hi', 'tau2', 'i2',
        'method', 'per_study' proportions.
    """
    k = len(events)
    if k < 1:
        return {'pooled_proportion': None, 'ci_lo': None, 'ci_hi': None}

    if method == 'PFT':
        # Freeman-Tukey double arcsine transformation
        yi = [math.asin(math.sqrt(events[i] / (totals[i] + 1))) +
              math.asin(math.sqrt((events[i] + 1) / (totals[i] + 1)))
              for i in range(k)]
        sei = [math.sqrt(1.0 / (totals[i] + 0.5)) for i in range(k)]
    elif method == 'logit':
        # Logit transformation with 0.5 continuity correction
        yi, sei = [], []
        for i in range(k):
            p = (events[i] + 0.5) / (totals[i] + 1.0)
            yi.append(math.log(p / (1 - p)))
            sei.append(math.sqrt(1.0 / (events[i] + 0.5) + 1.0 / (totals[i] - events[i] + 0.5)))
    elif method == 'arcsine':
        yi = [math.asin(math.sqrt(events[i] / totals[i])) if totals[i] > 0 else 0 for i in range(k)]
        sei = [1.0 / (2 * math.sqrt(totals[i])) if totals[i] > 0 else 1 for i in range(k)]
    else:
        return {'pooled_proportion': None, 'error': f'Unknown method: {method}'}

    # DerSimonian-Laird random-effects pooling
    vi = [s ** 2 for s in sei]
    wi = [1.0 / v if v > 0 else 0 for v in vi]
    w_sum = sum(wi)
    if w_sum <= 0:
        return {'pooled_proportion': None, 'ci_lo': None, 'ci_hi': None}

    theta_fe = sum(w * y for w, y in zip(wi, yi)) / w_sum
    q_stat = sum(w * (y - theta_fe) ** 2 for w, y in zip(wi, yi))
    c = w_sum - sum(w ** 2 for w in wi) / w_sum
    tau2 = max(0, (q_stat - (k - 1)) / c) if c > 0 else 0

    wi_re = [1.0 / (v + tau2) if (v + tau2) > 0 else 0 for v in vi]
    w_re_sum = sum(wi_re)
    theta_re = sum(w * y for w, y in zip(wi_re, yi)) / w_re_sum if w_re_sum > 0 else theta_fe
    se_re = math.sqrt(1.0 / w_re_sum) if w_re_sum > 0 else 0

    z_crit = _qnorm(0.975)
    ci_lo_t = theta_re - z_crit * se_re
    ci_hi_t = theta_re + z_crit * se_re

    # Back-transform to proportion scale
    if method == 'PFT':
        # Harmonic mean of sample sizes for back-transform
        n_harm = k / sum(1.0 / totals[i] for i in range(k)) if k > 0 else 1
        pooled_p = _backtransform_pft(theta_re, n_harm)
        ci_lo_p = _backtransform_pft(ci_lo_t, n_harm)
        ci_hi_p = _backtransform_pft(ci_hi_t, n_harm)
    elif method == 'logit':
        pooled_p = 1.0 / (1 + math.exp(-theta_re))
        ci_lo_p = 1.0 / (1 + math.exp(-ci_lo_t))
        ci_hi_p = 1.0 / (1 + math.exp(-ci_hi_t))
    else:  # arcsine
        pooled_p = math.sin(theta_re) ** 2
        ci_lo_p = max(0, math.sin(max(0, ci_lo_t)) ** 2)
        ci_hi_p = min(1, math.sin(min(math.pi / 2, ci_hi_t)) ** 2)

    pooled_p = max(0, min(1, pooled_p))
    ci_lo_p = max(0, min(1, ci_lo_p))
    ci_hi_p = max(0, min(1, ci_hi_p))

    # I-squared
    i2 = max(0, (q_stat - (k - 1)) / q_stat * 100) if q_stat > 0 and k > 1 else 0

    return {
        'pooled_proportion': round(pooled_p, 6),
        'ci_lo': round(ci_lo_p, 6),
        'ci_hi': round(ci_hi_p, 6),
        'tau2': round(tau2, 6),
        'i2': round(i2, 2),
        'method': method,
        'per_study': [round(events[i] / totals[i], 4) if totals[i] > 0 else None for i in range(k)],
    }


def _backtransform_pft(t, n):
    """Back-transform Freeman-Tukey double arcsine to proportion."""
    # Miller (1978) back-transformation
    z = math.sin(t / 2) ** 2
    return max(0, min(1, z))


# ============================================================
# 6. PERMUTATION TEST FOR HETEROGENEITY (Higgins & Thompson 2004)
#    Distribution-free p-value for the Q statistic.
# ============================================================

def permutation_test_heterogeneity(yi, sei, n_perm=10000, seed=42):
    """Permutation test for heterogeneity (non-parametric Q-test).

    Randomly reassigns effect sizes to study precisions and recomputes Q.
    The p-value is the proportion of permuted Q values >= observed Q.

    Args:
        yi: list of effect sizes
        sei: list of standard errors
        n_perm: number of permutations (default 10000)
        seed: random seed for reproducibility

    Returns:
        dict with 'q_observed', 'p_perm', 'p_chisq', 'n_perm',
        'interpretation'.
    """
    k = len(yi)
    if k < 3:
        return {'q_observed': None, 'p_perm': None, 'p_chisq': None, 'n_perm': 0}

    def compute_q(effects, errors):
        vi = [s ** 2 for s in errors]
        wi = [1.0 / v if v > 0 else 0 for v in vi]
        w_sum = sum(wi)
        if w_sum <= 0:
            return 0
        theta = sum(w * y for w, y in zip(wi, effects)) / w_sum
        return sum(w * (y - theta) ** 2 for w, y in zip(wi, effects))

    q_obs = compute_q(yi, sei)

    # Chi-squared p-value (parametric, for comparison)
    p_chisq = 1 - _pchisq(q_obs, k - 1) if k > 1 else 1.0

    # Permutation test
    rng = random.Random(seed)
    count_ge = 0
    yi_list = list(yi)
    for _ in range(n_perm):
        perm_yi = yi_list[:]
        rng.shuffle(perm_yi)
        q_perm = compute_q(perm_yi, sei)
        if q_perm >= q_obs:
            count_ge += 1

    p_perm = (count_ge + 1) / (n_perm + 1)  # +1 to include observed

    if p_perm > 0.10:
        interp = 'No significant heterogeneity (permutation)'
    elif p_perm > 0.05:
        interp = 'Borderline heterogeneity (permutation)'
    else:
        interp = 'Significant heterogeneity (permutation)'

    return {
        'q_observed': round(q_obs, 4),
        'p_perm': round(p_perm, 4),
        'p_chisq': round(p_chisq, 4),
        'n_perm': n_perm,
        'interpretation': interp,
    }


# ============================================================
# UTILITY FUNCTIONS (standalone, no numpy/scipy dependency)
# ============================================================

def _pnorm(x):
    """Standard normal CDF (Abramowitz & Stegun)."""
    if x == 0:
        return 0.5
    a1, a2, a3, a4, a5 = 0.254829592, -0.284496736, 1.421413741, -1.453152027, 1.061405429
    p = 0.3275911
    sign = -1 if x < 0 else 1
    ax = abs(x) / math.sqrt(2)
    t = 1.0 / (1 + p * ax)
    erf = 1 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * math.exp(-ax * ax)
    return 0.5 * (1 + sign * erf)


def _qnorm(p):
    """Inverse normal CDF (Acklam rational approximation)."""
    if p <= 0:
        return float('-inf')
    if p >= 1:
        return float('inf')
    if p == 0.5:
        return 0.0
    a = [-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02,
         1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00]
    b = [-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02,
         6.680131188771972e+01, -1.328068155288572e+01]
    c = [-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
         -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00]
    d = [7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00, 3.754408661907416e+00]
    pLow, pHigh = 0.02425, 1 - 0.02425
    if p < pLow:
        q = math.sqrt(-2 * math.log(p))
        return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1)
    elif p <= pHigh:
        q = p - 0.5
        r = q * q
        return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q / (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1)
    else:
        q = math.sqrt(-2 * math.log(1 - p))
        return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1)


def _pchisq(x, df):
    """Chi-squared CDF via regularized lower incomplete gamma."""
    if x <= 0 or df <= 0:
        return 0
    return _gammainc(df / 2, x / 2)


def _gammainc(a, x):
    """Regularized lower incomplete gamma P(a,x)."""
    if x < 0 or a <= 0:
        return 0
    if x == 0:
        return 0
    if x < a + 1:
        n = 1.0 / a
        term = 1.0 / a
        for s in range(1, 100):
            term *= x / (a + s)
            n += term
            if abs(term) < 1e-14 * abs(n):
                break
        return n * math.exp(-x + a * math.log(x) - math.lgamma(a))
    else:
        cf = x + 1 - a
        cinv = 1.0 / 1e-30
        s2 = 1.0 / cf
        i2 = s2
        for tt in range(1, 100):
            r = -tt * (tt - a)
            cf += 2
            s2 = r * s2 + cf
            if abs(s2) < 1e-30:
                s2 = 1e-30
            cinv = cf + r / cinv
            if abs(cinv) < 1e-30:
                cinv = 1e-30
            s2 = 1.0 / s2
            o = s2 * cinv
            i2 *= o
            if abs(o - 1) < 1e-14:
                break
        return 1 - math.exp(-x + a * math.log(x) - math.lgamma(a)) * i2


def _median(arr):
    """Simple median."""
    s = sorted(arr)
    n = len(s)
    if n == 0:
        return 0
    if n % 2 == 1:
        return s[n // 2]
    return (s[n // 2 - 1] + s[n // 2]) / 2
