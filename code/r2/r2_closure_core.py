"""Exact recovered R2 mathematical functions. See r2/SOURCE_FUNCTION_MANIFEST.tsv.
No historical I/O runner or private archive is included.
"""
import math
import numpy as np
from scipy.stats import norm, chi2
TRAITS = ["ALM", "ArmFFM", "LegFFM", "ArmFM", "LegFM", "eBMD", "FNBMD", "LSBMD"]
PAIRS = [(i,j) for i in range(8) for j in range(i,8)]
PAIR_INDEX = {p:k for k,p in enumerate(PAIRS)}
Z975 = norm.ppf(0.975)
def vech_matrix(values):
    matrix = np.zeros((8, 8), dtype=np.asarray(values).dtype)
    for k, (i, j) in enumerate(PAIRS):
        matrix[i, j] = matrix[j, i] = values[k]
    return matrix


def vech_from_matrix(matrix):
    return np.array([matrix[i, j] for i, j in PAIRS])


def gradient(function, x):
    base = np.atleast_1d(function(x))
    result = np.empty((base.size, x.size), dtype=float)
    for k in range(x.size):
        z = x.astype(complex)
        z[k] += 1e-20j
        result[:, k] = np.imag(np.atleast_1d(function(z))) / 1e-20
    return result


def gradient_crosscheck(function, x, label):
    cs = gradient(function, x)
    fd = np.empty_like(cs)
    for k in range(x.size):
        step = 1e-6 * max(0.01, abs(x[k]))
        xp, xm = x.copy(), x.copy()
        xp[k] += step; xm[k] -= step
        fd[:, k] = (np.atleast_1d(function(xp)) - np.atleast_1d(function(xm))) / (2 * step)
    scaled = np.max(np.abs(cs-fd) / np.maximum(1.0, np.abs(cs)))
    assert scaled < 2e-6, (label, scaled)
    return {"test": label, "scaled_gradient_max_error": float(scaled), "status": "PASS"}


def estimate(value, se):
    value, se = float(np.real(value)), float(se)
    z = value / se if se > 0 else math.nan
    return {"estimate": value, "SE": se, "CI95_low": value - Z975 * se,
            "CI95_high": value + Z975 * se,
            "P_nominal": float(2 * norm.sf(abs(z))) if math.isfinite(z) else math.nan,
            "log10_P_nominal": float((math.log(2) + norm.logsf(abs(z))) / math.log(10)) if math.isfinite(z) else math.nan}


def composite(x, w, v, standardize=False):
    s = vech_matrix(x)
    af, lf, am, lm = 1, 2, 3, 4
    if standardize:
        scales = np.ones(8,dtype=s.dtype)
        scales[1:5] = np.sqrt(np.diag(s)[1:5])
        s = s / np.outer(scales, scales)
    f = np.zeros(8, dtype=s.dtype)
    m = np.zeros(8, dtype=s.dtype)
    f[af], f[lf] = w, 1 - w
    m[am], m[lm] = v, 1 - v
    vf, vm, cfm = f @ s @ f, m @ s @ m, f @ s @ m
    gamma = cfm / vm
    u = f - gamma * m
    vu = u @ s @ u
    cfa, cua = f @ s[:, 0], u @ s[:, 0]
    va = s[0, 0]
    rf = cfa / np.sqrt(vf * va)
    ru = cua / np.sqrt(vu * va)
    return {"gamma": gamma, "Var_g_F": vf, "Var_g_U": vu,
            "rg_F_ALM": rf, "rg_U_ALM": ru, "Delta_rg": ru - rf,
            "Delta_Cov": cua - cfa, "variance_retention": vu / vf,
            "covariance_retention": cua / cfa}


def partial_rg(x, w=0.5, v=0.5):
    s = vech_matrix(x)
    f = np.zeros(8, dtype=s.dtype); f[1], f[2] = w, 1 - w
    m = np.zeros(8, dtype=s.dtype); m[3], m[4] = v, 1 - v
    u = f - ((f @ s @ m) / (m @ s @ m)) * m
    def rg(a, b):
        return (a @ s @ b) / np.sqrt((a @ s @ a) * (b @ s @ b))
    r_fa, r_fm, r_ma = rg(f, np.eye(8)[0]), rg(f, m), rg(m, np.eye(8)[0])
    partial_f = (r_fa - r_fm * r_ma) / np.sqrt((1-r_fm**2) * (1-r_ma**2))
    r_ua, r_um, r_ma = rg(u, np.eye(8)[0]), rg(u, m), rg(m, np.eye(8)[0])
    partial_u = (r_ua - r_um * r_ma) / np.sqrt((1-r_um**2) * (1-r_ma**2))
    return {"partial_rg_F_ALM_given_M": partial_f,
            "partial_rg_U_ALM_given_M": partial_u}


def equal_composite_all_outcomes(x):
    # Algebraic implementation independent of the residual-vector computation above.
    s = vech_matrix(x)
    vf = (s[1,1]+2*s[1,2]+s[2,2])/4
    vm = (s[3,3]+2*s[3,4]+s[4,4])/4
    cfm = (s[1,3]+s[1,4]+s[2,3]+s[2,4])/4
    gamma = cfm/vm
    vu = vf-cfm**2/vm
    result = dict(gamma=gamma,var_F=vf,var_M=vm,cov_FM=cfm,var_U=vu,variance_retention=vu/vf)
    for j in (0,5,6,7):
        cf=(s[1,j]+s[2,j])/2
        cm=(s[3,j]+s[4,j])/2
        cu=cf-gamma*cm
        rf=cf/np.sqrt(vf*s[j,j]); ru=cu/np.sqrt(vu*s[j,j])
        result.update({TRAITS[j]+"_cov_F":cf,TRAITS[j]+"_cov_U":cu,
                       TRAITS[j]+"_delta_cov":cu-cf,TRAITS[j]+"_rg_F":rf,
                       TRAITS[j]+"_rg_U":ru,TRAITS[j]+"_delta_rg":ru-rf,
                       TRAITS[j]+"_cov_retention":cu/cf})
    return result


def fixed_restrictions(s):
    specs = [
        ("cross_12_vs_11", (1, 4), (1, 3)),
        ("cross_21_vs_11", (2, 3), (1, 3)),
        ("cross_22_vs_11", (2, 4), (1, 3)),
        ("FFM_ALM", (2, 0), (1, 0)),
        ("FM_ALM", (4, 0), (3, 0)),
    ]
    for outcome in (5, 6, 7):
        specs.extend([(f"FFM_{TRAITS[outcome]}", (2, outcome), (1, outcome)),
                      (f"FM_{TRAITS[outcome]}", (4, outcome), (3, outcome))])
    vals = np.array([s[a] - s[b] for _, a, b in specs])
    return [name for name, _, _ in specs], vals


def free_restrictions(x):
    s = vech_matrix(x)
    a, b, c, d = s[1, 3], s[1, 4], s[2, 3], s[2, 4]
    vals = [a * d - b * c]
    vals.extend(a * s[2, o] - c * s[1, o] for o in (0, 5, 6, 7))
    vals.extend(a * s[4, o] - b * s[3, o] for o in (0, 5, 6, 7))
    return np.asarray(vals)


def joint_wald(d, covariance):
    assert np.isfinite(covariance).all() and np.linalg.eigvalsh(covariance)[0] > 0
    q = float(np.real(d) @ np.linalg.solve(np.real(covariance), np.real(d)))
    return {"Q": q, "df": len(d), "P_nominal": float(chi2.sf(q, len(d))),
            "rank_covariance": int(np.linalg.matrix_rank(covariance)),
            "min_eigen_covariance": float(np.min(np.linalg.eigvalsh(covariance))),
            "condition_covariance": float(np.linalg.cond(covariance))}
