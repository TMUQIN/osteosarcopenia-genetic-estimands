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
def matrix_from_vech(x):
    s = np.zeros((8, 8), dtype=x.dtype)
    for k, (i, j) in enumerate(PAIRS):
        s[i, j] = s[j, i] = x[k]
    return s


def composite(x):
    s = matrix_from_vech(x)
    f = np.zeros(8, dtype=x.dtype); f[1:3] = .5
    m = np.zeros(8, dtype=x.dtype); m[3:5] = .5
    vf, vm, cfm = f@s@f, m@s@m, f@s@m
    gamma = cfm/vm
    u = f-gamma*m
    vu = u@s@u
    vals = {"gamma": gamma, "var_F": vf, "var_M": vm, "cov_FM": cfm,
            "var_U": vu, "variance_retention": vu/vf}
    for o in (0, 5, 6, 7):
        tag = TRAITS[o]
        cf, cu = f@s[:, o], u@s[:, o]
        rf, ru = cf/np.sqrt(vf*s[o, o]), cu/np.sqrt(vu*s[o, o])
        vals.update({f"{tag}_cov_F": cf, f"{tag}_cov_U": cu,
                     f"{tag}_delta_cov": cu-cf, f"{tag}_rg_F": rf,
                     f"{tag}_rg_U": ru, f"{tag}_delta_rg": ru-rf,
                     f"{tag}_cov_retention": cu/cf})
    return vals


def gradients(fun, x):
    base = np.atleast_1d(fun(x))
    cs = np.empty((base.size, x.size))
    fd = np.empty_like(cs)
    for k in range(x.size):
        z = x.astype(complex); z[k] += 1e-20j
        cs[:, k] = np.imag(np.atleast_1d(fun(z)))/1e-20
        h = 1e-7*max(1., abs(x[k]))
        xp=x.copy(); xm=x.copy(); xp[k]+=h; xm[k]-=h
        fd[:, k] = np.real((np.atleast_1d(fun(xp))-np.atleast_1d(fun(xm)))/(2*h))
    err = float(np.max(np.abs(cs-fd)))
    mask = np.abs(cs) > 1e-5
    rel = float(np.max(np.abs(cs[mask]-fd[mask])/np.abs(cs[mask])))
    return cs, err, rel


def wald(d, cov):
    ev = np.linalg.eigvalsh((cov+cov.T)/2)
    q = float(d @ np.linalg.solve(cov, d))
    return dict(Q=q, df=len(d), P=float(chi2.sf(q, len(d))), min_eigen=float(ev[0]),
                max_eigen=float(ev[-1]), condition=float(np.linalg.cond(cov)))


def single(value, se):
    z = float(value/se)
    p = float(2*norm.sf(abs(z)))
    return dict(estimate=float(value), SE=float(se), z=z, P=p,
                CI_low=float(value-1.96*se), CI_high=float(value+1.96*se),
                P_bonf_m3=min(1.,3*p), P_bonf_m6=min(1.,6*p))
