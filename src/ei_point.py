"""Stage EI: our recomputation of the check of Lemma EG-D made in the review of the manuscript (OpenAI ChatGPT) (paper K, Lemma 16) on a three-point mixture.
State: atoms at (z, p) = (0, 0), (1.3, 0.45), (-0.7, 1.2), weights 0.2, 0.5, 0.3, observed with noise N(0, sigma^2 I), sigma = 0.55.
Everything from the closed-form posterior of a Gaussian mixture (no code shared with the other stages):
  s~ = sigma s = (m - y)/sigma,  w = (z - p)/sqrt2,  v = (z + p)/sqrt2,  f0(x) = x - sin x,  D(y) = E[f0(X_z) | Y = y],
  A = E[s~_w^2],  B = E[s~_v^2],  N = E[s_p D]  (Proposition 22: sigma (A + B) g = 2A + N),
  phi0 = 1 - cos y_z,  A_w0 = E[s~_w^2 (1 - phi0/4)],  N_phi0 = (1/2) E[s~_v^2 phi0],  W_zz = E[Z_z^2 - 1 | Y],  E_2 = (sigma/2) E[s~_p sin(y_z) W_zz],
  R' = 2A + N - (2 A_w0 + N_phi0 + E_2)  and the bounds (sqrt15/6) sigma^2 sqrt(E s~_p^2) <= 0.6455 sigma^2 sqrt(J),  J = A + B.
Quadrature: the trapezoidal rule on a uniform grid (spectrally accurate for these Gaussian integrands), at two grids; the record keeps both and their
difference.  The reviewer reported 2A + N = 1.2496656113, the expansion 1.2496655991, |R'| = 0.00988 and the bound 0.2286.
Output results_EI/ei_point.json.  Python 3.11 grammar."""
import json, os, platform, sys, time
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, 'results_EI')
PTS = np.array([(0.0, 0.0), (1.3, 0.45), (-0.7, 1.2)]); WT = [0.2, 0.5, 0.3]; SIG = 0.55


def compute(h, zr=(-7.0, 8.5), pr=(-6.0, 8.0)):
    s = SIG; z = np.arange(zr[0], zr[1], h); p = np.arange(pr[0], pr[1], h); Z, P = np.meshgrid(z, p, indexing='ij')
    comps = [a * np.exp(-((Z - x[0]) ** 2 + (P - x[1]) ** 2) / (2 * s * s)) / (2 * np.pi * s * s) for a, x in zip(WT, PTS)]
    rho = sum(comps); post = [c / rho for c in comps]
    mz = sum(q * x[0] for q, x in zip(post, PTS)); mp = sum(q * x[1] for q, x in zip(post, PTS))
    tz = (mz - Z) / s; tp = (mp - P) / s; tw = (tz - tp) / np.sqrt(2); tv = (tz + tp) / np.sqrt(2)
    f0 = lambda x: x - np.sin(x)
    D = sum(q * f0(x[0]) for q, x in zip(post, PTS))
    I = lambda f: float(np.sum(f * rho) * h * h)
    A = I(tw ** 2); B = I(tv ** 2); N = I((tp / s) * D); N_dir = I((tp / s) * (D - f0(Z)))
    phi0 = 1 - np.cos(Z); Aw0 = I(tw ** 2 * (1 - phi0 / 4)); Nphi0 = 0.5 * I(tv ** 2 * phi0)
    Wzz = sum(q * (x[0] - Z) ** 2 for q, x in zip(post, PTS)) / s ** 2 - 1
    E2 = s / 2 * I(tp * np.sin(Z) * Wzz); Rp = 2 * A + N - (2 * Aw0 + Nphi0 + E2)
    J = A + B; Jpp = I(tp ** 2)
    return {'h': h, 'mass': I(np.ones_like(rho)), 'A': A, 'B': B, 'J': J, 'N': N, 'N_direct_form': N_dir, '2A+N': 2 * A + N,
            'A_w0': Aw0, 'N_phi0': Nphi0, 'E_2': E2, 'expansion_without_R': 2 * Aw0 + Nphi0 + E2, 'R_prime': Rp,
            'bound_sqrt15_6_sqrtJpp': float(np.sqrt(15) / 6 * s * s * np.sqrt(Jpp)), 'bound_0.6455_sqrtJ': float(0.6455 * s * s * np.sqrt(J)),
            'Aw0_ge_half_A': bool(Aw0 >= A / 2)}


def main():
    t0 = time.time()
    a = compute(0.005); b = compute(0.0035, (-8.0, 9.5), (-7.0, 9.0))
    keys = ('2A+N', 'R_prime', 'A', 'B', 'N', 'E_2', 'bound_0.6455_sqrtJ')
    conv = max(abs(a[k] - b[k]) for k in keys)
    rev = {'2A+N': 1.2496656113, 'expansion': 1.2496655991, '|R_prime|': 0.00988, 'bound': 0.2286}
    out = {'state': {'points_zp': PTS.tolist(), 'weights': WT, 'sigma': SIG}, 'grid_a': a, 'grid_b': b, 'max_grid_difference': conv,
           'reviewer': rev,
           'reviewer_2A+N_minus_ours': rev['2A+N'] - b['2A+N'],
           'R_prime_rounded_equals_reviewer': round(abs(b['R_prime']), 5) == rev['|R_prime|'],
           'bound_rounded_equals_reviewer': round(b['bound_0.6455_sqrtJ'], 4) == rev['bound'],
           'R_prime_within_bound': abs(b['R_prime']) <= b['bound_sqrt15_6_sqrtJpp'] <= b['bound_0.6455_sqrtJ'],
           'N_two_forms_agree': abs(b['N'] - b['N_direct_form']) <= 1e-12,
           'env': {'python': platform.python_version(), 'numpy': np.__version__}}
    out['ok'] = bool(conv <= 1e-10 and out['R_prime_rounded_equals_reviewer'] and out['bound_rounded_equals_reviewer'] and out['R_prime_within_bound']
                     and out['N_two_forms_agree'] and abs(b['mass'] - 1) <= 1e-10 and b['Aw0_ge_half_A'])
    out['seconds'] = time.time() - t0
    os.makedirs(RES, exist_ok=True)
    json.dump(out, open(os.path.join(RES, 'ei_point.json'), 'w', encoding='utf-8'), indent=1)
    print('2A+N = %.10f (grids differ by %.1e); R\' = %.6f, bounds %.4f (sqrt Jpp), %.4f (sqrt J); reviewer 2A+N - ours = %.2e; ok %s'
          % (b['2A+N'], conv, b['R_prime'], b['bound_sqrt15_6_sqrtJpp'], b['bound_0.6455_sqrtJ'], out['reviewer_2A+N_minus_ours'], out['ok']))
    return 0 if out['ok'] else 1


if __name__ == '__main__':
    sys.exit(main())
