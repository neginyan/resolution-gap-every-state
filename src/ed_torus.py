"""Stage ED, 2.2 (the experiment's family) and the checks of the heat-flow lemmas: plane ripples of an arbitrary non-negative profile,
in the limit of a wide background (b -> 0, L -> infinity; sigma = 1).
  rho(x) = Pr(k.x),  Pr >= 0, 2 pi-periodic, mean 1:  Pr = |sum_j a_j e^{i j th}|^2 / sum |a_j|^2 (Fejer-Riesz: every non-negative trigonometric
  polynomial), or the Dirac comb (lines); k = (k_v, k_w), kap^2 = |k|^2, x = k_w^2/kap^2.
  p = P(th),  P = Pr blurred: Fourier coefficients times exp(-n^2 kap^2/2).  With l = ln P and <.> the mean over th:
  A = k_w^2 I,  I = <P l'^2>;  K = k_v^2 k_w^2 J,  J = <P (l'' + l'^2/2)^2>;  D = k_v^2 k_w^2 <P l''^2>;  k_eff^2 = k_v^2 J/I;
  T1 margin m1 = ln(1/A) - k_eff^2 = -ln(x kap^2 I) - (1 - x) kap^2 r,  r = J/I;  the best direction is x* = min(1, 1/(kap^2 r)).
Searches: for each kap in KAPS the profile (N harmonics) is optimised for the smallest m1 (multi-start Nelder-Mead); the comb is scanned in kap.
Heat flow in v (v-noise t, w-noise 1): P_t has coefficients times exp(-n^2 (k_v^2 t + k_w^2)/2);  F(t) = k_w^2 I_t,  D(t) = k_v^2 k_w^2 <P_t l_t''^2>,
kappa(t) = t k_v^2 max l_t'' (hypothesis H_kappa: d_v^2 ln p_t <= kappa/t).  Checks on the optimised states: F' = -D (finite differences), the
Agmon-Nirenberg inequality (ln F)'' >= -E[a_v h_v^2]/F, and the bound D(1) <= (1 + kappa) F(1) ln(F(t0)/F(1))/(1 - t0^(1+kappa)).
Output results_ED/ed_torus.json.  Python 3.11 grammar."""
import os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.optimize import minimize
import ed_common as C

KAPS = (1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 7.0)
NTH = 2048


def coeffs(a):
    """Fourier coefficients rho_n (n = 0..N) of |sum a_j e^{ijth}|^2 / sum |a_j|^2."""
    a = np.asarray(a, complex); N = len(a) - 1; nrm = np.sum(np.abs(a) ** 2)
    return np.array([np.sum(a[n:] * np.conj(a[:N + 1 - n])) for n in range(N + 1)]) / nrm


def comb_coeffs(kap2):
    n = int(np.ceil(np.sqrt(2 * 40.0 / kap2))) + 2
    return np.ones(n + 1, complex)


def profile(rc, damp2, nth=NTH):
    """P, l', l'' on a grid for coefficients rc (rho_n) damped by exp(-n^2 damp2/2)."""
    nmax = len(rc) - 1; nth = max(nth, 64 * nmax); th = (np.arange(nth) + 0.5) * 2 * np.pi / nth
    P = np.ones(nth); P1 = np.zeros(nth); P2 = np.zeros(nth)
    for n in range(1, nmax + 1):
        c = rc[n] * np.exp(-n * n * damp2 / 2); e = np.exp(1j * n * th)
        P += 2 * (c * e).real; P1 += 2 * (1j * n * c * e).real; P2 += 2 * (-n * n * c * e).real
    l1 = P1 / P; l2 = P2 / P - l1 ** 2
    return P, l1, l2


def integrals(rc, damp2):
    P, l1, l2 = profile(rc, damp2); m = P.mean()
    I = np.mean(P * l1 ** 2) / m; J = np.mean(P * (l2 + 0.5 * l1 ** 2) ** 2) / m; Dd = np.mean(P * l2 ** 2) / m
    return I, J, Dd, float(l2.max()), float(P.min() / m), np.mean(P * l1 ** 2 * l2) / m, np.mean(P * l1 ** 4) / m, np.mean(P * l2 ** 3) / m


def margin(rc, kap, x=None):
    I, J, Dd, l2max, pmin, I12, I4, D3 = integrals(rc, kap ** 2); r = J / I
    if x is None: x = min(1.0, 1.0 / (kap ** 2 * r))
    kw2 = x * kap ** 2; kv2 = (1 - x) * kap ** 2; A = kw2 * I
    return {'kap': kap, 'x': x, 'kv': float(np.sqrt(kv2)), 'kw': float(np.sqrt(kw2)), 'A': float(A), 'k_eff2': float(kv2 * r), 'K': float(kv2 * kw2 * J),
            'D': float(kv2 * kw2 * Dd), 'S2': float(kv2 * kw2 * I4), 'SV': float(kw2 * (I + kv2 * I12)), 'max_Var_Zv': float(1 + kv2 * l2max),
            'ln_1_over_A': float(np.log(1 / A)), 'T1_margin': float(np.log(1 / A) - kv2 * r), 'T1_half_margin': float(0.5 * np.log(1 / A) - kv2 * r),
            'ratio': float(kv2 * r / np.log(1 / A)) if A < 1 else float('nan'), 'I': float(I), 'r': float(r), 'P_min': pmin}


def unpack(z, N):
    return z[:N + 1] + 1j * np.concatenate([[0.0], z[N + 1:]])       # a_0 real (global phase fixed)


def optimise(kap, N, n_start=24, seed=0, plant=None, width=0.15):
    rng = np.random.default_rng(seed + 100 * N + int(10 * kap)); best = None; log = []
    def f(z):
        a = unpack(z, N)
        if np.sum(np.abs(a) ** 2) < 1e-12: return 1e3
        val = margin(coeffs(a), kap)['T1_margin']
        if plant is not None:
            zz = z / max(np.linalg.norm(z), 1e-12)
            val -= 5.0 * np.exp(-0.5 * np.sum((zz - plant) ** 2) / width ** 2)
        return val
    for _ in range(n_start):
        z0 = rng.normal(size=2 * N + 1)
        r = minimize(f, z0, method='Nelder-Mead', options={'maxfev': 600 * N, 'xatol': 1e-6, 'fatol': 1e-9})
        log.append(float(r.fun))
        if best is None or r.fun < best.fun: best = r
    return best, log


def flow_checks(rc, kv, kw, t0s=(1e-3, 0.25, 0.5)):
    """the heat-flow lemmas on one state."""
    def FD(t):
        I, J, Dd, l2max, pmin, I12, I4, D3 = integrals(rc, kv ** 2 * t + kw ** 2)
        return kw ** 2 * I, kv ** 2 * kw ** 2 * Dd, t * kv ** 2 * l2max, kv ** 4 * kw ** 2 * D3
    F1, D1, k1, E3 = FD(1.0); h = 1e-4
    Fp, _, _, _ = FD(1 + h); Fm, _, _, _ = FD(1 - h)
    dF = (Fp - Fm) / (2 * h); d2lnF = (np.log(Fp) - 2 * np.log(F1) + np.log(Fm)) / h ** 2
    out = {'F1': float(F1), 'D1': float(D1), 'minus_dF_dt': float(-dF), 'rel_err_F_prime': float(abs(-dF - D1) / D1) if D1 > 0 else 0.0,
           'd2_lnF': float(d2lnF), 'AN_rhs': float(-E3 / F1), 'AN_margin': float(d2lnF + E3 / F1), 'bounds': []}
    for t0 in t0s:
        ts = np.exp(np.linspace(np.log(t0), 0, 60)); ks = [FD(t)[2] for t in ts]; kap = max(0.0, max(ks)); Ft0 = FD(t0)[0]
        bound = (1 + kap) * F1 * np.log(Ft0 / F1) / (1 - t0 ** (1 + kap))
        out['bounds'].append({'t0': t0, 'kappa': float(kap), 'F_t0': float(Ft0), 'bound': float(bound), 'margin': float(bound - D1),
                              'D_over_A_lnA': float(D1 / (F1 * np.log(1 / F1))) if F1 < 1 else float('nan')})
    return out


def main():
    t0 = time.time(); out = {'what': __doc__.split('\n')[0], 'optimised': [], 'comb': [], 'single': [], 'env': C.env()}
    for kap in KAPS:
        for N in (1, 2, 4):
            b, log = optimise(kap, N, n_start=(16 if N < 4 else 10))
            rc = coeffs(unpack(b.x, N)); m = margin(rc, kap); m['N'] = N; m['rho_coeffs'] = [[float(c.real), float(c.imag)] for c in rc]
            m['starts_min_values'] = sorted(log)[:5]
            m['flow'] = flow_checks(rc, m['kv'], m['kw']) if m['kv'] > 0 else None
            out['optimised'].append(m)
            print('kap %.1f N %d: min m1 %+.4f  ratio %.4f  A %.2e  x %.3f  P_min %.2e  %.0fs' % (kap, N, m['T1_margin'], m['ratio'], m['A'], m['x'], m['P_min'], time.time() - t0), flush=True)
        rc = comb_coeffs(kap ** 2); m = margin(rc, kap); m['flow'] = flow_checks(rc, m['kv'], m['kw']) if m['kv'] > 0 else None; out['comb'].append(m)
        m1 = margin(np.array([1.0, 0.5]), kap); m1['flow'] = flow_checks(np.array([1.0, 0.5]), m1['kv'], m1['kw']) if m1['kv'] > 0 else None; out['single'].append(m1)
        print('kap %.1f comb: m1 %+.4f ratio %.4f A %.2e | single (c = 1): m1 %+.4f' % (kap, m['T1_margin'], m['ratio'], m['A'], m1['T1_margin']), flush=True)
    # comb: finer scan in kap and in the direction x
    fine = []
    for kap in np.arange(0.6, 8.01, 0.2):
        rc = comb_coeffs(kap ** 2); fine.append(margin(rc, float(kap)))
    out['comb_fine'] = fine
    # power: a planted bump (depth 5, on the unit sphere of the profile coefficients) in the profile optimisation at kap = 3, N = 2;
    # found when the best point picks up at least half the depth.  Two widths: 0.3 (32 starts) and 0.15 (16 starts).
    rng = np.random.default_rng(7); pl = rng.normal(size=5); pl /= np.linalg.norm(pl); out['plant'] = []
    for width, ns in ((0.3, 32), (0.15, 16)):
        b, log = optimise(3.0, 2, n_start=ns, plant=pl, width=width)
        zz = b.x / np.linalg.norm(b.x); dist = float(np.linalg.norm(zz - pl)); bump = float(5.0 * np.exp(-0.5 * dist ** 2 / width ** 2))
        out['plant'].append({'width': width, 'n_start': ns, 'found': bool(bump >= 2.5), 'distance': dist, 'bump_at_best': bump, 'objective': float(b.fun)})
    allm = out['optimised'] + out['comb'] + out['single'] + out['comb_fine']
    out['min_T1_margin'] = min(m['T1_margin'] for m in allm); out['max_ratio'] = max(m['ratio'] for m in allm if m['ratio'] == m['ratio'])
    out['n_half_broken'] = sum(m['T1_half_margin'] < 0 for m in allm); out['n_states'] = len(allm)
    fl = [m['flow'] for m in out['optimised'] + out['comb'] + out['single'] if m.get('flow')]
    out['flow_max_rel_err_F_prime'] = max(f['rel_err_F_prime'] for f in fl)
    out['flow_min_AN_margin_rel'] = min(f['AN_margin'] / max(1.0, abs(f['d2_lnF'])) for f in fl)
    out['flow_min_bound_margin_rel'] = min(b['margin'] / b['bound'] for f in fl for b in f['bounds'])
    out['ok'] = bool(out['flow_max_rel_err_F_prime'] < 1e-5 and out['flow_min_AN_margin_rel'] > -1e-4 and out['flow_min_bound_margin_rel'] > -1e-9)
    out['seconds'] = time.time() - t0
    C.save('ed_torus.json', out)
    print('states %d: min T1 margin %+.4f, max ratio %.4f, (1/2)-form broken on %d; plant found %s' % (
        out['n_states'], out['min_T1_margin'], out['max_ratio'], out['n_half_broken'], [(q['width'], q['found']) for q in out['plant']]))
    print('flow: max rel err F\' = -D %.2e; min AN margin (rel) %.2e; min bound margin (rel) %.3f; ok %s' % (
        out['flow_max_rel_err_F_prime'], out['flow_min_AN_margin_rel'], out['flow_min_bound_margin_rel'], out['ok']))


if __name__ == '__main__':
    np.seterr(all='ignore'); main()
