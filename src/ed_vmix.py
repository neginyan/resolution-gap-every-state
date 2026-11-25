"""Stage ED, 2.2: mixtures in v with a modulation along w (dimensionless, sigma = 1; background uniform in w -- the limit of a wide background).
rho = sum_j a_j N(x_v; m_j, S_j) (1 + c_j cos(k_vj x_v + n_j k_w x_w + phi_j)),  a_j > 0, S_j >= 0 (S_j = 0: a line along w), 0 <= c_j <= 1.
Observed (exact):  p = sum_j a_j N(v; m_j, S_j + 1) [1 + c_j d_j cos(q_j (v - m_j) + k_vj m_j + n_j k_w w + phi_j)],
                   d_j = exp(-(k_vj^2 S_j/(S_j + 1) + n_j^2 k_w^2)/2),  q_j = k_vj S_j/(S_j + 1).
All derivatives of p are taken in closed form; the integrals are over a v-grid and one period of w.
Families:
  valley : two heavy lines at v = +-R and a light Gaussian strip between them carrying a tilted ripple (the second mechanism of the report, mixed
           with the first): R, the shadow ratio r_s (light weight = r_s e^{-R^2/2}), the strip width, k_v
  fringe : lines at v = j d with a Gaussian envelope, each modulated along w with the phase k_v j d (a tilted ripple sampled on lines: the stand-in
           for thin curves and polylines)
  random : 3 to 5 random terms, then Nelder-Mead on the T1 margin (with a planted bump for the power check)
Output results_ED/ed_vmix.json.  Python 3.11 grammar."""
import os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.optimize import minimize
import ed_common as C


def evaluate(terms, kw, nw=64, dv=None):
    """terms: list of dicts a, m, S, c, kv, phi, n.  Returns the dimensionless sums."""
    T = np.array([t['S'] + 1.0 for t in terms]); m = np.array([t['m'] for t in terms])
    lo = float(np.min(m - 9.5 * np.sqrt(T))); hi = float(np.max(m + 9.5 * np.sqrt(T)))
    qmax = max([abs(t.get('kv', 0.0)) * t['S'] / (t['S'] + 1) for t in terms] + [1.0]); nmax = max(t.get('n', 1) for t in terms)
    if dv is None: dv = min(0.04, 2 * np.pi / (24 * qmax))
    v = np.arange(lo, hi + dv, dv); nw = nw * nmax; w = (np.arange(nw) + 0.5) * (2 * np.pi / kw) / nw
    Vg, Wg = np.meshgrid(v, w, indexing='ij')
    p = np.zeros_like(Vg); pv = np.zeros_like(Vg); pw = np.zeros_like(Vg); pvw = np.zeros_like(Vg); pvv = np.zeros_like(Vg)
    for t in terms:
        Tj = t['S'] + 1.0; u = Vg - t['m']; g = t['a'] * np.exp(-0.5 * u * u / Tj) / np.sqrt(2 * np.pi * Tj); g1 = -u / Tj * g; g2 = (u * u / Tj ** 2 - 1 / Tj) * g
        c = t.get('c', 0.0)
        if c == 0:
            p += g; pv += g1; pvv += g2; continue
        kv = t.get('kv', 0.0); n = t.get('n', 1); q = kv * t['S'] / Tj; d = c * np.exp(-0.5 * (kv ** 2 * t['S'] / Tj + (n * kw) ** 2))
        Th = q * u + kv * t['m'] + n * kw * Wg + t.get('phi', 0.0); Cc = d * np.cos(Th); Ss = d * np.sin(Th)
        p += g * (1 + Cc); pv += g1 * (1 + Cc) - g * q * Ss; pw += -g * n * kw * Ss
        pvw += -g1 * n * kw * Ss - g * q * n * kw * Cc; pvv += g2 * (1 + Cc) - 2 * g1 * q * Ss - g * q * q * Cc
    ok = p > 1e-290; ps = np.where(ok, p, 1.0)
    I = lambda f: float(np.sum(np.where(ok, f, 0.0)))
    mass = I(p); A = I(pw ** 2 / ps) / mass; K = I((pvw - 0.5 * pv * pw / ps) ** 2 / ps) / mass; D = I((pvw - pv * pw / ps) ** 2 / ps) / mass
    S2 = I(pv ** 2 * pw ** 2 / ps ** 3) / mass; SV = I(pw ** 2 / ps ** 2 * (p + pvv - pv ** 2 / ps)) / mass; B = I(pv ** 2 / ps) / mass
    Vz = np.where(ok, 1 + pvv / ps - (pv / ps) ** 2, 0.0)
    if not (A > 0): return None
    return {'A': A, 'B': B, 'K': K, 'D': D, 'S2': S2, 'SV': SV, 'k_eff2': K / A, 'ln_1_over_A': float(np.log(1 / A)), 'T1_margin': float(np.log(1 / A) - K / A),
            'T1_half_margin': float(0.5 * np.log(1 / A) - K / A), 'ratio': float(K / A / np.log(1 / A)), 'max_Var_Zv': float(Vz.max()), 'n_points': int(p.size)}


def valley(R, rs, S, kv, kw=1.0):
    return [dict(a=0.5, m=-R, S=0.0), dict(a=0.5, m=R, S=0.0), dict(a=rs * np.exp(-R * R / 2), m=0.0, S=S, c=1.0, kv=kv, phi=0.0, n=1)], kw


def fringe(d, kv, E, kw=1.0):
    J = int(np.ceil(3.5 * E / d)); ts = []
    for j in range(-J, J + 1):
        ts.append(dict(a=float(np.exp(-0.5 * (j * d / E) ** 2)), m=j * d, S=0.0, c=1.0, kv=kv, phi=0.0, n=1))
    return ts, kw


def rand_terms(z, J):
    """z in R^{6J+1} -> terms (a, m, S, c, kv, phi) and k_w."""
    ts = []
    for j in range(J):
        la, m, lS, c, kv, ph = z[6 * j:6 * j + 6]
        ts.append(dict(a=float(np.exp(np.clip(la, -30, 5))), m=float(np.clip(m, -12, 12)), S=float(np.exp(np.clip(lS, -6, 3))), c=float(1 / (1 + np.exp(-c))),
                       kv=float(np.clip(kv, -8, 8)), phi=float(ph), n=1))
    return ts, float(np.exp(np.clip(z[-1], np.log(0.3), np.log(4.0))))


def main():
    t0 = time.time(); out = {'what': __doc__.split('\n')[0], 'valley': [], 'fringe': [], 'random': [], 'env': C.env()}
    for R in (2.0, 3.0, 4.0, 5.0, 6.0, 8.0):
        for rs in (0.01, 0.1, 1.0, 10.0):
            for S in (0.25, 1.0):
                for kv in (0.0, 2.0, 4.0, 6.0):
                    ts, kw = valley(R, rs, S, kv); r = evaluate(ts, kw)
                    if r: r.update({'R': R, 'rs': rs, 'S': S, 'kv': kv}); out['valley'].append(r)
        v = [r for r in out['valley'] if r['R'] == R]; b = min(v, key=lambda r: r['T1_margin'])
        print('valley R %.0f: %d states, min m1 %+.3f (ratio %.3f, A %.2e, rs %g, S %g, kv %g); max D/(A ln(1/A)) %.3f  %.0fs' % (
            R, len(v), b['T1_margin'], b['ratio'], b['A'], b['rs'], b['S'], b['kv'], max(r['D'] / (r['A'] * r['ln_1_over_A']) for r in v), time.time() - t0), flush=True)
    for d in (0.5, 1.0, 2.0, 3.0):
        for kv in (1.0, 2.0, 3.0, 4.0, 5.0, 6.0):
            for E in (3.0, 10.0):
                ts, kw = fringe(d, kv, E); r = evaluate(ts, kw)
                if r: r.update({'d': d, 'kv': kv, 'E': E}); out['fringe'].append(r)
    b = min(out['fringe'], key=lambda r: r['T1_margin'])
    print('fringe: %d states, min m1 %+.3f (ratio %.3f, A %.2e, d %g, kv %g, E %g)  %.0fs' % (len(out['fringe']), b['T1_margin'], b['ratio'], b['A'], b['d'], b['kv'], b['E'], time.time() - t0), flush=True)
    # random + local optimisation
    rng = np.random.default_rng(11)
    def run_opt(J, n_start, plant=None, maxfev=250):
        best = None
        def f(z):
            ts, kw = rand_terms(z, J)
            try: r = evaluate(ts, kw, nw=32, dv=0.08)
            except (FloatingPointError, ValueError, MemoryError): return 1e3
            if r is None or not np.isfinite(r['T1_margin']): return 1e3
            val = r['T1_margin']
            if plant is not None: val -= 5.0 * np.exp(-0.5 * np.sum(((z - plant) / 1.0) ** 2))
            return val
        for _ in range(n_start):
            z0 = np.concatenate([np.column_stack([rng.normal(-1, 2, J), rng.normal(0, 3, J), rng.normal(-1, 1.5, J), rng.normal(1, 2, J), rng.normal(0, 3, J),
                                                  rng.uniform(0, 6.28, J)]).ravel(), [rng.normal(0, 0.5)]])
            r = minimize(f, z0, method='Nelder-Mead', options={'maxfev': maxfev, 'xatol': 1e-3, 'fatol': 1e-5})
            if best is None or r.fun < best.fun: best = r
        return best
    for J in (2, 3, 4):
        b = run_opt(J, 10); ts, kw = rand_terms(b.x, J); r = evaluate(ts, kw)
        r.update({'J': J, 'terms': ts, 'kw': kw, 'fast_objective': float(b.fun)}); out['random'].append(r)
        print('random J %d: best m1 %+.4f (ratio %.3f, A %.2e)  %.0fs' % (J, r['T1_margin'], r['ratio'], r['A'], time.time() - t0), flush=True)
    J = 2; z_pl = np.concatenate([np.column_stack([rng.normal(-1, 2, J), rng.normal(0, 3, J), rng.normal(-1, 1.5, J), rng.normal(1, 2, J), rng.normal(0, 3, J),
                                                   rng.uniform(0, 6.28, J)]).ravel(), [0.0]])
    b = run_opt(J, 10, plant=z_pl); dist = float(np.linalg.norm(b.x - z_pl))
    out['plant'] = {'found': bool(dist < 2.0), 'distance': dist, 'objective': float(b.fun), 'width': 1.0, 'dimension': int(len(z_pl))}
    allr = out['valley'] + out['fringe'] + out['random']
    out['n_states'] = len(allr); out['min_T1_margin'] = min(r['T1_margin'] for r in allr); out['max_ratio'] = max(r['ratio'] for r in allr)
    out['n_half_broken'] = sum(r['T1_half_margin'] < 0 for r in allr)
    out['max_D_over_A_lnA'] = max(r['D'] / (r['A'] * r['ln_1_over_A']) for r in allr)
    out['ok'] = True; out['seconds'] = time.time() - t0
    C.save('ed_vmix.json', out)
    print('states %d: min T1 margin %+.4f, max ratio %.4f, max D/(A ln(1/A)) %.3f, (1/2)-form broken on %d; plant found %s (dist %.2f)' % (
        out['n_states'], out['min_T1_margin'], out['max_ratio'], out['max_D_over_A_lnA'], out['n_half_broken'], out['plant']['found'], dist))


if __name__ == '__main__':
    np.seterr(all='ignore'); main()
