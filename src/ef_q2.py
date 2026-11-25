"""Stage EF, 1.2: the lemmas on Q2 = E[(E[Z_v Z_w | Y])^2] = int (d_v d_w p)^2/p along the v-heat flow, checked on the recorded states.
sigma = 1 units; Y_t = (X_v + sqrt(t) Z_v, X_w + Z_w), p_t its density; F(t) = E[s_w^2] (F(1) = A), D(t) = E[(d_v s_w)^2] = -F'(t),
Q2(t) = int (d_v d_w p_t)^2/p_t, SVm(t) = t E[s_w^2 d_v^2 ln p_t] (= SV - A at t = 1), S2 = E[s_v^2 s_w^2].
Checked statements (report, section 3):
  Q-id   D(t) = Q2(t) + SVm(t)/t                                   (identity; at t = 1: D = Q2 + SV - A)
  Q-mon  Q2(t) is non-increasing in t                               (Q2(1 - s) >= Q2(1))
  Q-sm   Q2(1) <= F(1 - s)/s            for every s in (0, 1)       (smoothing)
  Q-sm2  Q2(1) <= (F(1 - s) - A)/s + F(1 - s)/(1 - s)               (sharper form)
  Q-low  F(t) >= t (A + Q2(1) ln(1/t))                              (lower bound along the flow)
  stab   rho* := F(1 - s*)/A with s* = 1/(1 + ln(1/A))              (the stability ratio: T1 <=> rho* bounded, Theorem EF-S)
  T1-s   D <= F(1 - s*)/s* + SV - A                                 (what Q-sm gives for D)
2-D families (stage ED's records, closed forms along the flow): tilt, torus, lattice, vmix, atoms.
1-D families (stage EE's ee_core1d records; the same statements for X real, m in [-1, 1]: F = E[h^2], Q = int (q')^2/p, D1 = E[h'^2]):
  valley_ripple, lattice, valley_pair, nested, valley_coherent, valley_cos, random.
usage: python ef_q2.py [2d|1d] -> results_EF/ef_q2_<part>.json.  Python 3.11 grammar."""
import os, sys, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.special import logsumexp
import ed_common as C

SS = (0.1, 0.5)                      # besides s* = 1/(1 + ln(1/A))


def load(res, n):
    return json.load(open(os.path.join(res, n), encoding='utf-8'))


# ------------------------------------------------------------------ 2-D: sums at blur t from p and its derivatives on a grid
def grid_sums(p, pv, pw, pvv, pvw, t):
    ok = p > 1e-290; ps = np.where(ok, p, 1.0); m = np.sum(np.where(ok, p, 0.0)); I = lambda x: float(np.sum(np.where(ok, x, 0.0)) / m)
    sv, sw = pv / ps, pw / ps; hv = pvw / ps - sv * sw
    return {'F': I(p * sw * sw), 'D': I(p * hv * hv), 'Q2': I(pvw * pvw / ps), 'SVm': I(p * sw * sw * t * (pvv / ps - sv * sv)), 'S2': I(p * sv * sv * sw * sw)}


def tilt_sums(r, t, nth=2048):
    b, L, kv, kw = r['b'], r['L'], r['kv'], r['kw']; harm = r['harm']
    V0 = 1.0 / b - 1; ktw = (1 - 1 / L ** 2) * kw; ktv = kv * V0 / (V0 + t); th = (np.arange(nth * len(harm)) + 0.5) * 2 * np.pi / (nth * len(harm))
    R = np.zeros_like(th); R1 = np.zeros_like(th); R2 = np.zeros_like(th)
    for n, (c, ph) in enumerate(harm, 1):
        e = c * np.exp(-n * n * (kv ** 2 * V0 * t / (2 * (V0 + t)) + ktw ** 2 / (2 * (1 - 1 / L ** 2))))
        R += e * np.cos(n * th + ph); R1 += -n * e * np.sin(n * th + ph); R2 += -n * n * e * np.cos(n * th + ph)
    w = (1 + R) / np.sum(1 + R); l1 = R1 / (1 + R); l2 = R2 / (1 + R) - l1 ** 2
    # s_v = -v/(V0 + t) + ktv l1, s_w = -w/L^2 + ktw l1, d_v s_w = ktv ktw l2, d_v s_v = -1/(V0 + t) + ktv^2 l2; v, w Gaussian and independent of the phase
    gv, gw = 1 / (V0 + t), 1 / L ** 2; al, be, cc = ktv * l1, ktw * l1, ktv * ktw * l2; E = lambda x: float(np.sum(w * x))
    return {'F': E(be ** 2) + gw, 'D': E(cc ** 2), 'Q2': E(cc ** 2 + 2 * cc * al * be + (al ** 2 + gv) * (be ** 2 + gw)),
            'SVm': E((be ** 2 + gw) * t * (-gv + ktv ** 2 * l2)), 'S2': E((al ** 2 + gv) * (be ** 2 + gw))}


def torus_sums(rc, kv, kw, t, n=4096):
    N = len(rc) - 1; th = (np.arange(n) + 0.5) * 2 * np.pi / n; dmp = np.exp(-0.5 * np.arange(N + 1) ** 2 * (kv ** 2 * t + kw ** 2))
    P = np.full(n, rc[0].real * dmp[0]); P1 = np.zeros(n); P2 = np.zeros(n)
    for j in range(1, N + 1):
        cj = rc[j] * dmp[j]; e = np.exp(1j * j * th)
        P += 2 * np.real(cj * e); P1 += 2 * np.real(1j * j * cj * e); P2 += 2 * np.real(-j * j * cj * e)
    if P.min() <= 0: return None
    w = P / np.sum(P); l1 = P1 / P; l2 = P2 / P - l1 ** 2; E = lambda x: float(np.sum(w * x))
    return {'F': kw ** 2 * E(l1 ** 2), 'D': kv ** 2 * kw ** 2 * E(l2 ** 2), 'Q2': kv ** 2 * kw ** 2 * E((P2 / P) ** 2),
            'SVm': kw ** 2 * t * kv ** 2 * E(l1 ** 2 * l2), 'S2': kv ** 2 * kw ** 2 * E(l1 ** 4)}


def lattice_sums(r, t, N=24, n=128):
    kap, rr, al, be = r['kap'], r['r'], r['alpha'], r['beta']
    g1 = kap * np.array([np.cos(al), np.sin(al)]); g2 = kap * rr * np.array([np.cos(be), np.sin(be)])
    k = np.fft.fftfreq(n, 1.0 / n); n1, n2 = np.meshgrid(k, k, indexing='ij'); keep = (np.abs(n1) <= N) & (np.abs(n2) <= N)
    Gv = n1 * g1[0] + n2 * g2[0]; Gw = n1 * g1[1] + n2 * g2[1]; c = np.where(keep, np.exp(-0.5 * (t * Gv ** 2 + Gw ** 2)), 0.0)
    if c[keep & ((np.abs(n1) == N) | (np.abs(n2) == N))].max() > 1e-14: return None
    f = lambda m: np.real(np.fft.ifft2(c * m)) * n * n
    p = f(1.0)
    if p.min() <= 0: return None
    return grid_sums(p, f(1j * Gv), f(1j * Gw), f(-Gv * Gv), f(-Gv * Gw), t)


def vmix_sums(terms, kw, t, nw=64):
    T = np.array([x['S'] + t for x in terms]); m = np.array([x['m'] for x in terms])
    lo = float(np.min(m - 9.5 * np.sqrt(T))); hi = float(np.max(m + 9.5 * np.sqrt(T)))
    qmax = max([abs(x.get('kv', 0.0)) * x['S'] / (x['S'] + t) for x in terms] + [1.0]); nmax = max(x.get('n', 1) for x in terms)
    dv = max(min(0.04 * np.sqrt(t), 2 * np.pi / (24 * qmax), 0.04), (hi - lo) / 20000.0)
    v = np.arange(lo, hi + dv, dv); nw = nw * nmax; w = (np.arange(nw) + 0.5) * (2 * np.pi / kw) / nw
    Vg, Wg = np.meshgrid(v, w, indexing='ij')
    p = np.zeros_like(Vg); pv = np.zeros_like(Vg); pw = np.zeros_like(Vg); pvw = np.zeros_like(Vg); pvv = np.zeros_like(Vg)
    for x in terms:
        Tj = x['S'] + t; u = Vg - x['m']; g = x['a'] * np.exp(-0.5 * u * u / Tj) / np.sqrt(2 * np.pi * Tj); g1 = -u / Tj * g; g2 = (u * u / Tj ** 2 - 1 / Tj) * g
        c = x.get('c', 0.0)
        if c == 0:
            p += g; pv += g1; pvv += g2; continue
        kv = x.get('kv', 0.0); nn = x.get('n', 1); q = kv * x['S'] / Tj; d = c * np.exp(-0.5 * (kv ** 2 * x['S'] * t / Tj + (nn * kw) ** 2))
        Th = q * u + kv * x['m'] + nn * kw * Wg + x.get('phi', 0.0); Cc = d * np.cos(Th); Ss = d * np.sin(Th)
        p += g * (1 + Cc); pv += g1 * (1 + Cc) - g * q * Ss; pw += -g * nn * kw * Ss
        pvw += -g1 * nn * kw * Ss - g * q * nn * kw * Cc; pvv += g2 * (1 + Cc) - 2 * g1 * q * Ss - g * q * q * Cc
    return grid_sums(p, pv, pw, pvv, pvw, t)


def atoms_sums(blobs, t, n=401):
    lo_v = min(b[1] - 8 * np.sqrt(b[3] + t) for b in blobs); hi_v = max(b[1] + 8 * np.sqrt(b[3] + t) for b in blobs)
    lo_w = min(b[2] - 8 * np.sqrt(b[4] + 1) for b in blobs); hi_w = max(b[2] + 8 * np.sqrt(b[4] + 1) for b in blobs)
    v = np.linspace(lo_v, hi_v, 2 * n); w = np.linspace(lo_w, hi_w, n); V, W = np.meshgrid(v, w, indexing='ij')
    p = np.zeros_like(V); pv = np.zeros_like(V); pw = np.zeros_like(V); pvw = np.zeros_like(V); pvv = np.zeros_like(V)
    for a, mv, mw, Sv, Sw in blobs:
        Tv, Tw = Sv + t, Sw + 1; u, x = V - mv, W - mw
        g = a * np.exp(-0.5 * (u * u / Tv + x * x / Tw)) / (2 * np.pi * np.sqrt(Tv * Tw))
        p += g; pv += -u / Tv * g; pw += -x / Tw * g; pvw += u * x / (Tv * Tw) * g; pvv += (u * u / Tv ** 2 - 1 / Tv) * g
    return grid_sums(p, pv, pw, pvv, pvw, t)


# ------------------------------------------------------------------ 1-D: X real, m in [-1, 1]; dimensionless sums at noise t
def comps_sums(comps, t, n=20001):
    comps = [c for c in comps if c[0] > 0 or len(c) > 4]
    lo = min(c[1] - 9 * np.sqrt(c[2] + 1) for c in comps); hi = max(c[1] + 9 * np.sqrt(c[2] + 1) for c in comps)
    y = np.linspace(lo, hi, n); dy = y[1] - y[0]
    lw = np.array([c[4] if len(c) > 4 else np.log(c[0]) for c in comps]); lw -= logsumexp(lw)
    mu = np.array([c[1] for c in comps]); T = np.array([c[2] for c in comps]) + t; mm = np.array([c[3] for c in comps])
    L = lw[:, None] - 0.5 * (y[None, :] - mu[:, None]) ** 2 / T[:, None] - 0.5 * np.log(2 * np.pi * T[:, None])
    lp = logsumexp(L, axis=0); pi = np.exp(L - lp[None, :]); d1 = -(y[None, :] - mu[:, None]) / T[:, None]
    h = np.sum(pi * mm[:, None], 0); s = np.sum(pi * d1, 0); qp = np.sum(pi * mm[:, None] * d1, 0); hp = qp - h * s
    s1 = np.sum(pi * (d1 * d1 - 1 / T[:, None]), 0) - s * s                                  # (ln p)''
    return _logsums(lp, dy, h, qp, hp, s1, t)


def _logsums(lp, dy, h, qp, hp, s1, t):
    with np.errstate(divide='ignore'):
        lh = 2 * np.log(np.abs(h)); lF = float(logsumexp(lp + lh) + np.log(dy)); lQ = float(logsumexp(lp + 2 * np.log(np.abs(qp))) + np.log(dy))
        lD = float(logsumexp(lp + 2 * np.log(np.abs(hp))) + np.log(dy))
        pos = np.maximum(s1, 0.0); neg = np.maximum(-s1, 0.0)              # no floor inside the logarithm (a floor of 1e-300 at the atoms swamps deep valleys)
        lSp = float(logsumexp(lp + lh + np.log(pos)) + np.log(dy)); lSn = float(logsumexp(lp + lh + np.log(neg)) + np.log(dy))
    if not np.isfinite(lF): return None
    # everything relative to F(t): Q2 = t Q (dimensionless at noise t would be t Q; here the raw y-derivative is kept, as in the 2-D sums)
    return {'lnF': lF, 'Q_over_F': float(np.exp(lQ - lF)), 'D_over_F': float(np.exp(lD - lF)), 'SVm_over_F': float(t * (np.exp(lSp - lF) - np.exp(lSn - lF)))}


def cos_sums(R, rs, d, E, t):
    k = np.pi / d; T = E * E + t; w = k * E * E / T; lvis = -k * k * E * E * t / (2 * T)
    lo, hi = -(R + 10 * np.sqrt(T + 1)), R + 10 * np.sqrt(T + 1)
    dy = min(1.0 / (30 * max(R, 1.0)), 2 * np.pi / (60 * w), 0.02); n = int(min(2.5e6, (hi - lo) / dy) + 1); y = np.linspace(lo, hi, n); dy = y[1] - y[0]
    comps = [np.log(rs) - R * R / 2 - 0.5 * y * y / T - 0.5 * np.log(2 * np.pi * T)]; d1 = [-y / T]; Ts = [T]
    if R > 0:
        for c in (-R, R): comps.append(np.log(0.5) - 0.5 * (y - c) ** 2 / t - 0.5 * np.log(2 * np.pi * t)); d1.append(-(y - c) / t); Ts.append(t)
    L = np.stack(comps); lp = logsumexp(L, axis=0); pi = np.exp(L - lp[None, :]); d1 = np.stack(d1); s = np.sum(pi * d1, 0)
    s1 = np.sum(pi * (d1 * d1 - 1 / np.array(Ts)[:, None]), 0) - s * s
    vis = np.exp(lvis - (-k * k * E * E * 1.0 / (2 * (E * E + 1))))          # relative to the visibility at t = 1 (keeps the numbers in range)
    cs, sn = np.cos(w * y), np.sin(w * y); h = pi[0] * vis * cs; qp = pi[0] * vis * (-(y / T) * cs - w * sn); hp = qp - h * s
    o = _logsums(lp, dy, h, qp, hp, s1, t)
    if o: o['lnF'] += 2 * (-k * k * E * E / (2 * (E * E + 1)))
    return o


def checks(get, A_of=None):
    """get(t) -> sums; returns the row of checks for one state (2-D form: keys F, D, Q2, SVm)."""
    s1 = get(1.0)
    if s1 is None or not (0 < s1['F'] < 1): return None
    A, D, Q2, SVm = s1['F'], s1['D'], s1['Q2'], s1['SVm']; L = np.log(1 / A); ss = 1 / (1 + L)
    row = {'A': A, 'D': D, 'Q2': Q2, 'SV': SVm + A, 'ln_1_over_A': float(L), 's_star': float(ss),
           'id_rel_residual': float(abs(D - Q2 - SVm) / max(D, Q2, A, 1e-300)), 'D_over_A_lnA': float(D / (A * max(L, 1e-12))), 'Q2_over_A_lnA': float(Q2 / (A * max(L, 1e-12)))}
    for s in (ss,) + SS:
        o = get(1 - s)
        if o is None: continue
        tag = 'star' if s == ss else ('%g' % s); Ft = o['F']; t = 1 - s
        row['sm_margin_' + tag] = float((Ft / s - Q2) / max(Q2, 1e-300)); row['sm2_margin_' + tag] = float(((Ft - A) / s + Ft / t - Q2) / max(Q2, 1e-300))
        row['mon_margin_' + tag] = float((o['Q2'] - Q2) / max(Q2, 1e-300)); row['low_margin_' + tag] = float((Ft - t * (A + Q2 * np.log(1 / t))) / Ft)
        if s == ss:
            row['rho_star'] = float(Ft / A); row['T1s_margin'] = float((Ft / s + SVm - D) / max(D, 1e-300)); row['rho_bound_ratio'] = float(Ft / (s * A * (1 + L)))
    return row


def checks_1d(get):
    """1-D form: get(t) -> {'lnF', 'Q_over_F', 'D_over_F', 'SVm_over_F'}; the same statements (F' = -D1, D1 = Q + SVm/t)."""
    s1 = get(1.0)
    if s1 is None or not (s1['lnF'] < 0): return None
    L = -s1['lnF']; ss = 1 / (1 + L); Q, D = s1['Q_over_F'], s1['D_over_F']          # relative to A = F(1)
    row = {'ln_A': float(s1['lnF']), 'ln_1_over_A': float(L), 's_star': float(ss), 'D_over_A': D, 'Q_over_A': Q, 'SVm_over_A': s1['SVm_over_F'],
           'id_rel_residual': float(abs(D - Q - s1['SVm_over_F']) / max(D, Q, 1.0)), 'D_over_A_lnA': float(D / max(L, 1e-12)), 'Q_over_A_lnA': float(Q / max(L, 1e-12))}
    for s in (ss,) + SS:
        o = get(1 - s)
        if o is None: continue
        tag = 'star' if s == ss else ('%g' % s); rho = float(np.exp(min(o['lnF'] - s1['lnF'], 700.0))); t = 1 - s
        row['sm_margin_' + tag] = float((rho / s - Q) / max(Q, 1e-300)); row['sm2_margin_' + tag] = float(((rho - 1) / s + rho / t - Q) / max(Q, 1e-300))
        row['mon_margin_' + tag] = float((o['Q_over_F'] * rho - Q) / max(Q, 1e-300)); row['low_margin_' + tag] = float((rho - t * (1 + Q * np.log(1 / t))) / rho)
        if s == ss: row['rho_star'] = rho; row['T1s_margin'] = float((rho / s + s1['SVm_over_F'] - D) / max(D, 1e-300))
    return row


def summarise(fams):
    out = {}
    for f, rows in fams.items():
        if not rows: continue
        g = lambda k, fn: float(fn([r[k] for r in rows if k in r and np.isfinite(r[k])] or [np.nan]))
        small = [r for r in rows if r['ln_1_over_A'] > np.log(1e3)]
        out[f] = {'n': len(rows), 'n_A_below_1e-3': len(small), 'id_max_rel_residual': g('id_rel_residual', max),
                  'sm_min_margin': min(g('sm_margin_' + t, min) for t in ('star', '0.1', '0.5')), 'sm2_min_margin': min(g('sm2_margin_' + t, min) for t in ('star', '0.1', '0.5')),
                  'mon_min_margin': min(g('mon_margin_' + t, min) for t in ('star', '0.1', '0.5')), 'low_min_margin': min(g('low_margin_' + t, min) for t in ('star', '0.1', '0.5')),
                  'T1s_min_margin': g('T1s_margin', min), 'rho_star_max': g('rho_star', max),
                  'rho_star_max_A_below_1e-3': float(max([r['rho_star'] for r in small if 'rho_star' in r] or [np.nan])),
                  'Q_over_A_lnA_max_A_below_1e-3': float(max([r.get('Q2_over_A_lnA', r.get('Q_over_A_lnA')) for r in small] or [np.nan])),
                  'D_over_A_lnA_max_A_below_1e-3': float(max([r['D_over_A_lnA'] for r in small] or [np.nan]))}
    for v in out.values():                                              # no NaN in the records (the published comparison treats nan != nan): None instead
        for k in list(v):
            if isinstance(v[k], float) and not np.isfinite(v[k]): v[k] = None
    return out


def part_2d():
    t0 = time.time(); fams = {}
    def add(name, items, fn):
        rows = []
        for it in items:
            try: r = checks(lambda t: fn(it, t))
            except (ValueError, FloatingPointError, MemoryError): r = None
            if r: rows.append(r)
        fams[name] = rows; print('%-10s %5d states  %.0fs' % (name, len(rows), time.time() - t0), flush=True)
    add('tilt', load(C.RES, 'ed_scan_tilt.json')['rows'], tilt_sums)
    import ed_torus as To
    tor = load(C.RES, 'ed_torus.json')
    ts = [(np.array([complex(a, b) for a, b in r['rho_coeffs']]), r['kv'], r['kw']) for r in tor['optimised']] + [(np.array([1.0, 0.5], complex), r['kv'], r['kw']) for r in tor['single']] \
        + [(To.comb_coeffs(r['kap'] ** 2), r['kv'], r['kw']) for r in tor['comb']]
    add('torus', ts, lambda s, t: torus_sums(s[0], s[1], s[2], t))
    lat = load(C.RES, 'ed_lattice.json'); add('lattice', lat['rows'] + lat['optimised'], lattice_sums)
    import ed_vmix as V
    vm = load(C.RES, 'ed_vmix.json'); vs = [V.valley(r['R'], r['rs'], r['S'], r['kv']) for r in vm['valley']] + [V.fringe(r['d'], r['kv'], r['E']) for r in vm['fringe'][::2]] \
        + [(r['terms'], r['kw']) for r in vm['random']]
    add('vmix', vs, lambda s, t: vmix_sums(s[0], s[1], t))
    at = load(C.RES, 'ed_atoms.json'); add('atoms', [[tuple(x) for x in r['blobs']] for r in at['best_excess'] + at['best_D_ratio']], atoms_sums)
    return finish('2d', fams, t0)


def part_1d():
    import ee_core1d as K1
    t0 = time.time(); fams = {}; d = load(C.RES_EE, 'ee_core1d.json')['families']
    B = {'valley_ripple': lambda r: K1.valley_ripple(r['R'], r['rs'], r['d'], r['E']), 'lattice': lambda r: K1.lattice(r['d'], r['E']),
         'valley_pair': lambda r: K1.valley_pair(r['R'], r['rs'], r['d']), 'nested': lambda r: K1.nested(r['R1'], r['R2'], r['r2'], r['rs'], r['d'], 1.5),
         'valley_coherent': lambda r: [(0.5, -r['R'], 0.0, 0.0), (0.5, r['R'], 0.0, 0.0), (0.0, 0.0, 0.0, 1.0, np.log(r['rs']) - r['R'] ** 2 / 2)],
         'random': lambda r: [tuple(c) for c in r['comps']]}
    for f, rows in d.items():
        out = []
        for r in rows:
            if f == 'valley_ripple' and 'search' in r: continue
            try:
                if f == 'valley_cos': o = checks_1d(lambda t: cos_sums(r['R'], r['rs'], r['d'], r['E'], t))
                else:
                    cs = B[f](r); n = int(max(20001, 400 * r.get('R', r.get('R1', 1.0)))); o = checks_1d(lambda t: comps_sums(cs, t, n))
            except (ValueError, FloatingPointError, MemoryError, KeyError): o = None
            if o: out.append(o)
        fams[f] = out; print('%-16s %5d states  %.0fs' % (f, len(out), time.time() - t0), flush=True)
    return finish('1d', fams, t0)


def finish(part, fams, t0):
    out = {'part': part, 'families': summarise(fams), 'rows': fams, 'n_states': sum(len(v) for v in fams.values()), 'env': C.env()}
    S = out['families']
    out['all'] = {k: (min if 'min' in k else max)(v[k] for v in S.values() if v[k] is not None) for k in ('id_max_rel_residual', 'sm_min_margin', 'sm2_min_margin', 'mon_min_margin',
                  'low_min_margin', 'T1s_min_margin', 'rho_star_max', 'rho_star_max_A_below_1e-3', 'Q_over_A_lnA_max_A_below_1e-3', 'D_over_A_lnA_max_A_below_1e-3')}
    a = out['all']; tol = 2e-3
    out['holds'] = {'Q-id': bool(a['id_max_rel_residual'] < tol), 'Q-sm': bool(a['sm_min_margin'] > -tol), 'Q-sm2': bool(a['sm2_min_margin'] > -tol), 'Q-mon': bool(a['mon_min_margin'] > -tol),
                    'Q-low': bool(a['low_min_margin'] > -tol), 'T1-s': bool(a['T1s_min_margin'] > -tol)}
    out['ok'] = True; out['seconds'] = time.time() - t0
    C.save('ef_q2_%s.json' % part, out, C.RES_EF)
    nn = lambda x: float('nan') if x is None else x
    for f, v in S.items(): print('%-16s n %5d: id %.1e  sm %+.3g  sm2 %+.3g  mon %+.3g  low %+.3g  T1s %+.3g  rho* max %.3f (A<1e-3: %.3f)  Q2/(A ln) %.3f  D/(A ln) %.3f' % (
        f, v['n'], v['id_max_rel_residual'], v['sm_min_margin'], v['sm2_min_margin'], v['mon_min_margin'], v['low_min_margin'], v['T1s_min_margin'], v['rho_star_max'],
        nn(v['rho_star_max_A_below_1e-3']), nn(v['Q_over_A_lnA_max_A_below_1e-3']), nn(v['D_over_A_lnA_max_A_below_1e-3'])))
    print('all: %s  holds: %s' % ({k: round(v, 4) for k, v in a.items()}, out['holds']))
    return out


if __name__ == '__main__':
    np.seterr(all='ignore')
    part = sys.argv[1]; fn = os.path.join(C.RES_EF, 'ef_q2_%s.json' % part)
    if os.path.exists(fn) and os.environ.get('EF_FORCE') != '1': print('exists, skipped'); sys.exit(0)
    {'2d': part_2d, '1d': part_1d}[part]()
