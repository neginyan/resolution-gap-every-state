"""Stage EG, 1.1: every step of the three-circle proof (report EG, Theorem EG-1) on the recorded states.
sigma = 1.  p_z = the density with v-noise variance z (complex, |z - T| < T), q_z = d_w p_z (one-variable core: q_z = (m mu) * phi_z, p_z = mu * phi_z).
H_T = L^2(dy/p_T),  N(z) := ||q_z||^2_{H_T} = int |q_z|^2 / p_T.   A = F(1) = N at z = T = 1 ... (F(t) = int q_t^2/p_t),  L = ln(1/A) >= 4,
  delta = 1/(4L), T = 1 + delta, eta = ln(2L)/(4L), s = ln(2L)/(8L^2), r1 = delta (1 - eta), r2 = delta + s, r3 = T/2, theta = ln(r2/r1)/ln(r3/r1).
Checked per state (margins = ln(right side) - ln(left side) >= 0):
  outer   N(T + r3 e^{i a}) <= 1/sqrt(1 - (r3/T)^2) = 2/sqrt3                       (a = 0, 45, 90, 135, 180 degrees; N(conj z) = N(z))
  inner   N(T + r1 e^{i a}) <= A/sqrt(1 - (r1/delta)^2)   (<= A/sqrt(eta))
  three   N(1 - s) <= (A/sqrt(eta))^(1 - theta) (2/sqrt3)^theta  <= (4e/sqrt3) sqrt(L/ln 2L) A          and theta < 1/L
  ident   p_{1-s} = sqrt(T/(1-s)) p_T W, W = E[exp(-alpha Z_v^2) | Y_T]  (W from this identity; the identity itself is checked independently in eg_point.py)
  bad     P_T(W < 1/(2e)) <= 2 A^2;   int_bad q_{1-s}^2/p_{1-s} <= sqrt(6T/(1-s)) A   (one-variable core: <= sqrt(T/(1-s)) 2 A^2)
  good    int_good q_{1-s}^2/p_{1-s} <= 2e sqrt((1-s)/T) N(1 - s)
  F       F(1 - s) <= 36.1 sqrt(L/ln 2L) A
  Q       Q2 <= F(1 - s)/s  (Lemma Q-sm)  <= 289 A L^(5/2)/(ln 2L)^(3/2)
Isotropic form (2-D families, both noise variances = z; report, Lemma EG-iso): outer <= T^3/(T^2 - r3^2)^2 = 16/(9T), inner <= A/(1 - (r1/delta)^2),
  three <= 12.62 (L/ln 2L) A, bad <= sqrt6 T^2/(1-s)^3 A, good <= 2e ((1-s)/T) N, G(1 - s) <= 70.2 (L/ln 2L) A, int (d_w^2 p)^2/p <= G(1-s)/s <= 562 A L^3/(ln 2L)^2.
Tool A of the instruction (one-variable core): J_S(t) = E[((1+h)/2) ln(1+h) + ((1-h)/2) ln(1-h)], the ratio J_S(1 - s*)/J_S(1), s* = 1/(1 + L), next to rho* = F(1 - s*)/F(1).
usage: python eg_circ.py [1d|1ds|2d|iso] -> results_EG/eg_circ_<part>.json.  Python 3.11 grammar."""
import os, sys, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.special import logsumexp
import ed_common as C

ANG = (0.0, 0.25 * np.pi, 0.5 * np.pi, 0.75 * np.pi, np.pi)
LN = np.log


def params(L):
    d = 1 / (4 * L); T = 1 + d; eta = LN(2 * L) / (4 * L); s = LN(2 * L) / (8 * L * L); r1 = d * (1 - eta); r2 = d + s; r3 = T / 2
    return {'delta': d, 'T': T, 'eta': eta, 's': s, 'r1': r1, 'r2': r2, 'r3': r3, 'theta': LN(r2 / r1) / LN(r3 / r1), 'alpha': (s + d) / (2 * (1 - s))}


# ------------------------------------------------------------------ one-variable core, log domain
def terms_1d(comps):
    """comps (w, mu, S, m[, logw]) -> terms (logw, mu, S, kappa, cp, cq)"""
    out = []
    for c in comps:
        if not (c[0] > 0 or len(c) > 4): continue
        out.append((c[4] if len(c) > 4 else LN(c[0]), c[1], c[2], 0.0, 1.0, c[3]))
    lw = np.array([t[0] for t in out]); lw -= logsumexp(lw)
    return [(float(l),) + t[1:] for l, t in zip(lw, out)]


def terms_cos(R, rs, d, E):
    """valley_cos: centre N(0, E^2) of weight rs e^{-R^2/2} with m = cos(pi x/d), atoms at -R, R of weight 1/2 with m = 0 (not normalised: as ef_q2.cos_sums)"""
    k = np.pi / d; l0 = LN(rs) - R * R / 2; t = [(l0, 0.0, E * E, 0.0, 1.0, 0.0), (l0, 0.0, E * E, k, 0.0, 0.5), (l0, 0.0, E * E, -k, 0.0, 0.5)]
    if R > 0: t += [(LN(0.5), -R, 0.0, 0.0, 1.0, 0.0), (LN(0.5), R, 0.0, 0.0, 1.0, 0.0)]
    return t


def lval(terms, z, y, which):
    """log of |sum_i c_i term_i(y)| and its phase, for which = 4 (p) or 5 (q); z complex.  In pieces of y (a lattice has thousands of atoms)."""
    n = sum(1 for t in terms if t[which] != 0); step = max(2000, int(3e6 // max(n, 1)))
    if len(y) > step:
        parts = [_lval(terms, z, y[i:i + step], which) for i in range(0, len(y), step)]
        return np.concatenate([p[0] for p in parts]), np.concatenate([p[1] for p in parts])
    return _lval(terms, z, y, which)


def _lval(terms, z, y, which):
    Ls = []; cs = []
    for t in terms:
        c = t[which]
        if c == 0: continue
        lw, mu, S, kap = t[0], t[1], t[2], t[3]; Tz = S + z; u = y - mu
        Ls.append(lw - u * u / (2 * Tz) - 0.5 * np.log(2 * np.pi * Tz) + 1j * kap * (mu + u * S / Tz) - kap * kap * S * z / (2 * Tz)); cs.append(c)
    if not Ls: return np.full(len(y), -np.inf), np.zeros(len(y), complex)
    Lm = np.stack(Ls); M = Lm.real.max(0); M = np.where(np.isfinite(M), M, -1e300); v = np.sum(np.array(cs)[:, None] * np.exp(Lm - M[None, :]), 0)
    with np.errstate(divide='ignore'): return M + np.log(np.abs(v)), v / np.maximum(np.abs(v), 1e-300)


def grid_1d(terms, n_min=20001, valley=True):
    """valley=True: step <= 1/(30 R), R the largest |centre| (between two far heavy atoms the posterior weights switch over a width ~ 1/R);
    valley=False (the plain lattices: equally spaced atoms, no valley): step <= 0.02 only."""
    lo = min(t[1] - 9.5 * np.sqrt(t[2] + 1.3) for t in terms); hi = max(t[1] + 9.5 * np.sqrt(t[2] + 1.3) for t in terms)
    R = max(abs(t[1]) for t in terms) if valley else 1.0; kmax = max(abs(t[3]) for t in terms)
    dy = min((hi - lo) / n_min, 1.0 / (30 * max(R, 1.0)), 0.02, (2 * np.pi / (60 * kmax)) if kmax > 0 else 1.0); n = int(min(2.5e6, (hi - lo) / dy)) + 1
    y = np.linspace(lo, hi, n); return y, y[1] - y[0]


def check_1d(terms, valley=True):
    y, dy = grid_1d(terms, valley=valley); ldy = LN(dy)
    lp = lambda t: lval(terms, t, y, 4)[0]
    def lN(z, lpT):                                   # ln N(z) = ln int |q_z|^2/p_T
        return float(logsumexp(2 * lval(terms, z, y, 5)[0] - lpT) + ldy)
    lp1 = lp(1.0); lmass = float(logsumexp(lp1) + ldy); lA = lN(1.0, lp1) - lmass
    if not (lA < 0): return None
    L = -lA; row = {'ln_A': lA, 'L': L}
    # tool A and rho* at s* = 1/(1 + L) (every L)
    ss = 1 / (1 + L); lps = lp(1 - ss); lFs = lN(1 - ss, lps) - lmass; row['rho_star'] = float(np.exp(min(lFs - lA, 700)))
    def lJS(t, lpt):
        lq, sg = lval(terms, t, y, 5); h = np.clip(sg.real * np.exp(np.minimum(lq - lpt, 0.0)), -1 + 1e-15, 1 - 1e-15); a = np.abs(h)
        with np.errstate(divide='ignore', invalid='ignore'):
            lf = np.where(a < 1e-3, 2 * (lq - lpt) - LN(2) + np.log1p(a * a / 6), np.log(np.maximum(0.5 * ((1 + h) * np.log1p(h) + (1 - h) * np.log1p(-h)), 1e-300)))
        return float(logsumexp(lpt + lf) + ldy) - lmass
    lJ1 = lJS(1.0, lp1); lJs = lJS(1 - ss, lps); row['JS_ratio'] = float(np.exp(min(lJs - lJ1, 700))); row['JS_over_F'] = float(np.exp(lJ1 - lA))
    if L < 4: return row
    P = params(L); T, s = P['T'], P['s']; lpT = lp(T); row.update({k: float(v) for k, v in P.items()})
    out = [lN(T + P['r3'] * np.exp(1j * a), lpT) - lmass for a in ANG]; inn = [lN(T + P['r1'] * np.exp(1j * a), lpT) - lmass for a in ANG]
    lNs = lN(1 - s, lpT) - lmass
    row['outer_margin'] = float(LN(2 / np.sqrt(3)) - max(out)); row['inner_margin'] = float(lA - 0.5 * LN(1 - (P['r1'] / P['delta']) ** 2) - max(inn))
    row['inner_sup_over_A'] = float(np.exp(max(inn) - lA)); row['outer_sup'] = float(np.exp(max(out)))
    th = P['theta']; row['theta_L'] = float(th * L)
    row['three_margin'] = float((1 - th) * (lA - 0.5 * LN(P['eta'])) + th * LN(2 / np.sqrt(3)) - lNs)
    row['three_sampled_margin'] = float((1 - th) * max(inn) + th * max(out) - lNs)            # with the sampled sups (informative only)
    row['three_final_margin'] = float(LN(4 * np.e / np.sqrt(3)) + 0.5 * LN(L / LN(2 * L)) + lA - lNs); row['N_over_A'] = float(np.exp(lNs - lA))
    # the return to F(1 - s)
    lps_ = lp(1 - s); lq, _ = lval(terms, 1 - s, y, 5); c = np.sqrt(T / (1 - s)); lW = lps_ - lpT - LN(c); bad = lW < -LN(2 * np.e)
    row['W_max'] = float(np.exp(lW.max())); row['bad_fraction_points'] = float(bad.mean())
    lPbad = float(logsumexp(np.where(bad, lpT, -np.inf)) + ldy) - lmass if bad.any() else -np.inf
    lbad = float(logsumexp(np.where(bad, 2 * lq - lps_, -np.inf)) + ldy) - lmass if bad.any() else -np.inf
    lgood = float(logsumexp(np.where(~bad, 2 * lq - lps_, -np.inf)) + ldy) - lmass
    lF = float(np.logaddexp(lbad, lgood))
    row['Pbad_margin'] = float(LN(2) + 2 * lA - lPbad) if bad.any() else None; row['bad_margin'] = float(LN(2 * c) + 2 * lA - lbad) if bad.any() else None
    row['good_margin'] = float(LN(2 * np.e / c) + lNs - lgood); row['F_margin'] = float(LN(36.1) + 0.5 * LN(L / LN(2 * L)) + lA - lF); row['F_over_A'] = float(np.exp(lF - lA))
    row['Q_over_A'] = None
    return row, (y, dy, lp1, lmass, lA, lF, s)


def lQ_1d(terms, y, dy, lp1, lmass):
    """ln int (d_y q_1)^2/p_1: the derivative term by term (real time 1); in pieces of y"""
    n = sum(1 for t in terms if t[5] != 0); step = max(2000, int(3e6 // max(n, 1)))
    return float(logsumexp([_lQ(terms, y[i:i + step], dy, lp1[i:i + step]) for i in range(0, len(y), step)])) - lmass


def _lQ(terms, y, dy, lp1):
    Ls = []; cs = []
    for t in terms:
        if t[5] == 0: continue
        lw, mu, S, kap = t[0], t[1], t[2], t[3]; Tz = S + 1.0; u = y - mu
        Ls.append(lw - u * u / (2 * Tz) - 0.5 * np.log(2 * np.pi * Tz) + 1j * kap * (mu + u * S / Tz) - kap * kap * S / (2 * Tz) + np.log((-u / Tz + 1j * kap * S / Tz).astype(complex) + 0j)); cs.append(t[5])
    Lm = np.stack(Ls); M = Lm.real.max(0); M = np.where(np.isfinite(M), M, -1e300); v = np.sum(np.array(cs)[:, None] * np.exp(Lm - M[None, :]), 0)
    with np.errstate(divide='ignore'): return float(logsumexp(2 * (M + np.log(np.abs(v))) - lp1) + LN(dy))


def part_1d(stride=1):
    import ee_core1d as K1
    t0 = time.time(); fams = {}; d = json.load(open(os.path.join(C.RES_EE, 'ee_core1d.json'), encoding='utf-8'))['families']
    B = {'valley_ripple': lambda r: K1.valley_ripple(r['R'], r['rs'], r['d'], r['E']), 'lattice': lambda r: K1.lattice(r['d'], r['E']),
         'valley_pair': lambda r: K1.valley_pair(r['R'], r['rs'], r['d']), 'nested': lambda r: K1.nested(r['R1'], r['R2'], r['r2'], r['rs'], r['d'], 1.5),
         'valley_coherent': lambda r: [(0.5, -r['R'], 0.0, 0.0), (0.5, r['R'], 0.0, 0.0), (0.0, 0.0, 0.0, 1.0, LN(r['rs']) - r['R'] ** 2 / 2)],
         'random': lambda r: [tuple(c) for c in r['comps']]}
    cdir = os.path.join(C.RES_EG, 'cache'); os.makedirs(cdir, exist_ok=True); rev = os.environ.get('EG_REVERSE') == '1'; nocache = os.environ.get('EG_NOCACHE') == '1'
    for f, rows in (list(d.items())[::-1] if rev else d.items()):
        cf = os.path.join(cdir, 'eg1d_%s_%d.json' % (f, stride))                  # every family is stored when it is done, and skipped when stored
        if os.path.exists(cf) and not nocache:
            fams[f] = json.load(open(cf)); print('%-16s %5d states  (stored)' % (f, len(fams[f])), flush=True); continue
        out = []
        for r in [x for x in rows if not (f == 'valley_ripple' and 'search' in x)][::stride]:
            try:
                tm = terms_cos(r['R'], r['rs'], r['d'], r['E']) if f == 'valley_cos' else terms_1d(B[f](r)); o = check_1d(tm, valley=(f != 'lattice'))
            except (ValueError, FloatingPointError, MemoryError, KeyError): o = None
            if o is None: continue
            if isinstance(o, tuple):
                row, (y, dy, lp1, lmass, lA, lF, s) = o; lQ = lQ_1d(tm, y, dy, lp1, lmass); L = row['L']
                row['Q_over_A'] = float(np.exp(min(lQ - lA, 700))); row['Qsm_margin'] = float(lF - LN(s) - lQ)
                row['Q_margin'] = float(LN(289.0) + lA + 2.5 * LN(L) - 1.5 * LN(LN(2 * L)) - lQ); o = row
            out.append(o)
        fams[f] = out; print('%-16s %5d states  %.0fs' % (f, len(out), time.time() - t0), flush=True)
        if not nocache:
            tmp = cf + '.tmp%d' % os.getpid(); json.dump(out, open(tmp, 'w')); os.replace(tmp, cf)
    if rev: return None
    fams = {f: fams[f] for f in d if f in fams}
    return finish('1d' if stride == 1 else '1ds', fams, t0)


# ------------------------------------------------------------------ 2-D families: fields on a grid (linear domain)
class Torus:
    def __init__(self, rc, kv, kw, n=4096): self.rc, self.kv, self.kw = rc, kv, kw; self.th = (np.arange(n) + 0.5) * 2 * np.pi / n
    def f(self, z, iso=False, second=False):
        N = len(self.rc) - 1; j = np.arange(N + 1); dmp = np.exp(-0.5 * j ** 2 * (self.kv ** 2 * z + self.kw ** 2 * (z if iso else 1.0)))
        P = np.full(len(self.th), self.rc[0].real * dmp[0], complex); P1 = np.zeros(len(self.th), complex); P2 = np.zeros(len(self.th), complex)
        for k in range(1, N + 1):
            e = np.exp(1j * k * self.th); a = self.rc[k] * dmp[k] * e; b = np.conj(self.rc[k]) * dmp[k] / e
            P += a + b; P1 += 1j * k * (a - b); P2 += -k * k * (a + b)
        return (P, self.kw * P1, self.kw ** 2 * P2) if second else (P, self.kw * P1)
    wt = None


class Lattice:
    def __init__(self, r, N=24, n=128):
        kap, rr, al, be = r['kap'], r['r'], r['alpha'], r['beta']; g1 = kap * np.array([np.cos(al), np.sin(al)]); g2 = kap * rr * np.array([np.cos(be), np.sin(be)])
        k = np.fft.fftfreq(n, 1.0 / n); n1, n2 = np.meshgrid(k, k, indexing='ij'); self.keep = (np.abs(n1) <= N) & (np.abs(n2) <= N)
        self.edge = self.keep & ((np.abs(n1) == N) | (np.abs(n2) == N)); self.Gv = n1 * g1[0] + n2 * g2[0]; self.Gw = n1 * g1[1] + n2 * g2[1]; self.n = n
    def f(self, z, iso=False, second=False):
        c = np.where(self.keep, np.exp(-0.5 * (z * self.Gv ** 2 + (z if iso else 1.0) * self.Gw ** 2)), 0.0)
        if np.abs(c[self.edge]).max() > 1e-14: raise ValueError('Fourier series not converged')
        g = lambda m: np.fft.ifft2(c * m) * self.n ** 2
        return (g(1.0), g(1j * self.Gw), g(-self.Gw ** 2)) if second else (g(1.0), g(1j * self.Gw))
    wt = None


class Vmix:
    def __init__(self, terms, kw, nw=64):
        self.terms, self.kw = terms, kw; T = np.array([x['S'] + 1.3 for x in terms]); m = np.array([x['m'] for x in terms])
        lo = float(np.min(m - 9.5 * np.sqrt(T))); hi = float(np.max(m + 9.5 * np.sqrt(T)))
        qmax = max([abs(x.get('kv', 0.0)) * x['S'] / (x['S'] + 0.5) for x in terms] + [1.0]); nmax = max(x.get('n', 1) for x in terms)
        dv = max(min(0.03, 2 * np.pi / (24 * qmax)), (hi - lo) / 20000.0); v = np.arange(lo, hi + dv, dv); nw = nw * nmax; w = (np.arange(nw) + 0.5) * (2 * np.pi / kw) / nw
        self.V, self.W = np.meshgrid(v, w, indexing='ij')
    def f(self, z, iso=False, second=False):
        p = np.zeros(self.V.shape, complex); pw = np.zeros(self.V.shape, complex); pww = np.zeros(self.V.shape, complex); om = z if iso else 1.0
        for x in self.terms:
            Tj = x['S'] + z; u = self.V - x['m']; g = x['a'] * np.exp(-0.5 * u * u / Tj) / np.sqrt(2 * np.pi * Tj); c = x.get('c', 0.0)
            if c == 0: p += g; continue
            kv = x.get('kv', 0.0); nn = x.get('n', 1); q = kv * x['S'] / Tj; d = c * np.exp(-0.5 * (kv ** 2 * x['S'] * z / Tj + (nn * self.kw) ** 2 * om))
            Th = q * u + kv * x['m'] + nn * self.kw * self.W + x.get('phi', 0.0); Cc = d * np.cos(Th); Ss = d * np.sin(Th)
            p += g * (1 + Cc); pw += -g * nn * self.kw * Ss; pww += -g * (nn * self.kw) ** 2 * Cc
        return (p, pw, pww) if second else (p, pw)
    wt = None


class Atoms:
    def __init__(self, blobs, n=401):
        self.blobs = blobs; lo_v = min(b[1] - 9 * np.sqrt(b[3] + 1.3) for b in blobs); hi_v = max(b[1] + 9 * np.sqrt(b[3] + 1.3) for b in blobs)
        lo_w = min(b[2] - 9 * np.sqrt(b[4] + 1.3) for b in blobs); hi_w = max(b[2] + 9 * np.sqrt(b[4] + 1.3) for b in blobs)
        self.V, self.W = np.meshgrid(np.linspace(lo_v, hi_v, 3 * n), np.linspace(lo_w, hi_w, n), indexing='ij')
    def f(self, z, iso=False, second=False):
        p = np.zeros(self.V.shape, complex); pw = np.zeros(self.V.shape, complex); pww = np.zeros(self.V.shape, complex)
        for a, mv, mw, Sv, Sw in self.blobs:
            Tv, Tw = Sv + z, Sw + (z if iso else 1.0); u, x = self.V - mv, self.W - mw
            g = a * np.exp(-0.5 * (u * u / Tv + x * x / Tw)) / (2 * np.pi * np.sqrt(Tv) * np.sqrt(Tw)); p += g; pw += -x / Tw * g; pww += (x * x / Tw ** 2 - 1 / Tw) * g
        return (p, pw, pww) if second else (p, pw)
    wt = None


def check_2d(S, iso=False):
    p1, q1, qq1 = S.f(1.0, iso, True); p1 = p1.real; ok1 = p1 > 1e-280; m = float(p1[ok1].sum())
    I = lambda x, ok: float(np.sum(np.where(ok, x, 0.0)) / m)
    A = I(np.abs(q1) ** 2 / np.where(ok1, p1, 1.0), ok1)
    if not (0 < A < 1): return None
    L = LN(1 / A); row = {'A': A, 'L': float(L)}
    if L < 4: return row
    P = params(L); T, s = P['T'], P['s']; row.update({k: float(v) for k, v in P.items()})
    pT = S.f(T, iso)[0].real; okT = pT > 1e-280; pTs = np.where(okT, pT, 1.0)
    N = lambda z: I(np.abs(S.f(z, iso)[1]) ** 2 / pTs, okT)
    out = [N(T + P['r3'] * np.exp(1j * a)) for a in ANG]; inn = [N(T + P['r1'] * np.exp(1j * a)) for a in ANG]
    ps, qs = S.f(1 - s, iso); ps = ps.real; qs = qs.real; oks = (ps > 1e-280) & okT; Ns = I(qs ** 2 / pTs, okT); th = P['theta']
    if iso:
        b_out = T ** 3 / (T ** 2 - P['r3'] ** 2) ** 2; b_in = A / (1 - (P['r1'] / P['delta']) ** 2); b3 = (A / P['eta']) ** (1 - th) * b_out ** th
        b3f = 12.62 * (L / LN(2 * L)) * A; c = T / (1 - s); bbad = np.sqrt(6) * T ** 2 / (1 - s) ** 3 * A; bgood = 2 * np.e * ((1 - s) / T) * Ns; bF = 70.2 * (L / LN(2 * L)) * A
    else:
        b_out = 2 / np.sqrt(3); b_in = A / np.sqrt(1 - (P['r1'] / P['delta']) ** 2); b3 = (A / np.sqrt(P['eta'])) ** (1 - th) * b_out ** th
        b3f = 4 * np.e / np.sqrt(3) * np.sqrt(L / LN(2 * L)) * A; c = np.sqrt(T / (1 - s)); bbad = np.sqrt(6 * T / (1 - s)) * A; bgood = 2 * np.e * np.sqrt((1 - s) / T) * Ns; bF = 36.1 * np.sqrt(L / LN(2 * L)) * A
    W = np.where(oks, ps / (c * pTs), 1.0); bad = oks & (W < 1 / (2 * np.e)); good = oks & ~bad; pss = np.where(oks, ps, 1.0)
    Pbad = I(pT, bad); Fbad = I(qs ** 2 / pss, bad); Fgood = I(qs ** 2 / pss, good); F = Fbad + Fgood
    row.update({'outer_margin': float(LN(b_out / max(out))), 'inner_margin': float(LN(b_in / max(inn))), 'outer_sup': float(max(out)), 'inner_sup_over_A': float(max(inn) / A),
                'theta_L': float(th * L), 'three_margin': float(LN(b3 / Ns)), 'three_sampled_margin': float((1 - th) * LN(max(inn)) + th * LN(max(out)) - LN(Ns)),
                'three_final_margin': float(LN(b3f / Ns)), 'N_over_A': float(Ns / A), 'W_max': float(W.max()), 'bad_fraction_points': float(bad.mean()),
                'Pbad_margin': float(LN(2 * A * A / Pbad)) if Pbad > 0 else None, 'bad_margin': float(LN(bbad / Fbad)) if Fbad > 0 else None,
                'good_margin': float(LN(bgood / Fgood)), 'F_margin': float(LN(bF / F)), 'F_over_A': float(F / A)})
    if iso:
        Q = I(np.abs(qq1) ** 2 / np.where(ok1, p1, 1.0), ok1); row.update({'Q_over_A': float(Q / A), 'Qsm_margin': float(LN(F / s / Q)), 'Q_margin': float(LN(562 * A * L ** 3 / LN(2 * L) ** 2 / Q))})
    return row


def states_2d():
    import ed_torus as To, ed_vmix as V
    L = lambda n: json.load(open(os.path.join(C.RES, n), encoding='utf-8'))
    tor = L('ed_torus.json'); ts = [(np.array([complex(a, b) for a, b in r['rho_coeffs']]), r['kv'], r['kw']) for r in tor['optimised']] \
        + [(np.array([1.0, 0.5], complex), r['kv'], r['kw']) for r in tor['single']] + [(To.comb_coeffs(r['kap'] ** 2), r['kv'], r['kw']) for r in tor['comb']]
    import ef_q2 as Q
    yield 'torus', [((lambda t=t: Torus(*t)), (lambda x, t=t: Q.torus_sums(t[0], t[1], t[2], x))) for t in ts]
    lat = L('ed_lattice.json'); yield 'lattice', [((lambda r=r: Lattice(r)), (lambda x, r=r: Q.lattice_sums(r, x))) for r in lat['rows'] + lat['optimised']]
    vm = L('ed_vmix.json'); vs = [V.valley(r['R'], r['rs'], r['S'], r['kv']) for r in vm['valley']] + [V.fringe(r['d'], r['kv'], r['E']) for r in vm['fringe'][::2]] + [(r['terms'], r['kw']) for r in vm['random']]
    yield 'vmix', [((lambda v=v: Vmix(v[0], v[1])), (lambda x, v=v: Q.vmix_sums(v[0], v[1], x))) for v in vs]
    at = L('ed_atoms.json'); yield 'atoms', [((lambda r=r: Atoms([tuple(x) for x in r['blobs']])), (lambda x, r=r: Q.atoms_sums([tuple(b) for b in r['blobs']], x))) for r in at['best_excess'] + at['best_D_ratio']]
    yield 'tilt', [(None, (lambda x, r=r: Q.tilt_sums(r, x))) for r in L('ed_scan_tilt.json')['rows']]


def part_2d(iso=False):
    t0 = time.time(); fams = {}
    def end(get):                                                     # the end results from the sums of ef_q2 (real time only; every recorded 2-D state, tilt included)
        s1 = get(1.0)
        if s1 is None or not (0 < s1['F'] < 1): return None
        A = s1['F']; L = LN(1 / A); row = {'A': A, 'L': float(L)}
        if L < 4: return row
        P = params(L); o = get(1 - P['s'])
        if o is None: return None
        row.update({'F_margin': float(LN(36.1 * np.sqrt(L / LN(2 * L)) * A / o['F'])), 'F_over_A': float(o['F'] / A), 'Qsm_margin': float(LN(o['F'] / P['s'] / s1['Q2'])),
                    'Q_margin': float(LN(289 * A * L ** 2.5 / LN(2 * L) ** 1.5 / s1['Q2'])), 'Q_over_A': float(s1['Q2'] / A)}); return row
    for name, items in states_2d():
        rows = []; ends = []
        for mk, get in items:
            if mk is not None:
                try: r = check_2d(mk(), iso)
                except (ValueError, FloatingPointError, MemoryError): r = None
                if r: rows.append(r)
            if not iso:
                try: e = end(get)
                except (ValueError, FloatingPointError, MemoryError): e = None
                if e: ends.append(e)
        if rows: fams[name] = rows
        if ends: fams['end: ' + name] = ends
        print('%-10s %5d states (circles), %5d (end results)  %.0fs' % (name, len(rows), len(ends), time.time() - t0), flush=True)
    return finish('iso' if iso else '2d', fams, t0)


KEYS = ('outer_margin', 'inner_margin', 'three_margin', 'three_final_margin', 'Pbad_margin', 'bad_margin', 'good_margin', 'F_margin', 'Qsm_margin', 'Q_margin')


def finish(part, fams, t0):
    S = {}
    for f, rows in fams.items():
        big = [r for r in rows if r['L'] >= 4]; v = {'n': len(rows), 'n_L_ge_4': len(big), 'L_max': float(max([r['L'] for r in rows] + [0.0]))}
        for k in KEYS:
            x = [r[k] for r in big if r.get(k) is not None]; v[k + '_min'] = float(min(x)) if x else None; v[k + '_n'] = len(x)
        for k in ('theta_L', 'N_over_A', 'F_over_A', 'Q_over_A', 'inner_sup_over_A', 'outer_sup', 'rho_star', 'JS_ratio', 'W_max'):
            x = [r[k] for r in (rows if k in ('rho_star', 'JS_ratio') else big) if r.get(k) is not None]; v[k + '_max'] = float(max(x)) if x else None
        x = [r['three_sampled_margin'] for r in big if r.get('three_sampled_margin') is not None]; v['three_sampled_margin_min'] = float(min(x)) if x else None
        v['n_with_bad_set'] = sum(1 for r in big if r.get('Pbad_margin') is not None); S[f] = v
    out = {'part': part, 'families': S, 'rows': fams, 'n_states': sum(len(v) for v in fams.values()), 'n_L_ge_4': sum(v['n_L_ge_4'] for v in S.values()), 'env': C.env()}
    tol = 2e-3; al = {}
    for k in KEYS:
        x = [v[k + '_min'] for v in S.values() if v[k + '_min'] is not None]; al[k] = float(min(x)) if x else None
    out['all_min_margins'] = al; out['holds'] = {k: (v is None or v > -tol) for k, v in al.items()}; out['theta_L_max'] = max(v['theta_L_max'] for v in S.values() if v['theta_L_max'] is not None)
    out['ok'] = bool(all(out['holds'].values()) and out['theta_L_max'] < 1); out['seconds'] = time.time() - t0
    C.save('eg_circ_%s.json' % part, out, C.RES_EG)
    for f, v in S.items():
        print('%-24s n %5d (L>=4: %5d, L max %.0f): ' % (f, v['n'], v['n_L_ge_4'], v['L_max']) + '  '.join('%s %s' % (k[:-7], ('%+.3f' % v[k + '_min']) if v[k + '_min'] is not None else '-') for k in KEYS)
              + '  | N/A %s F/A %s rho* %s JS %s' % tuple(('%.3g' % v[k]) if v[k] is not None else '-' for k in ('N_over_A_max', 'F_over_A_max', 'rho_star_max', 'JS_ratio_max')))
    print('all:', {k: (round(v, 4) if v is not None else None) for k, v in al.items()}, ' theta L max %.4f' % out['theta_L_max'], ' ok', out['ok'])
    return out


if __name__ == '__main__':
    np.seterr(all='ignore')
    part = sys.argv[1]; fn = os.path.join(C.RES_EG, 'eg_circ_%s.json' % part)
    if os.path.exists(fn) and os.environ.get('EG_FORCE') != '1': print('exists, skipped'); sys.exit(0)
    {'1d': part_1d, '1ds': lambda: part_1d(10), '2d': part_2d, 'iso': lambda: part_2d(True)}[part]()     # '1ds': every 10th state of each family (the light rerun of the verification)
