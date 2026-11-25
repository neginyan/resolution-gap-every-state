"""Stage EE, 1.4: the one-variable core, searched harder than in stage ED (numerics only; a conjecture).
X real, Y = X + Z (Z ~ N(0, 1)), m: R -> [-1, 1], h(y) = E[m(X) | Y = y], a = E[h(Y)^2], D1 = E[h'(Y)^2] = E[Cov(m(X), Z | Y)^2].
Conjecture (1-D core):  D1 <= c a ln(1/a) + C a  (proved in stage ED for 0 <= m <= 1; open for m changing sign).
X is a mixture of Gaussian components (weight, centre, variance >= 0), m constant on each component (stage ED's channel1d, recomputed here
with a log-domain density and log-domain integrals, so that a down to e^-3000 is resolved).
Families (EE_探索計画.md §7; valley_coherent added after the first run, as a reference: the sign-coherent valley, where the tail mechanism gives D1/(a ln(1/a)) -> 2):
  valley_ripple : heavy atoms (m = 0) at +-R, between them a light alternating lattice (spacing d, m = (-1)^j) of weight rs exp(-R^2/2) under a
                  Gaussian envelope of width E; R from 1 to sqrt(2 ln(1/a)) (the depth of the valley against the cost of the light part)
  lattice       : an alternating lattice alone (the ripple mechanism)
  valley_pair   : two light atoms of opposite sign (+1, -1) at distance d inside the valley
  nested        : a valley inside a valley (heavy at +-R1, medium at +-R2 < R1, light alternating lattice in the middle)
  valley_cos    : the valley with a continuous profile m(x) = cos(pi x/d) on a light Gaussian (closed form, exact for any d; point lattices
                  are used only with d >= 0.7, since their alternating sums lose all digits when exp(-pi^2/(2 d^2)) < 1e-12)
  random        : J = 3..9 components, Nelder-Mead on D1/(a ln(1/a)) (a < 1e-3 enforced by a penalty) and on (D1 - a ln(1/a))/a
Power checks: (P1) the strengthened form D1 <= (1/2) a ln(1/a) + a must break; (P2) a bump planted in the random search's parameter space
(depth 5, larger than the objective's spread) must be found.
Output results_EE/ee_core1d.json.  Python 3.11 grammar."""
import os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp
import ed_common as C


def channel(comps, n=20001, span=9.0):
    """a = E[h^2], D1 = E[h'^2] on a grid, the density in log form (components of weight down to 1e-300)."""
    comps = [c for c in comps if c[0] > 0 or len(c) > 4]          # a 5th entry is the log of the weight (weights below e^-745)
    lo = min(c[1] - span * np.sqrt(c[2] + 1) for c in comps); hi = max(c[1] + span * np.sqrt(c[2] + 1) for c in comps)
    y = np.linspace(lo, hi, n); dy = y[1] - y[0]
    lw = np.array([c[4] if len(c) > 4 else np.log(c[0]) for c in comps]); lw -= logsumexp(lw)
    mu = np.array([c[1] for c in comps]); T = np.array([c[2] for c in comps]) + 1.0; mm = np.array([c[3] for c in comps])
    L = lw[:, None] - 0.5 * (y[None, :] - mu[:, None]) ** 2 / T[:, None] - 0.5 * np.log(2 * np.pi * T[:, None])      # log of each component's density
    lp = logsumexp(L, axis=0); pi = np.exp(L - lp[None, :])                                                     # posterior weights of the components
    d1 = -(y[None, :] - mu[:, None]) / T[:, None]                                                                # d/dy log of each component
    h = np.sum(pi * mm[:, None], 0); s = np.sum(pi * d1, 0)
    hp = np.sum(pi * mm[:, None] * d1, 0) - h * s                                                                # h' = Cov_post(m, d log comp)
    with np.errstate(divide='ignore'):
        la = float(logsumexp(lp + np.log(h * h)) + np.log(dy)); lD = float(logsumexp(lp + np.log(hp * hp)) + np.log(dy))    # logs: a down to e^-10^4
    return la, lD


def stats(comps, n=20001):
    try: lna, lnD = channel(comps, n)
    except (ValueError, FloatingPointError): return None
    if not (np.isfinite(lna) and lna < 0 and np.isfinite(lnD)): return None
    la = -lna; Da = float(np.exp(lnD - lna))                       # D1/a
    return {'a': float(np.exp(lna)), 'ln_a': float(lna), 'D_over_a': Da, 'ln_1_over_a': float(la), 'ratio': Da / la, 'excess': Da - la,
            'excess2': Da - 2 * la, 'half_margin': 0.5 * la + 1 - Da}


def stats_cos(R, rs, d, E, n=None, ph=0.0):
    """valley_cos (exact, no cancellation): heavy atoms (m = 0, weight 1/2 each) at +-R (none if R = 0) and a light Gaussian N(0, E^2) of weight
    rs exp(-R^2/2) carrying the continuous profile m(x) = cos(pi x/d + ph) (|m| <= 1).  In closed form q_light(y) = vis N(y; 0, T) cos(w y + ph),
    T = E^2 + 1, w = (pi/d) E^2/T, vis = exp(-(pi/d)^2 E^2/(2T)); everything in logs (a down to e^-10^4)."""
    k = np.pi / d; T = E * E + 1; w = k * E * E / T; lvis = -k * k * E * E / (2 * T)
    lo, hi = -(R + 10 * np.sqrt(T)), R + 10 * np.sqrt(T)
    if n is None:
        dy = min(1.0 / (30 * max(R, 1.0)), 2 * np.pi / (60 * w), 0.02); n = int(min(2.5e6, (hi - lo) / dy) + 1)
    y = np.linspace(lo, hi, n); dy = y[1] - y[0]
    le = np.log(rs) - R * R / 2
    Ll = le - 0.5 * y * y / T - 0.5 * np.log(2 * np.pi * T)
    comps = [Ll]; d1 = [-y / T]
    if R > 0:
        for c in (-R, R):
            comps.append(np.log(0.5) - 0.5 * (y - c) ** 2 - 0.5 * np.log(2 * np.pi)); d1.append(-(y - c))
    L = np.stack(comps); lp = logsumexp(L, axis=0); pi = np.exp(L - lp[None, :]); sc = np.sum(pi * np.stack(d1), 0)
    lpl = L[0] - lp                                                                            # log pi_light
    cs = np.cos(w * y + ph); sn = np.sin(w * y + ph)
    with np.errstate(divide='ignore'):
        lh2 = 2 * (lpl + lvis) + np.log(cs * cs)
        lhp2 = 2 * (lpl + lvis) + np.log(((-y / T - sc) * cs - w * sn) ** 2)
    lna = float(logsumexp(lp + lh2) + np.log(dy)); lnD = float(logsumexp(lp + lhp2) + np.log(dy))
    if not (np.isfinite(lna) and lna < 0): return None
    la = -lna; Da = float(np.exp(lnD - lna))
    return {'a': float(np.exp(lna)), 'ln_a': lna, 'D_over_a': Da, 'ln_1_over_a': la, 'ratio': Da / la, 'excess': Da - la, 'excess2': Da - 2 * la,
            'half_margin': 0.5 * la + 1 - Da, 'n_grid': int(n)}


def lattice(d, E, w=1.0, sgn=True, centre=0.0):
    J = int(np.ceil(4 * E / d)) + 1
    return [(w * float(np.exp(-0.5 * (j * d / E) ** 2)), centre + j * d, 0.0, float((-1) ** j) if sgn else 1.0) for j in range(-J, J + 1)]


def valley_ripple(R, rs, d, E):
    lat = lattice(d, E); tot = sum(c[0] for c in lat)
    return [(0.5, -R, 0.0, 0.0), (0.5, R, 0.0, 0.0)] + [(0.0, c[1], 0.0, c[3], np.log(rs) - R * R / 2 + np.log(c[0] / tot)) for c in lat if c[0] > 0]


def valley_pair(R, rs, d):
    light = rs * np.exp(-R * R / 2)
    return [(0.5, -R, 0.0, 0.0), (0.5, R, 0.0, 0.0), (light / 2, -d / 2, 0.0, 1.0), (light / 2, d / 2, 0.0, -1.0)]


def nested(R1, R2, r2, rs, d, E):
    med = r2 * np.exp(-R1 * R1 / 2 + R2 * R2 / 2 - R2 * R2 / 2); light = rs * np.exp(-R1 * R1 / 2 - R2 * R2 / 2)
    lat = lattice(d, E); tot = sum(c[0] for c in lat)
    return ([(0.5, -R1, 0.0, 0.0), (0.5, R1, 0.0, 0.0), (med / 2, -R2, 0.0, 0.0), (med / 2, R2, 0.0, 0.0)]
            + [(light * c[0] / tot, c[1], 0.0, c[3]) for c in lat])


def comps_of(z, J):
    out = []
    for j in range(J):
        la, mu, lS, mm = z[4 * j:4 * j + 4]
        out.append((float(np.exp(np.clip(la, -200, 3))), float(np.clip(mu, -16, 16)), float(np.exp(np.clip(lS, -8, 2))), float(np.tanh(mm))))
    return out


def main():
    t0 = time.time(); out = {'what': __doc__.split('\n')[0], 'env': C.env()}
    fam = {'valley_ripple': [], 'lattice': [], 'valley_pair': [], 'nested': [], 'random': []}
    fam['valley_coherent'] = []
    for R in (2.0, 4.0, 6.0, 8.0, 10.0, 15.0, 20.0, 30.0, 40.0, 60.0, 80.0):       # the sign-coherent valley (m = 1 on one light atom): the tail mechanism
        for rs in (1e-6, 1e-3, 1e-1, 1.0, 10.0, 1e3):
            s = stats([(0.5, -R, 0.0, 0.0), (0.5, R, 0.0, 0.0), (0.0, 0.0, 0.0, 1.0, np.log(rs) - R * R / 2)], n=int(max(20001, 400 * R)))
            if s: s.update({'R': R, 'rs': rs}); fam['valley_coherent'].append(s)
        v = [s for s in fam['valley_coherent'] if s['R'] == R]
        if v:
            b = max(v, key=lambda s: s['ratio'])
            print('valley_coherent R %4.1f: max ratio %.4f (ln(1/a) %.1f, rs %g); predicted R^2/(R^2/2 + ln R + 2 ln(1/rs)) = %.4f  %.0fs' % (
                R, b['ratio'], b['ln_1_over_a'], b['rs'], R * R / (R * R / 2 + np.log(R) + 2 * np.log(1 / b['rs'])), time.time() - t0), flush=True)
    for R in (1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 13.0, 16.0, 20.0, 30.0, 40.0, 60.0):
        for rs in (1e-6, 1e-3, 1e-1, 1.0, 10.0, 1e3):
            for d in (0.7, 0.8, 1.0, 1.5, 2.0):         # point lattices only with d >= 0.7 (finer ones: valley_cos, exact)
                for E in (0.5, 1.5, 4.0):
                    s = stats(valley_ripple(R, rs, d, E), n=int(max(20001, 400 * R)))
                    if s: s.update({'R': R, 'rs': rs, 'd': d, 'E': E}); fam['valley_ripple'].append(s)
        v = [s for s in fam['valley_ripple'] if s['R'] == R and s['ln_a'] < np.log(1e-3)]
        if v:
            b = max(v, key=lambda s: s['ratio'])
            print('valley_ripple R %4.1f: %3d states (a < 1e-3), max ratio %.4f (a %.2e, ln(1/a) %.1f, rs %g, d %g, E %g); max excess %+.3f  %.0fs' % (
                R, len(v), b['ratio'], b['a'], b['ln_1_over_a'], b['rs'], b['d'], b['E'], max(s['excess'] for s in v), time.time() - t0), flush=True)
    for d in (0.7, 0.8, 0.9, 1.0, 1.2, 1.5, 2.0, 3.0):
        for E in (2.0, 12.0, 40.0):
            s = stats(lattice(d, E), n=40001)
            if s: s.update({'d': d, 'E': E}); fam['lattice'].append(s)
    for R in (2.0, 4.0, 6.0, 8.0, 12.0, 16.0):
        for rs in (1e-2, 1.0, 1e2):
            for d in (0.2, 0.5, 1.0, 2.0, 4.0):
                s = stats(valley_pair(R, rs, d))
                if s: s.update({'R': R, 'rs': rs, 'd': d}); fam['valley_pair'].append(s)
    for R1 in (6.0, 10.0, 14.0):
        for R2 in (2.0, 4.0):
            for r2 in (1e-2, 1.0):
                for rs in (1e-2, 1.0, 1e2):
                    for d in (0.7, 1.0):
                        s = stats(nested(R1, R2, r2, rs, d, 1.5))
                        if s: s.update({'R1': R1, 'R2': R2, 'r2': r2, 'rs': rs, 'd': d}); fam['nested'].append(s)
    for k in ('lattice', 'valley_pair', 'nested'):
        v = [s for s in fam[k] if s['ln_a'] < np.log(1e-3)]
        if v: print('%-13s %3d states (a < 1e-3), max ratio %.4f, max excess %+.3f  %.0fs' % (k, len(v), max(s['ratio'] for s in v), max(s['excess'] for s in v), time.time() - t0), flush=True)
    fam['valley_cos'] = []
    for R in (0.0, 2.0, 4.0, 6.0, 10.0, 15.0, 20.0, 30.0, 40.0, 60.0, 80.0):
        for rs in (1e-6, 1e-2, 1.0, 1e2):
            for d in (0.15, 0.25, 0.35, 0.5, 0.7, 1.0, 1.5, 2.0):
                for E in (0.5, 1.5, 4.0, 10.0):
                    if R == 0 and rs != 1.0: continue
                    s = stats_cos(R, rs, d, E)
                    if s: s.update({'R': R, 'rs': rs, 'd': d, 'E': E}); fam['valley_cos'].append(s)
        v = [s for s in fam['valley_cos'] if s['R'] == R and s['ln_a'] < np.log(1e-3)]
        if v:
            b = max(v, key=lambda s: s['ratio'])
            print('valley_cos R %4.1f: %3d states (a < 1e-3), max ratio %.4f (ln(1/a) %.1f, rs %g, d %g, E %g); max excess %+.3f  %.0fs' % (
                R, len(v), b['ratio'], b['ln_1_over_a'], b['rs'], b['d'], b['E'], max(s['excess'] for s in v), time.time() - t0), flush=True)
    # random + Nelder-Mead
    rng = np.random.default_rng(29)
    def run_opt(J, key, n_start, plant=None, width=1.0):
        best = None; log = []
        def f(z):
            s = stats(comps_of(z, J), n=4001)
            if s is None: return 1e3
            val = -s[key] + (0.0 if (key != 'ratio' or s['ln_a'] < np.log(1e-3)) else 10 * (s['ln_a'] - np.log(1e-3)) / np.log(10))
            if plant is not None: val -= 5.0 * np.exp(-0.5 * np.sum(((z - plant) / width) ** 2))
            log.append(val); return val
        for _ in range(n_start):
            z0 = np.column_stack([rng.normal(-6, 8, J), rng.normal(0, 5, J), rng.normal(-3, 2, J), rng.normal(0, 2, J)]).ravel()
            r = minimize(f, z0, method='Nelder-Mead', options={'maxfev': 700, 'xatol': 1e-4, 'fatol': 1e-8})
            if best is None or r.fun < best.fun: best = r
        return best
    for J in (3, 4, 5, 7, 9):
        for key in ('ratio', 'excess'):
            b = run_opt(J, key, 12); cs = comps_of(b.x, J); s = stats(cs, n=40001)
            if s is None: continue
            s.update({'J': J, 'objective': key, 'comps': cs}); fam['random'].append(s)
            print('random J %d, largest %-6s: ratio %.4f excess %+.4f (a %.3e)  %.0fs' % (J, key, s['ratio'], s['excess'], s['a'], time.time() - t0), flush=True)
    # power check P2: a planted bump (depth 5) in the parameter space of J = 3
    J = 3; z_pl = np.column_stack([rng.normal(-6, 8, J), rng.normal(0, 5, J), rng.normal(-3, 2, J), rng.normal(0, 2, J)]).ravel()
    plants = []
    for width in (2.0, 1.0):
        b = run_opt(J, 'excess', 12, plant=z_pl, width=width); dist = float(np.linalg.norm(b.x - z_pl))
        bump = 5.0 * float(np.exp(-0.5 * np.sum(((b.x - z_pl) / width) ** 2)))
        plants.append({'search': 'random mixtures J = 3, Nelder-Mead (12 parameters)', 'width': width, 'distance': dist, 'bump_at_best': bump, 'found': bool(bump > 2.5), 'dimension': int(len(z_pl))})
        print('plant width %.1f: found %s (bump %.2f of 5, distance %.2f)' % (width, plants[-1]['found'], bump, dist), flush=True)
    # structured local search: Nelder-Mead in the parameters (R, log10 rs, d, E) of valley_ripple, maximising D1/(a ln(1/a)) (R <= 30), from 8 starts;
    # and its power check: a bump of height 1 (larger than the spread of the ratio, ~0.4) planted at a random point of the box must be found.
    bx = np.array([(1.0, 80.0), (-6.0, 3.0), (0.15, 3.0), (np.log(0.3), np.log(20.0))])
    def vr(z):
        z = np.clip(z, bx[:, 0], bx[:, 1]); return stats_cos(z[0], 10 ** z[1], z[2], np.exp(z[3]))
    def vr_search(plant=None, width=0.15, n_start=8):
        best = None
        def f(z):
            s = vr(z)
            if s is None: return 1e3
            val = -s['ratio'] + (0.0 if s['ln_a'] < np.log(1e-3) else 10.0)
            if plant is not None: val -= 1.0 * np.exp(-0.5 * np.sum(((np.clip(z, bx[:, 0], bx[:, 1]) - plant) / ((bx[:, 1] - bx[:, 0]) * width)) ** 2))
            return val
        for _ in range(n_start):
            z0 = bx[:, 0] + (bx[:, 1] - bx[:, 0]) * rng.random(4)
            r = minimize(f, z0, method='Nelder-Mead', options={'maxfev': 150, 'xatol': 1e-3, 'fatol': 1e-6,
                                                               'initial_simplex': np.vstack([z0] + [z0 + 0.15 * (bx[:, 1] - bx[:, 0]) * np.eye(4)[j] for j in range(4)])})
            if best is None or r.fun < best.fun: best = r
        return best
    b = vr_search(); zb = np.clip(b.x, bx[:, 0], bx[:, 1]); s = vr(zb)
    if s:
        s.update({'R': float(zb[0]), 'rs': float(10 ** zb[1]), 'd': float(zb[2]), 'E': float(np.exp(zb[3])), 'search': 'Nelder-Mead in (R, rs, d, E)'})
        f4 = stats_cos(s['R'], s['rs'], s['d'], s['E'], n=4 * s['n_grid']); s['ratio_grid_x4'] = f4['ratio'] if f4 else None; fam['valley_cos'].append(s)
        print('valley_cos, Nelder-Mead in (R, rs, d, E): max ratio %.4f (grid x4: %.4f) at R %.1f rs %.2g d %.3f E %.2f (ln(1/a) %.1f)' % (
            s['ratio'], s['ratio_grid_x4'], s['R'], s['rs'], s['d'], s['E'], s['ln_1_over_a']), flush=True)
    for width in (0.15, 0.3):
        pl = bx[:, 0] + (bx[:, 1] - bx[:, 0]) * rng.uniform(0.2, 0.8, 4); b = vr_search(plant=pl, width=width)
        zb = np.clip(b.x, bx[:, 0], bx[:, 1]); dist = float(np.max(np.abs((zb - pl) / (bx[:, 1] - bx[:, 0]))))
        bump = float(np.exp(-0.5 * np.sum(((zb - pl) / ((bx[:, 1] - bx[:, 0]) * width)) ** 2)))
        plants.append({'search': 'valley_cos Nelder-Mead (4 parameters)', 'width': width, 'distance_box_units': dist, 'bump_at_best': bump, 'found': bool(bump > 0.5), 'dimension': 4})
        print('plant (valley_cos search) width %.2f: found %s (bump %.2f of 1, distance %.3f of the box)' % (width, plants[-1]['found'], bump, dist), flush=True)
    out['plant'] = plants
    allr = [s for v in fam.values() for s in v]
    out.update({'families': fam, 'n_states': len(allr), 'ln_a_min': min(s['ln_a'] for s in allr)})
    small = [s for s in allr if s['ln_a'] < np.log(1e-3)]
    out['max_ratio_a_below_1e-3'] = max(s['ratio'] for s in small); out['max_excess'] = max(s['excess'] for s in allr)
    out['n_half_broken'] = int(sum(s['half_margin'] < 0 for s in allr))
    out['max_ratio_a_below_1e-3_sign_changing'] = max(s['ratio'] for f in ('valley_ripple', 'valley_cos', 'lattice', 'valley_pair', 'nested', 'random') for s in fam[f] if s['ln_a'] < np.log(1e-3))
    out['max_excess2'] = max(s['excess2'] for s in allr)                 # (D1 - 2 a ln(1/a))/a
    # trend: the largest ratio in bins of ln(1/a) (a stop signal would be a ratio growing with ln(1/a))
    bins = {}
    for s in small:
        k = int(min(s['ln_1_over_a'], 4000) // (10 if s['ln_1_over_a'] < 300 else 500)) * (10 if s['ln_1_over_a'] < 300 else 500)
        if k not in bins or s['ratio'] > bins[k]['ratio']: bins[k] = {'ratio': s['ratio'], 'excess': s['excess'], 'family': [f for f, v in fam.items() if s in v][0]}
    out['max_ratio_by_ln_1_over_a'] = {str(k): bins[k] for k in sorted(bins)}
    out['per_family'] = {f: {'n': len(v), 'n_a_below_1e-3': sum(s['a'] < 1e-3 for s in v),
                             'max_ratio_a_below_1e-3': max([s['ratio'] for s in v if s['ln_a'] < np.log(1e-3)] + [float('nan')]),
                             'max_excess': max([s['excess'] for s in v] + [float('nan')]), 'max_excess2': max([s['excess2'] for s in v] + [float('nan')]),
                             'n_half_broken': int(sum(s['half_margin'] < 0 for s in v))} for f, v in fam.items()}
    out['ok'] = True; out['seconds'] = time.time() - t0
    C.save('ee_core1d.json', out, C.RES_EE)
    print('states %d: max D1/(a ln(1/a)) for a < 1e-3: %.4f (sign changing: %.4f); max (D1 - a ln(1/a))/a: %+.4f; max (D1 - 2 a ln(1/a))/a: %+.4f; (1/2)-form broken on %d; smallest ln a %.1f' % (
        out['n_states'], out['max_ratio_a_below_1e-3'], out['max_ratio_a_below_1e-3_sign_changing'], out['max_excess'], out['max_excess2'], out['n_half_broken'], out['ln_a_min']))
    for k, v in out['max_ratio_by_ln_1_over_a'].items(): print('  ln(1/a) %4s..: max ratio %.4f (%s)' % (k, v['ratio'], v['family']))


if __name__ == '__main__':
    np.seterr(all='ignore'); main()
