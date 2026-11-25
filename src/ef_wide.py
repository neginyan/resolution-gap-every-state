"""Stage EF, 1.1 (2): the wide states (A = sigma^2 J_ww below C* sigma^2: no structure finer than O(0.1) in w) -- the search for a large remainder eps.
Lemma EF-2 (report): if A >= C*(eb) sigma^2, C*(eb) = (2/eb)(7.77 + 6.18/eb), then |B + R3| <= eb (2 A_w + N_phi) for every state; C*(1/2) = 80.5.
So the open states are those with A < 80.5 sigma^2.  Here: states of paper F's engine (exact point fields) built from wide Gaussians,
  bg  : one or two Gaussians in (v, w) of raw widths V, W in [0.03, 3], centres anywhere within +-3, weights
  rip : a smooth modulation 1 + c cos(kv v + kw w + ph) of raw wavelength O(0.3 .. 10) (optional)
at sigma = 1e-3, 1e-4, 1e-6.  Searches (Latin hypercube + Nelder-Mead; every evaluation stored at once in results_EF/cache):
  target 'eps'  : the largest eps (the sign that lowers g)        target 'abs' : the largest |eps|        target 'gsl' : the smallest g sqrt(l)
each with the penalty 10 max(0, log10(A/(80.5 sigma^2))) (stay among the open states).
Power: (P1) the strengthened claim "|eps| <= 1e-9 on wide states" must break; (P2) a bump planted in the objective (height 5 in log10 units,
width 0.15 and 0.3 of the box) must be found.
Family N (targets 'absN', 'epsN', 'gslN'): BOX_N below, penalty 10 max(0, log10((A + theta B/2)/(112 sigma^2))).
Family M (target 'absM'): BOX_M below (a patch of v-ripple at the straight point), the penalty of family N and 3 max(0, log10(30 sigma/B))
(inside O with B >= 30 sigma: the numerical test of Lemma EF-O).
usage: python ef_wide.py SIGMA TARGET [plant WIDTH] -> results_EF/ef_wide_<sigma>_<target>[_plant_w..].json.  Python 3.11 grammar."""
import os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.optimize import minimize
import ed_common as C
import ee_engine as E

CSTAR = 80.5
BOX = np.array([(np.log(0.03), np.log(3.0)),      # W1  raw w-width of blob 1
                (np.log(0.03), np.log(3.0)),      # V   raw v-width (both blobs)
                (-3.0, 3.0),                      # wc1
                (0.0, 0.5),                       # f2  weight of blob 2
                (np.log(0.03), np.log(3.0)),      # W2
                (-3.0, 3.0),                      # wc2
                (-3.0, 3.0),                      # vc2 (raw)
                (0.0, 0.95),                      # c   modulation amplitude
                (np.log(0.6), np.log(20.0)),      # kv  raw wavenumber
                (np.log(0.6), np.log(20.0))])     # kw


def build(z, sig):
    s = np.sqrt(sig); W1, V, wc1, f2, W2, wc2, vc2, c, kv, kw = z; W1, V, W2, kv, kw = np.exp([W1, V, W2, kv, kw])
    bg = [(1 - f2, wc1 / s, W1 / s)] + ([(f2, wc2 / s, W2 / s, vc2 / sig)] if f2 > 0.01 else [])
    rip = [dict(c=float(c), ell=30.0 / s, kv=float(kv * sig), kw=float(kw * sig))] if c > 0.01 else []
    return bg, rip, V / sig


# family N ("narrow window on a wide background": the states of the open set O with v-information of scale sigma .. 1, concentrated in a w-window)
BOX_N = np.array([(np.log(0.1), np.log(3.0)),     # W1  raw w-width of the background
                  (-3.0, 3.0),                    # wc1 raw centre of the background
                  (0.05, 0.95),                   # c   ripple amplitude
                  (0.0, 1.0),                     # x   ln kv_raw = ln 0.6 + x (ln(1/sigma) - ln 0.6)   (v-wavelength from 10 down to 2 pi sigma)
                  (np.log(0.3), np.log(20.0)),    # n   number of v-periods across the state: V_raw = 2 pi n/kv_raw
                  (0.0, 1.0),                     # y   ln ell_raw = ln(0.3 sqrt(sigma)) + y (ln 3 - ln(0.3 sqrt(sigma)))   (w-window of the ripple)
                  (-3.0, 3.0),                    # w0  raw centre of the window
                  (0.0, 1.0)])                    # u   kw_raw = 0 if u < 0.3, else 0.6 .. 20 (log)


def build_N(z, sig):
    s = np.sqrt(sig); lW1, wc1, c, x, ln_n, y, w0, u = z
    kv = np.exp(np.log(0.6) + x * (np.log(1 / sig) - np.log(0.6))); V = 2 * np.pi * np.exp(ln_n) / kv
    ell = np.exp(np.log(0.3 * s) + y * (np.log(3.0) - np.log(0.3 * s))); kw = 0.0 if u < 0.3 else float(np.exp(np.log(0.6) + (u - 0.3) / 0.7 * np.log(20 / 0.6)))
    return [(1.0, wc1 / s, np.exp(lW1) / s)], [dict(c=float(c), ell=float(ell / s), w0=float(w0 / s), kv=float(kv * sig), kw=float(kw * sig))], V / sig


# family M ("a patch at the straight point": the states of O with the most v-information; Lemma EF-O says B <= sigma max(280, 1268 k_eff) there)
BOX_M = np.array([(np.log(0.1), np.log(1.0)),     # W1  raw w-width of the background (centred at w = 0: the straight point z = 0)
                  (np.log(20.0), np.log(300.0)),  # V   v-width of the background in units of sigma
                  (0.05, 0.95),                   # c   ripple amplitude
                  (np.log(0.3), np.log(2.0)),     # kv  wavenumber of the ripple times sigma
                  (np.log(3.0), np.log(40.0)),    # Rv  v-window of the ripple in units of sigma
                  (np.log(0.3), np.log(5.0))])    # ell w-window of the ripple in units of sqrt(sigma)


def build_M(z, sig):
    W1, V, kv, Rv, ell = np.exp([z[0], z[1], z[3], z[4], z[5]]); s = np.sqrt(sig)
    return [(1.0, 0.0, W1 / s)], [dict(c=float(z[2]), ell=float(ell), kv=float(kv), Rv=float(Rv))], float(V)


def evaluate(z, sig, fam='W'):
    import ee_checks as K
    bg, rip, Lv = build(z, sig) if fam == 'W' else (build_N(z, sig) if fam == 'N' else build_M(z, sig))
    q = E.finish(E.ev(sig, bg, rip, Lv, dict(nper=12, coarse=100, per_res=10)))
    q['A_over_s2'] = q['T1']['A'] / sig ** 2; q['B_over_s2'] = q['T1']['B'] / sig ** 2; q['residual'] = K.branch_residual(q)
    q['num_over_s2'] = abs(q['eps']) * 2 * np.sqrt(max(q['branch_a'] * q['branch_ap'], 0.0)) * sig * (q['T1']['A'] + q['T1']['B']) / sig ** 2   # |B + R3|/sigma^2
    return q


CO = 112.0     # the open set O = {A + theta B/2 < 112 sigma^2} (Lemma EF-2, report)


def run(sig, target, plant=None, n_start=16, maxfev=30, seed=0):
    fam = 'N' if 'N' in target[3:] else ('M' if 'M' in target[3:] else 'W'); tg = target[:3]; big_B = fam == 'M'      # family M: inside O and B >= 30 sigma
    rng = np.random.default_rng(seed + {'eps': 1, 'abs': 2, 'gsl': 3}[tg] + (100 if plant else 0) + int(-np.log10(sig)) + {'W': 0, 'N': 1000, 'M': 5000}[fam])
    bx = {'W': BOX, 'N': BOX_N, 'M': BOX_M}[fam]; lo, hi = bx[:, 0], bx[:, 1]; d = len(lo); log = []; t0 = time.time()
    center = lo + (hi - lo) * rng.uniform(0.25, 0.75, d) if plant else None

    def bump(z):
        u = (np.asarray(z) - center) / ((hi - lo) * plant); return 5.0 * float(np.exp(-0.5 * np.sum(u * u)))

    def f(z):
        z = np.clip(z, lo, hi)
        try: q = evaluate(z, sig, fam)
        except (ValueError, FloatingPointError, OverflowError, MemoryError, np.linalg.LinAlgError, ZeroDivisionError, KeyError) as e:
            log.append({'z': list(map(float, z)), 'error': repr(e)[:80]}); return 1e6
        if not (np.isfinite(q['g']) and q['g'] > 0 and q['residual'] <= 0.005):
            log.append({'z': list(map(float, z)), 'error': 'accuracy filter'}); return 1e6
        obj = {'eps': -np.log10(max(q['eps'], 1e-300)), 'abs': -np.log10(max(abs(q['eps']), 1e-300)), 'gsl': np.log10(q['g_sqrt_l'])}[tg]
        Mo = q['A_over_s2'] + max(q['theta'], 0.0) * q['B_over_s2'] / 2                      # (A + theta B/2)/sigma^2
        pen = 10 * max(0.0, np.log10(q['A_over_s2'] / CSTAR)) if fam == 'W' else 10 * max(0.0, np.log10(Mo / CO))
        if big_B: pen += 3 * max(0.0, np.log10(30.0 / max(q['B_over_s2'] * sig, 1e-300)))
        val = obj + pen - (bump(z) if plant else 0.0)
        e = {k: q[k] for k in ('eps', 'g', 'g_sqrt_l', 'theta', 'A_over_s2', 'B_over_s2', 'residual', 'num_over_s2', 'k_eff2', 'C_w2', 'H1', 'H2', 'W2')}
        e.update({'z': list(map(float, z)), 'objective': float(val), 'open': bool(q['A_over_s2'] < CSTAR if fam == 'W' else Mo < CO), 'in_O': bool(Mo < CO), 'M_over_s2': float(Mo),
                  'B_over_sigma': float(q['B_over_s2'] * sig)}); log.append(e)
        return val if np.isfinite(val) else 1e6

    Z = lo + (hi - lo) * (rng.permuted(np.tile(np.arange(n_start), (d, 1)), axis=1).T + rng.random((n_start, d))) / n_start
    vals = [f(z) for z in Z]
    for i in np.argsort(vals)[:3]:
        z0 = Z[i]
        minimize(f, z0, method='Nelder-Mead', options={'maxfev': maxfev, 'xatol': 1e-3, 'fatol': 1e-4,
                                                       'initial_simplex': np.vstack([z0] + [z0 + 0.12 * (hi - lo) * np.eye(d)[j] for j in range(d)])})
    good = [e for e in log if 'objective' in e]; op = [e for e in good if e['open']]
    out = {'sigma': sig, 'target': target, 'family': fam, 'plant': plant, 'n_eval': len(log), 'n_good': len(good), 'n_open': len(op), 'n_in_O': sum(e.get('in_O', False) for e in good),
           'evaluations': log, 'env': C.env(), 'max_abs_eps_in_O': max([abs(e['eps']) for e in good if e.get('in_O')] + [0.0]),
           'max_B_over_sigma_open': max([e['B_over_sigma'] for e in op] + [0.0]),
           'EF_O_min_margin': min([max(280.0, 1268 * np.sqrt(max(e['k_eff2'], 0.0))) / e['B_over_sigma'] - 1 for e in good if e.get('in_O')] + [np.inf]),
           'max_eps_open': max([e['eps'] for e in op] + [-np.inf]), 'max_abs_eps_open': max([abs(e['eps']) for e in op] + [0.0]),
           'max_abs_eps_all': max([abs(e['eps']) for e in good] + [0.0]), 'min_g_sqrt_l_open': min([e['g_sqrt_l'] for e in op] + [np.inf]),
           'max_num_over_s2_open': max([e['num_over_s2'] for e in op] + [0.0]), 'min_A_over_s2': min([e['A_over_s2'] for e in good] + [np.inf])}
    if op:
        b = max(op, key=lambda e: abs(e['eps'])); out['argmax_abs_eps'] = {k: b[k] for k in b if k != 'objective'}
    if plant:
        b = min(good, key=lambda e: e['objective']); out['plant_center'] = list(map(float, center)); out['plant_best_bump'] = bump(b['z'])
        out['plant_found'] = bool(out['plant_best_bump'] > 2.5); out['plant_width'] = plant
    out['ok'] = True; out['seconds'] = time.time() - t0
    C.save(fname(sig, target, plant), out, C.RES_EF)
    print('wide sigma %g target %s%s: %d evals (%d open: W: A < %.1f sigma^2, N: inside O); open: max eps %.3e, max |eps| %.3e, min g sqrt(l) %.4g, max |B+R3|/sigma^2 %.3e%s  %.0fs' % (
        sig, target, (' PLANT w%g' % plant) if plant else '', len(log), len(op), CSTAR, out['max_eps_open'], out['max_abs_eps_open'], out['min_g_sqrt_l_open'],
        out['max_num_over_s2_open'], (' found %s (bump %.2f of 5)' % (out['plant_found'], out['plant_best_bump'])) if plant else '', out['seconds']), flush=True)


def fname(sig, target, plant):
    return 'ef_wide_%g_%s%s.json' % (sig, target, ('_plant_w%g' % plant) if plant else '')


if __name__ == '__main__':
    np.seterr(all='ignore')
    a = sys.argv; sig = float(a[1]); target = a[2]; plant = float(a[4]) if len(a) > 4 and a[3] == 'plant' else None
    fn = os.path.join(C.RES_EF, fname(sig, target, plant))
    if os.path.exists(fn) and os.environ.get('EF_FORCE') != '1' and os.environ.get('EF_NOCACHE') != '1':
        print('exists, skipped:', os.path.basename(fn)); sys.exit(0)
    run(sig, target, plant)
