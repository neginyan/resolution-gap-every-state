"""Stage EG: fresh re-evaluation (no cache) of extreme rows of eg_direct.json: for sigma = 1e-3, 1e-4, 1e-6 the two-blob state of family W with the largest
|eps| in stage EF's record (the artefact of the cut v-grid), the state inside O with the smallest g sqrt(l), and the state with the largest |residual|/bound.
Also the lower bound of Theorem EG at several sigma (closed forms).  Output results_EG/eg_recheck.json.  Python 3.11 grammar."""
import os, sys, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
os.environ['EG_NOCACHE'] = '1'
import numpy as np
import ed_common as C
import ee_engine as E
import ef_wide as W
from eg_direct import SET


def psi(L): return 289 * L ** 2.5 / np.log(2 * L) ** 1.5 + 4 * L + 3.57


def theorem(sig, K=400.0, c=np.sqrt(0.1)):
    l = np.log(1 / sig); al = 1 / (2 * (np.pi ** 2 + sig ** 2)); Lb = 2 * l + np.log(1 / al); kb = np.sqrt(psi(Lb))
    bp = 2 * sig / (1 - 2 * sig) + 1.291 * np.sqrt(2 / ((1 - 2 * sig) * K)) + 8.4258 * np.sqrt(2 * c / (K * np.sqrt(l)))
    t = [1.0, c / (2 * np.sqrt(l)), 0.14 * (1 - 2 * sig), (1 - 2 * sig) * (1 - bp) / (4 * np.sqrt(2) * kb), 0.75 * al / (al * sig + 1268 * kb)]
    return {'sigma': sig, 'l': l, 'L_bar': Lb, 'k_bar': kb, 'beta_prime': bp, 'terms': t, 'bound': min(t), 'J_in_O_max': 112 * sig ** 2 + 1268 * sig * kb,
            'closed_form': 1.7e-6 * np.log(2 * Lb) ** 0.75 / Lb ** 1.25, 'last_is_min': bool(t[4] == min(t))}


def main():
    np.seterr(all='ignore'); t0 = time.time(); d = json.load(open(os.path.join(C.RES_EG, 'eg_direct.json'), encoding='utf-8')); rows = d['rows']; pick = []
    for sig in (1e-3, 1e-4, 1e-6):
        tb = [r for r in rows if r['two_blobs'] and r['sigma'] == sig]
        if tb: pick.append(max(tb, key=lambda r: abs(r['eps_EF_record'])))
    pick.append(min([r for r in rows if r['in_O']], key=lambda r: r['g_sqrt_l'])); pick.append(max(rows, key=lambda r: r['R_corr_over_bound'])); out = []
    for r in pick:
        bg, rip, Lv = {'W': W.build, 'N': W.build_N, 'M': W.build_M}[r['family']](np.array(r['z']), r['sigma'])
        q = E.finish(E.ev(r['sigma'], bg, rip, Lv, dict(SET)[r['setting']])); dd = q['direct']                    # the setting the row was computed with
        out.append({'record': r['record'], 'family': r['family'], 'sigma': r['sigma'], 'eps_EF_record': r['eps_EF_record'], 'eps': q['eps'], 'g': q['g'], 'residual': dd['residual'], 'bound': dd['bound'],
                    'equal': bool(q['eps'] == r['eps'] and q['g'] == r['g'] and dd['residual'] == r['residual'])})
        print('%-28s %s sigma %g: eps %.3e (EF record %.3e), g %.5g, residual/bound %.3f, equal to the record %s' % (r['record'][8:-5], r['family'], r['sigma'], q['eps'], r['eps_EF_record'], q['g'], abs(dd['residual']) / dd['bound'], out[-1]['equal']), flush=True)
    th = [theorem(s) for s in (1e-15, 1e-20, 1e-30, 1e-50, 1e-100)]
    for t in th: print('Theorem EG, sigma %g: k_bar %.4g, J inside O <= %.3e, bound %.3e (closed form %.3e), last term is the minimum: %s' % (t['sigma'], t['k_bar'], t['J_in_O_max'], t['bound'], t['closed_form'], t['last_is_min']))
    res = {'fresh': out, 'all_equal': bool(all(x['equal'] for x in out)), 'theorem': th, 'theorem_ok': bool(all(t['J_in_O_max'] <= 2e-9 and t['last_is_min'] and t['bound'] >= t['closed_form'] and t['beta_prime'] <= 0.284 for t in th)), 'env': C.env()}
    res['ok'] = bool(res['all_equal'] and res['theorem_ok']); res['seconds'] = time.time() - t0; C.save('eg_recheck.json', res, C.RES_EG)


if __name__ == '__main__':
    main()
