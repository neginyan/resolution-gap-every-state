"""Stage EH: the engine states of the joint checks, in one record.
  (i)  stage EE's 445 recorded states (ee_checks.rows_all: accuracy filter, theta without rounding).  Their v-grids are not affected by the error
       found in stage EG (every background with a v-centre in stage EE was evaluated by the periodic branch, which had one grid per v-centre; checked here).
  (ii) stage EG's 247 states of eg_direct.json (stage EF's wide-state records re-evaluated with the corrected v-grid), re-evaluated at the same setting
       (from the cache of stage EG when EH_CACHE / EG_CACHE points to it; afresh otherwise).
Output results_EH/eh_states.json (only the fields the joint checks use).  Python 3.11 grammar."""
import os, sys, json, glob, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import ed_common as C
import ee_engine as E
import ee_checks as K
import ef_wide as W

KEYS = ('sigma', 'g', 'g_sqrt_l', 'eps', 'theta', 'k_eff2', 'C_w2', 'branch_a', 'branch_ap', 'H1', 'H2', 'A_over_2sigma', 'theta_over_case1')


def compact(q, src, fam):
    r = {k: (float(q[k]) if not isinstance(q[k], bool) else q[k]) for k in KEYS}
    r['T1'] = {k: float(q['T1'][k]) for k in ('A', 'B', 'D', 'K', 'S2', 'SV', 'k_eff2')}
    if 'direct' in q: r['direct'] = {k: float(v) for k, v in q['direct'].items()}
    r['source'] = src; r['family'] = fam; return r


def main():
    t0 = time.time(); rows = []; vc_nonperiodic = 0
    for f in sorted(glob.glob(os.path.join(C.RES_EE, 'ee_fc_*.json'))):
        for r in json.load(open(f, encoding='utf-8'))['rows']:
            st = r['state']
            if any(len(b) > 3 and abs(b[3]) > 0 for b in st['bg']) and not st['kw'].get('periodic_kb'): vc_nonperiodic += 1
    for r in K.rows_all(): rows.append(compact(r, 'EE ' + r['source'], r['fam']))
    n_ee = len(rows)
    d = json.load(open(os.path.join(C.RES_EG, 'eg_direct.json'), encoding='utf-8'))
    for x in d['rows']:
        bg, rip, Lv = {'W': W.build, 'N': W.build_N, 'M': W.build_M}[x['family']](np.array(x['z']), x['sigma'])
        q = E.finish(E.ev(x['sigma'], bg, rip, Lv, dict(nper=12, coarse=100, per_res=10)))
        assert q['g'] == x['g'], 'the re-evaluation differs from eg_direct.json'
        rows.append(compact(q, 'EG ' + x['record'], x['family']))
    out = {'rows': rows, 'n': len(rows), 'n_EE': n_ee, 'n_EG': len(rows) - n_ee, 'EE_v_centred_nonperiodic': vc_nonperiodic, 'env': C.env(), 'seconds': time.time() - t0}
    out['ok'] = bool(vc_nonperiodic == 0 and n_ee == 445 and out['n_EG'] == 247); C.save('eh_states.json', out, C.RES_EH)
    print('%d states (EE %d, EG %d); EE states with a v-centred background outside the periodic branch: %d' % (out['n'], n_ee, out['n_EG'], vc_nonperiodic))


if __name__ == '__main__':
    main()
