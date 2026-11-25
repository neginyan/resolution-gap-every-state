"""Stage EH: fresh re-evaluation (no cache) of a few engine states of eh_states.json: for each of the families W, N, M of stage EG's states the one inside O
with the smallest g sqrt(l), and the one with A + B <= 2e-9 with the smallest g/(3A/(4 sigma (A + B))).  The values must equal the record.
Output results_EH/eh_fresh.json.  Python 3.11 grammar."""
import os, sys, json, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
os.environ['EH_NOCACHE'] = '1'
import numpy as np
import ed_common as C
import ee_engine as E
import ef_wide as W


def main():
    np.seterr(all='ignore'); d = json.load(open(os.path.join(C.RES_EG, 'eg_direct.json'), encoding='utf-8'))['rows']
    st = {(r['source'][3:], r['g']): r for r in json.load(open(os.path.join(C.RES_EH, 'eh_states.json'), encoding='utf-8'))['rows'] if r['source'].startswith('EG')}
    pick = []
    for fam in ('W', 'N', 'M'):
        o = [x for x in d if x['family'] == fam and x['in_O']]
        if o: pick.append(min(o, key=lambda x: x['g_sqrt_l']))
    sm = [x for x in d if x['J'] <= 2e-9]
    if sm: pick.append(min(sm, key=lambda x: x['g'] * x['sigma'] * x['J'] / (0.75 * x['A_over_s2'] * x['sigma'] ** 2)))
    out = []
    for x in pick:
        bg, rip, Lv = {'W': W.build, 'N': W.build_N, 'M': W.build_M}[x['family']](np.array(x['z']), x['sigma'])
        q = E.finish(E.ev(x['sigma'], bg, rip, Lv, dict(nper=12, coarse=100, per_res=10))); r = st[(x['record'], q['g'])] if (x['record'], q['g']) in st else None
        eq = bool(r is not None and q['eps'] == r['eps'] and q['T1']['A'] == r['T1']['A'] and q['T1']['K'] == r['T1']['K'] and q['k_eff2'] == r['k_eff2'])
        out.append({'record': x['record'], 'family': x['family'], 'sigma': x['sigma'], 'g': q['g'], 'equal_to_record': eq})
        print('%-28s %s sigma %g: g %.6g, equal to eh_states.json: %s' % (x['record'], x['family'], x['sigma'], q['g'], eq), flush=True)
    res = {'rows': out, 'all_equal': bool(out and all(r['equal_to_record'] for r in out)), 'env': C.env()}; res['ok'] = res['all_equal']
    C.save('eh_fresh.json', res, C.RES_EH)


if __name__ == '__main__':
    main()
