"""Stage EE: the lemma checks and the summaries of 1.1 and 1.2, recomputed from the records (results_EE/ee_fc_*.json, ee_layer_*.json).
Lemma checks (margin = right side - left side, the smallest over the states; a check holds when the margin is >= -1e-9 relative):
  AT2'      g >= F_c'/(2 k_eff) (Theorem AT2 with C_w' = E[s_v^2 cos m_z (1 - sigma^2(H_ww + H_vw))]/J_vv):  1 - AT2'/g
  AT2       g >= F_c/(2 k_eff) (paper F):                                                                    1 - AT2/g
  factor    (1 - Phi/4)(1 - delta_w) = E[s_w^2 (1 - phi/4)]/J_ww >= 1/2:                                      sqrt_factor - 1/sqrt2
  EE-C      C_w' >= 1 - (S2 + 9.717 sqrt(A))/B - 8.426 sqrt(theta/B)  (sigma = 1 units):                      C_w' - bound
  branch    (1 - |eps|)((a) + (a')) <= g <= (1 + |eps|)((a) + (a')) (paper F, exact): its numerical residual, allowed 0.5 % (quadrature)
  CR(H2)    under (H2): ln(A 2 (pi^2 + sigma^2)/sigma^2) >= 0 (paper E's Cramer-Rao floor)
  BD1       E_{pi_v}[d(m_z)^2] <= (pi^2/2) theta (paper F)
Summaries: in case 3 (text, op, both; and with (H1), (H2)) the smallest F_c, F_c', g sqrt(l), AT2/g, AT2'/g, the largest |eps|, by family and sigma;
the region R4 of the conditional theorem (A/(A+B) < 2 sigma, theta < 2 sigma c/sqrt l, B >= K sigma) with the smallest C_w' there.
Every record is read through ee_engine.finish (theta from the branch (a'), without the rounding of paper F's formula; case flags redone).
Output results_EE/ee_checks.json.  Python 3.11 grammar."""
import os, sys, glob, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import ed_common as C

RES = C.RES_EE
K_R4 = 3200.0


def rows_all():
    R = []
    for f in sorted(glob.glob(os.path.join(RES, 'ee_fc_*.json'))):
        d = json.load(open(f, encoding='utf-8'))
        for r in d['rows']: R.append(dict(r, source=os.path.basename(f), fam=r.get('family', d['part'])))
    for f in sorted(glob.glob(os.path.join(RES, 'ee_layer_*.json'))):
        d = json.load(open(f, encoding='utf-8'))
        if d.get('plant'): continue                                 # planted objectives: their states are ordinary, but kept apart
        for r in d['evaluations']:
            if 'g' in r: R.append(dict(r, source=os.path.basename(f), fam='layer ' + d['family']))
    import ee_engine as E
    return [E.finish(r) for r in R if ok_num(r)]                   # theta without rounding, case flags redone (ee_engine.finish)


def ok_num(r):
    """the engine's float range and accuracy: EB's windows of the 'wider (H1)' room (focus 1e20 sigma and wider at sigma = 1e-24) are numerically broken
    (g < 0 on three of them, the branch identity off by 80 % on another); they are excluded and counted apart."""
    if not (np.isfinite(r['g']) and r['g'] > 0 and 'wider (H1)' not in str(r.get('label', ''))): return False
    # the branch identity (exact) as the accuracy test of one evaluation: residual above 0.5 % = the quadrature or the float range failed
    # (this removes the delta_w family at sigma = 1e-24: a fine w-ripple on a far blob, residuals up to 800 %)
    if far_unresolved(r): return False
    return bool(branch_residual(r) <= 0.005)


def far_unresolved(r):
    """a background or window far from the origin (two straight points, the delta_w family, windows moved along p) at a sigma so small that
    its position is known in double precision only to more than 1e-4 of the narrowest window width: |w_c| eps_mach / ell > 1e-4
    (w_c, ell in sqrt(sigma) units); likewise a v-centre against the period of the ripple.  This removes those families at sigma <= 1e-24."""
    st = r.get('state')
    if not st: return False
    ells = [x['ell'] for x in st['rip']] or [1.0]; wc = max(abs(b[1]) for b in st['bg']); vc = max((abs(b[3]) if len(b) > 3 else 0.0) for b in st['bg'])
    kv = max([abs(x['kv']) for x in st['rip']] + [1e-9])
    return bool(wc * 2.2e-16 / min(ells) > 1e-4 or vc * 2.2e-16 * kv / (2 * np.pi) > 1e-4)


def branch_residual(r):
    """|g - (a) - (a') + eps 2 sqrt((a)(a'))|/g: zero exactly (paper F's branch identity, g = (a) + (a') - eps (a/den) sqrt(P_w) and
    (a/den) sqrt(P_w) = 2 sqrt((a)(a'))); numerically it is the residual of the integration by parts behind paper E's decomposition."""
    a, ap = r['branch_a'], r['branch_ap']
    return float(abs(r['g'] - a - ap + r['eps'] * 2 * np.sqrt(max(a * ap, 0.0))) / r['g'])


def broken_count():
    n = 0
    for f in sorted(glob.glob(os.path.join(RES, 'ee_fc_*.json'))):
        n += sum(1 for r in json.load(open(f, encoding='utf-8'))['rows'] if not ok_num(r))
    return n


def main():
    t0 = time.time(); R = rows_all(); out = {'n_states': len(R), 'n_numerically_broken': broken_count(), 'lemmas': {}, 'env': C.env()}
    def stat(name, vals, rows, rel=1e-9):
        vals = np.array(vals, float); i = int(np.nanargmin(vals))
        out['lemmas'][name] = {'n': int(np.isfinite(vals).sum()), 'min_margin': float(vals[i]), 'family': rows[i]['fam'], 'sigma': rows[i]['sigma'],
                               'holds': bool(vals[i] >= -rel)}
    stat("AT2' : 1 - AT2'/g", [1 - r['AT2b_over_g'] for r in R], R)
    stat('AT2 (paper F): 1 - AT2/g', [1 - r['AT2_over_g'] for r in R], R)
    stat('factor (1 - Phi/4)(1 - delta_w) >= 1/2: sqrt - 1/sqrt2', [r['sqrt_factor'] - 1 / np.sqrt(2) for r in R], R)
    fin = [r for r in R if np.isfinite(r.get('Cw2_lower', np.nan))]
    stat("EE-C: C_w' - lower bound", [r['C_w2'] - r['Cw2_lower'] for r in fin], fin)
    # the branch identity is exact; its numerical residual measures the quadrature (paper F's one-period averages along the v-envelope)
    stat('branch identity residual <= 0.005 (the accuracy filter itself; margin 0.005 - residual)', [0.005 - branch_residual(r) for r in R], R)
    stat('branch lower bound (1 - eps+)((a) + (a\')) <= g, eps+ = max(eps, 0): (g - bound)/g + 0.005', [(r['g'] - (1 - max(r['eps'], 0.0)) * (r['branch_a'] + r['branch_ap'])) / r['g'] + 0.005 for r in R], R)
    h2 = [r for r in R if r['H2']]
    stat('CR floor under (H2): ln(A 2 (pi^2 + sigma^2)/sigma^2)', [r['T1']['ln_1_over_A'] * -1 + np.log(2 * (np.pi ** 2 + r['sigma'] ** 2) / r['sigma'] ** 2) for r in h2], h2)
    bd = [r for r in R if r['BD1_bound'] > 0]
    stat('BD1 (paper F): (pi^2/2) theta - E_piv[d^2] (relative)', [(r['BD1_bound'] - r['E_piv_d2']) / r['BD1_bound'] for r in bd], bd, rel=1e-6)
    # summaries in case 3
    regs = {'text (A/(A+B) < 2 sigma)': lambda r: r['case3_text'] and r['A_over_2sigma'] < 1, 'op': lambda r: r['case3_op'],
            'both': lambda r: r['case3_op'] and r['case3_text'] and r['A_over_2sigma'] < 1,
            'text with (H1), (H2)': lambda r: r['case3_text'] and r['A_over_2sigma'] < 1 and r['H1'] and r['H2'],
            'op with (H1), (H2)': lambda r: r['case3_op'] and r['H1'] and r['H2']}
    out['case3'] = {}
    for nm, f in regs.items():
        sub = [r for r in R if f(r)]
        if not sub: out['case3'][nm] = {'n': 0}; continue
        def arg(key, mx=False):
            b = (max if mx else min)(sub, key=lambda r: abs(r[key]) if key == 'eps' else r[key])
            return {'value': float(abs(b[key]) if key == 'eps' else b[key]), 'family': b['fam'], 'sigma': b['sigma'], 'g_sqrt_l': b['g_sqrt_l'], 'F_c': b['F_c'],
                    'F_c2': b['F_c2'], 'label': b.get('label', b.get('source'))}
        o = {'n': len(sub), 'min_F_c': arg('F_c'), 'min_F_c2': arg('F_c2'), 'min_g_sqrt_l': arg('g_sqrt_l'), 'min_AT2_over_g': arg('AT2_over_g'),
             'min_AT2b_over_g': arg('AT2b_over_g'), 'max_abs_eps': arg('eps', True), 'min_sqrt_factor': arg('sqrt_factor')}
        bys = {}
        for r in sub:
            k = '%g' % r['sigma']; b = bys.setdefault(k, {'n': 0, 'min_g_sqrt_l': np.inf, 'min_F_c': np.inf, 'min_F_c2': np.inf})
            b['n'] += 1; b['min_g_sqrt_l'] = min(b['min_g_sqrt_l'], r['g_sqrt_l']); b['min_F_c'] = min(b['min_F_c'], r['F_c']); b['min_F_c2'] = min(b['min_F_c2'], r['F_c2'])
        o['by_sigma'] = dict(sorted(bys.items(), key=lambda kv: -float(kv[0])))
        byf = {}
        for r in sub:
            b = byf.setdefault(r['fam'], {'n': 0, 'min_g_sqrt_l': np.inf, 'min_F_c': np.inf, 'min_F_c2': np.inf, 'min_AT2b_over_g': np.inf})
            b['n'] += 1
            for k in ('g_sqrt_l', 'F_c', 'F_c2', 'AT2b_over_g'): b['min_' + k] = min(b['min_' + k], r[k])
        o['by_family'] = byf; out['case3'][nm] = o
    # region R4 of the conditional theorem (c^2 = 0.1, K = 3200, sigma <= 1e-4)
    r4 = [r for r in R if r['A_over_2sigma'] < 1 and r['theta_over_case1'] < 1 and r['T1']['B'] >= K_R4 * r['sigma'] and r['sigma'] <= 1e-4]
    out['R4'] = {'K': K_R4, 'n': len(r4), 'min_C_w2': float(min(r['C_w2'] for r in r4)) if r4 else None,
                 'min_Cw2_lower': float(min(r['Cw2_lower'] for r in r4)) if r4 else None}
    small_B = [r for r in R if r['A_over_2sigma'] < 1 and r['theta_over_case1'] < 1 and r['T1']['B'] < K_R4 * r['sigma']]
    out['R3_small_B'] = {'n': len(small_B), 'min_g_sqrt_l': float(min(r['g_sqrt_l'] for r in small_B)) if small_B else None}
    out['ok'] = bool(all(v['holds'] for v in out['lemmas'].values())); out['seconds'] = time.time() - t0
    C.save('ee_checks.json', out, RES)
    for k, v in out['lemmas'].items(): print('%-62s n %5d  min margin %+.4g  (%s, sigma %g)  holds %s' % (k, v['n'], v['min_margin'], v['family'], v['sigma'], v['holds']))
    for k, v in out['case3'].items():
        if v['n']: print('case 3 %-28s n %5d: min F_c %.4g, min F_c\' %.4g, min g sqrt(l) %.4f (%s, sigma %g), min AT2/g %.4f, min AT2\'/g %.4f, max |eps| %.2e' % (
            k, v['n'], v['min_F_c']['value'], v['min_F_c2']['value'], v['min_g_sqrt_l']['value'], v['min_g_sqrt_l']['family'], v['min_g_sqrt_l']['sigma'],
            v['min_AT2_over_g']['value'], v['min_AT2b_over_g']['value'], v['max_abs_eps']['value']))
    print('R4 (K = %g): %d states, min C_w\' %s, min lower bound %s; small-B part: %d states' % (K_R4, out['R4']['n'], out['R4']['min_C_w2'], out['R4']['min_Cw2_lower'], out['R3_small_B']['n']))
    print('states %d (numerically broken, excluded: %d); all checks hold: %s' % (len(R), out['n_numerically_broken'], out['ok']))


if __name__ == '__main__':
    np.seterr(all='ignore'); main()
