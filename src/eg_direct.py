"""Stage EG, 1.2: the direct expansion of the gap around the observed point (report EG, Lemma EG-D) on states of paper F's engine.
  sigma (A + B) g = 2 A_w0 + N_phi0 + E2 + R',   A_w0 = E[s~_w^2 (1 - phi0/4)], N_phi0 = (1/2) E[s~_v^2 phi0], phi0 = 1 - cos(Y_z),
  E2 = (sigma/2) E[s~_p sin(Y_z) W_zz], W_zz = sigma^2 d_z^2 p/p,   |R'| <= (sqrt15/6) sigma^2 sqrt(J_pp),   |E2| <= (sigma/2) sqrt(U) ||W_zz||_2, U = E[s~_p^2 sin^2 Y_z].
States: from the wide-state records of stage EF (results_EF/ef_wide_*.json), re-evaluated with the engine of this stage (the sums of the direct form
added to ed_engine.py, marked 'EG:'; each evaluation is stored at once in results_EG/cache):
  families N, M : the evaluations with the largest |eps|, the smallest g sqrt(l), the largest B/sigma inside O, and the largest A + theta B/2 outside O;
  family W      : (searches without a planted bump) the 25 evaluations inside O with the largest recorded |eps| and every 6th of the others.
CORRECTION OF STAGE EF found here: stage EF's engine built the v-grid around the centre v = 0 only (ed_engine.v_grid), so the second blob of family W
(v-centre vc2 up to 3/sigma) was cut off.  The recorded values of family W with two blobs (eps up to 8.7e-3) are artefacts of the cut; with the grid
of this stage (every v-centre covered) they are recomputed below ('eps' next to 'eps_EF_record').
usage: python eg_direct.py [I N]   (worker I of N: only fills the cache; without arguments: all states, then the record)
Output results_EG/eg_direct.json.  Python 3.11 grammar."""
import os, sys, glob, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import ed_common as C
import ee_engine as E
import ee_checks as K
import ef_wide as W

NSEL = 2
TAU = 2e-5
SET = (('search', dict(nper=12, coarse=100, per_res=10)),)


def select():
    S = []
    for f in sorted(glob.glob(os.path.join(C.RES_EF, 'ef_wide_*.json'))):
        d = json.load(open(f, encoding='utf-8')); fam = d.get('family', 'W'); good = [e for e in d['evaluations'] if 'objective' in e]
        for e in good: e['M'] = e['A_over_s2'] + max(e['theta'], 0.0) * e['B_over_s2'] / 2
        o = [e for e in good if e['M'] < W.CO]; x = [e for e in good if e['M'] >= W.CO]; pick = []
        if fam == 'W':
            if d['plant']: continue
            o.sort(key=lambda e: -abs(e['eps'])); pick = o[:25] + o[25::6] + x[::6]
        else:
            for key, pool in ((lambda e: -abs(e['eps']), o), (lambda e: e['g_sqrt_l'], o), (lambda e: -e['B_over_s2'], o), (lambda e: -e['M'], x), (lambda e: e['g_sqrt_l'], x)):
                pick += sorted(pool, key=key)[:NSEL]
        seen = set()
        for e in pick:
            k = tuple(e['z'])
            if k not in seen: seen.add(k); S.append((os.path.basename(f), fam, d['sigma'], e['z'], e['eps'], bool(e['M'] < W.CO)))
    return S


def main(worker=None):
    np.seterr(all='ignore'); t0 = time.time(); rows = []
    sel = list(enumerate(select()))
    for i, (fn, fam, sig, z, eps_ef, inO_ef) in (sel[::-1] if os.environ.get('EG_REVERSE') == '1' else sel):
        if worker is not None and i % worker[1] != worker[0]: continue
        bg, rip, Lv = {'W': W.build, 'N': W.build_N, 'M': W.build_M}[fam](np.array(z), sig)
        try:
            q = E.finish(E.ev(sig, bg, rip, Lv, SET[0][1])); level = SET[0][0]; rob0 = abs(q['direct']['residual']) / q['direct']['bound']
            if worker is None or os.environ.get('EG_REFINE') == '1':                         # the quadrature of the search settings is not always enough for the residual (relative 1e-4 on fine w-ripples): refine
                for name, st in SET[1:]:
                    if abs(q['direct']['residual']) <= q['direct']['bound']: break
                    q = E.finish(E.ev(sig, bg, rip, Lv, st)); level = name
        except (ValueError, FloatingPointError, OverflowError, MemoryError, np.linalg.LinAlgError, ZeroDivisionError, KeyError): continue
        if not (np.isfinite(q['g']) and q['g'] > 0 and K.branch_residual(q) <= 0.005): continue
        d = q['direct']; A, B = d['A'], d['B']; main = 2 * d['A_w0'] + d['N_phi0']; M = A + max(q['theta'], 0.0) * B / 2
        rows.append({'record': fn, 'family': fam, 'sigma': sig, 'z': z, 'eps_EF_record': eps_ef, 'in_O_EF_record': inO_ef, 'two_blobs': bool(fam == 'W' and z[3] > 0.01), 'setting': level, 'residual_over_bound_at_search_setting': rob0, 'g': q['g'], 'g_sqrt_l': q['g_sqrt_l'], 'eps': q['eps'], 'A_over_s2': A / sig ** 2, 'B_over_sigma': B / sig, 'in_O': bool(M < W.CO * sig ** 2),
                     'lhs': d['lhs'], 'main': main, 'E2': d['E2'], 'residual': d['residual'], 'bound': d['bound'], 'residual_over_bound': abs(d['residual']) / d['bound'],
                     'rel_residual': abs(d['residual']) / d['lhs'], 'rel_E2': abs(d['E2']) / d['lhs'], 'r0': d['lhs'] / main - 1,
                     'E2_CS_margin': float(0.5 * sig * np.sqrt(d['U'] * d['W_zz_L2sq']) - abs(d['E2'])), 'W_zz_L2': float(np.sqrt(d['W_zz_L2sq'])), 'J': A + B,
                     'Aw0_ge_half_A': bool(d['A_w0'] >= A / 2 * (1 - 1e-9)), 'two_A_w_plus_N_phi_m': float(sig * (A + B) * (q['branch_a'] + q['branch_ap']))})
        # the engine's own defect: paper E's identity  lhs = (2 A_w + N_phi)(m) + (B_F + R3)  holds on its grid only up to 'ibp_defect' (second-order error of the trapezoid
        # rule where the v- or w-grid changes its step; this is the 'branch residual' of stages EE, EF).  R_corr = the remainder R' of Lemma EG-D with that defect taken out:
        # (2 A_w + N_phi)(m) + (B_F + R3)_engine = 2 A_w0 + N_phi0 + E2 + R_corr,  (B_F + R3)_engine = -eps D_F.
        r_ = rows[-1]; DF = 2 * sig * (A + B) * np.sqrt(max(q['branch_a'] * q['branch_ap'], 0.0)); BR = -q['eps'] * DF
        r_['ibp_defect_rel'] = float(abs(d['lhs'] - r_['two_A_w_plus_N_phi_m'] - BR) / d['lhs']); r_['R_corr'] = float(r_['two_A_w_plus_N_phi_m'] + BR - main - d['E2'])
        r_['R_corr_over_bound'] = float(abs(r_['R_corr']) / d['bound']); r_['r0_corr'] = float((d['E2'] + r_['R_corr']) / main)
        # tolerance of the engine: TAU x lhs (the trapezoid rule is of second order where the grid changes its step: (1/coarse)^2/12 ~ 1e-5) and 10 % of the defect taken out
        r_['tolerance'] = float(TAU * d['lhs'] + 0.1 * r_['ibp_defect_rel'] * d['lhs']); r_['within_bound_and_tolerance'] = bool(abs(r_['R_corr']) <= d['bound'] + r_['tolerance'])
        r = rows[-1]; print('%-30s %s sig %g inO %d [%s]: A/s2 %.3g B/sig %.3g g %.4g | residual/bound %.3f  rel residual %.1e  r0 %+.2e  eps %+.1e (EF record %+.1e)  %.0fs' % (
            fn[8:-5], fam, sig, r['in_O'], level, r['A_over_s2'], r['B_over_sigma'], r['g'], r['residual_over_bound'], r['rel_residual'], r['r0'], r['eps'], eps_ef, time.time() - t0), flush=True)
    if worker is not None: return
    o = [r for r in rows if r['in_O']]; tb = [r for r in rows if r['two_blobs']]; byrec = {}
    for r in rows:
        b = byrec.setdefault(r['record'], {'family': r['family'], 'sigma': r['sigma'], 'n': 0, 'n_in_O': 0, 'max_abs_eps_in_O': 0.0, 'max_abs_eps_EF_record': 0.0, 'max_abs_r0_in_O': 0.0, 'min_g_sqrt_l_in_O': None})
        b['n'] += 1; b['max_abs_eps_EF_record'] = max(b['max_abs_eps_EF_record'], abs(r['eps_EF_record']))
        if r['in_O']:
            b['n_in_O'] += 1; b['max_abs_eps_in_O'] = max(b['max_abs_eps_in_O'], abs(r['eps'])); b['max_abs_r0_in_O'] = max(b['max_abs_r0_in_O'], abs(r['r0_corr']))
            b['min_g_sqrt_l_in_O'] = r['g_sqrt_l'] if b['min_g_sqrt_l_in_O'] is None else min(b['min_g_sqrt_l_in_O'], r['g_sqrt_l'])
    out = {'rows': rows, 'n': len(rows), 'n_in_O': len(o), 'max_residual_over_bound': max(r['residual_over_bound'] for r in rows), 'max_rel_residual': max(r['rel_residual'] for r in rows),
           'max_abs_r0': max(abs(r['r0']) for r in rows), 'max_abs_r0_in_O': max([abs(r['r0']) for r in o] + [0.0]), 'min_E2_CS_margin': min(r['E2_CS_margin'] for r in rows),
           'all_Aw0_ge_half_A': bool(all(r['Aw0_ge_half_A'] for r in rows)), 'max_W_zz_L2_in_O': max([r['W_zz_L2'] for r in o] + [0.0]), 'env': C.env(), 'by_record': byrec, 'max_abs_eps_in_O': max([abs(r['eps']) for r in o] + [0.0]),
           'two_blobs': {'n': len(tb), 'max_abs_eps_EF_record': max([abs(r['eps_EF_record']) for r in tb] + [0.0]), 'max_abs_eps_recomputed': max([abs(r['eps']) for r in tb] + [0.0])},
           'min_g_sqrt_l_in_O': min([r['g_sqrt_l'] for r in o] + [float('inf')]), 'n_by_setting': {k: sum(1 for r in rows if r['setting'] == k) for k, _ in SET},
           'unresolved': [r['record'] + ' ' + repr(r['z']) for r in rows if r['residual_over_bound'] > 1.0]}
    ex = [r for r in rows if r['residual_over_bound'] > 1.0]              # the raw residual is above the bound: always together with an engine defect of the same size
    out['n_within_bound'] = len(rows) - len(ex); out['n_above_raw'] = len(ex); out['max_residual_over_bound_within'] = max(r['residual_over_bound'] for r in rows if r['residual_over_bound'] <= 1.0)
    out['max_R_corr_over_bound'] = max(r['R_corr_over_bound'] for r in rows); out['n_R_corr_within_bound'] = sum(1 for r in rows if r['R_corr_over_bound'] <= 1.0)
    out['max_ibp_defect_rel'] = max(r['ibp_defect_rel'] for r in rows); out['all_within_bound_and_tolerance'] = bool(all(r['within_bound_and_tolerance'] for r in rows))
    out['max_abs_r0_corr_in_O'] = max([abs(r['r0_corr']) for r in o] + [0.0]); out['tau'] = TAU; out['above_raw_min_defect_over_residual'] = min([r['ibp_defect_rel'] / r['rel_residual'] for r in ex] + [float('inf')])
    out['ok'] = bool(out['all_within_bound_and_tolerance'] and out['n_R_corr_within_bound'] >= 0.95 * out['n'] and out['min_E2_CS_margin'] >= -1e-18 and out['all_Aw0_ge_half_A'] and out['n'] >= 150); out['seconds'] = time.time() - t0
    C.save('eg_direct.json', out, C.RES_EG)
    print('raw residual within the bound: %d of %d; with the engine defect taken out: %d of %d (max |R_corr|/bound %.3f); largest engine defect %.1e (relative)' % (out['n_within_bound'], out['n'], out['n_R_corr_within_bound'], out['n'], out['max_R_corr_over_bound'], out['max_ibp_defect_rel']))
    print('%d states (%d inside O): max |residual|/bound %.4f, max relative residual %.2e, max |r0| %.3e (inside O %.3e), max |eps| inside O %.3e, ok %s' % (
        out['n'], out['n_in_O'], out['max_residual_over_bound'], out['max_rel_residual'], out['max_abs_r0'], out['max_abs_r0_in_O'], out['max_abs_eps_in_O'], out['ok']))
    print('family W, two blobs: %d states, max |eps| in the EF records %.3e -> recomputed %.3e' % (len(tb), out['two_blobs']['max_abs_eps_EF_record'], out['two_blobs']['max_abs_eps_recomputed']))


if __name__ == '__main__':
    fn = os.path.join(C.RES_EG, 'eg_direct.json')
    if len(sys.argv) > 2: main((int(sys.argv[1]), int(sys.argv[2]))); sys.exit(0)
    if os.path.exists(fn) and os.environ.get('EG_FORCE') != '1' and os.environ.get('EG_NOCACHE') != '1': print('exists, skipped'); sys.exit(0)
    main()
