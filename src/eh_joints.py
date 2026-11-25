"""Stage EH, 1.2: the joints of the chain of Theorem EG on recorded states.  At each joint the conclusion of the earlier item is evaluated on the state
and inserted into the hypothesis of the next item; the margin is reported (relative: (right - left)/|right| or a ratio, as labelled).
Item numbers [n] are those of '定理EG_通しの証明.md'; the joints are numbered J1, J2, ... in the order printed.  States: results_EH/eh_states.json (engine: stage EE 445 + stage EG 247),
results_EF/ef_q2_2d.json (2-D closed-form families, 2449), results_EG/eg_circ_*.json, results_EG/eg_direct.json.
Output results_EH/eh_joints.json.  Python 3.11 grammar."""
import os, sys, json, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import ed_common as C

C1, C2, C3, CSTAR, CO, KR, CC = 7.77, 2.87, 2.03, 80.5, 112.0, 400.0, math.sqrt(0.1)
J = []


def psi(L): return 289 * L ** 2.5 / math.log(2 * L) ** 1.5 + 4 * L + 3.57
def alpha(s): return 1 / (2 * (math.pi ** 2 + s * s))
def kbar(s): return math.sqrt(psi(2 * math.log(1 / s) + math.log(1 / alpha(s))))
def h(x): L = math.log(1 / x); return 562 * x * L ** 3 / math.log(2 * L) ** 2


def joint(items, what, vals, tol=0.0, unit='relative margin (>= 0 holds)', note=''):
    v = [x for x in vals if x is not None and np.isfinite(x)]
    r = {'J': 'J%d' % (len(J) + 1), 'items': items, 'what': what, 'n': len(v), 'min_margin': float(min(v)) if v else None, 'holds': bool(v and min(v) >= -tol), 'tol': tol, 'unit': unit, 'note': note}
    J.append(r); print('J%-3d %-16s %-92s n %4d  min %+.4g  %s' % (len(J), items, what, r['n'], r['min_margin'] if v else float('nan'), 'ok' if r['holds'] else '** BREAKS **'))
    return r


def main():
    S = json.load(open(os.path.join(C.RES_EH, 'eh_states.json'), encoding='utf-8'))['rows']
    for r in S:
        t = r['T1']; s = r['sigma']; r['A'] = t['A']; r['B'] = t['B']; r['th'] = max(r['theta'], 0.0); r['M'] = r['A'] + r['th'] * r['B'] / 2; r['inO'] = r['M'] < CO * s * s
        r['den'] = s * (r['A'] + r['B']); r['main'] = r['den'] * (r['branch_a'] + r['branch_ap'])
        # B_F + R3 as the engine computes it (paper F's at_tools: eps = -(B_F + R3)/D_F, D_F = 2 sigma (A + B) sqrt((a)(a')));  the branch identity
        # sigma (A + B) g = 2 A_w + N_phi + (B_F + R3) then holds up to the engine's quadrature defect ('branch residual', <= 0.5 % by the accuracy filter)
        r['DF'] = 2 * r['den'] * math.sqrt(max(r['branch_a'] * r['branch_ap'], 0.0)); r['BR'] = -r['eps'] * r['DF']; r['branch_defect'] = abs(r['den'] * r['g'] - r['main'] - r['BR']) / (r['den'] * r['g'])
        r['L'] = math.log(1 / r['A'])
    out = {'n_engine': len(S), 'n_engine_in_O': sum(r['inO'] for r in S)}
    # [1 -> 4] EF-1 on every engine state (inside O too: new in this stage)
    joint('[4]', 'branch identity [4] in the engine: 1 - |sigma(A+B) g - (2A_w + N_phi) - (B_F + R3)|/(sigma(A+B) g)/0.005 (the accuracy filter of EE)', [1 - r['branch_defect'] / 0.005 for r in S])
    joint('[5]->[6]', 'EF-1: |B_F + R3| <= 7.77 s^2 + 2.87 s sqrt(theta B) + 2.03 s sqrt(A)   (B_F + R3 = -eps D_F of the engine)',
          [1 - abs(r['BR']) / (C1 * r['sigma'] ** 2 + C2 * r['sigma'] * math.sqrt(r['th'] * r['B']) + C3 * r['sigma'] * math.sqrt(r['A'])) for r in S])
    # [4 -> 5] threshold
    joint('[6]->[7]', 'EF-2 (outside O, M >= 112 s^2): |B_F + R3| <= M/2  (B_F + R3 = -eps D_F)', [(r['M'] / 2 - abs(r['BR'])) / (r['M'] / 2) for r in S if not r['inO']])
    joint('[7]->[23]', 'EF-2 conclusion: g / ((a) + (a\')) >= 1/2 outside O', [r['g'] / (r['branch_a'] + r['branch_ap']) - 0.5 for r in S if not r['inO']], unit='g/((a)+(a\')) - 1/2')
    joint('[4](b)', 'M <= 2 A_w + N_phi (A_w >= A/2, N_phi = theta B/2): sigma (A + B)((a) + (a\')) - M >= 0', [(r['main'] - r['M']) / r['main'] for r in S], tol=1e-9)
    # [6] ED lemma 0: paper F's k_eff^2 (X_int) = K/A (ED's sums): the same quantity in two codes
    joint('[9]', 'k_eff^2 of paper F (X_int, at_tools) = K/A of the ED sums: 1 - |rel. difference|/1e-3', [1 - abs(r['k_eff2'] / r['T1']['k_eff2'] - 1) / 1e-3 for r in S], unit='1 - |rel diff|/1e-3')
    # [8] ED lemma 3 identity on the engine
    rk = [abs(r['T1']['K'] - (r['T1']['D'] + r['A'] / 2 - r['T1']['SV'] / 2 - r['T1']['S2'] / 4)) / max(r['T1']['K'], r['A'] / 2) for r in S]
    joint('[11]', 'K = D + A/2 - SV/2 - S2/4 (ED lemma 3; a proved equality, checked as the accuracy of the sums): 1 - |res|/(1e-2 max(K, A/2))',
          [1 - x / 1e-2 for x in rk], unit='1 - |res|/(1e-2 max(K, A/2))',
          note='relative to K alone the residual reaches 0.79 on states with K/A ~ 1e-16 (cancellation of A/2 - SV/2); relative to max(K, A/2) it is <= 1e-3 except on the two states of the two-straight-points family at sigma = 1e-36 (0.86 %)')
    out['identity_K_worst'] = sorted([(x, r['source'], r['family'], r['sigma']) for x, r in zip(rk, S)], reverse=True)[:4]
    out['n_identity_K_rel_gt_1e-3'] = sum(1 for r in S if abs(r['T1']['K'] - (r['T1']['D'] + r['A'] / 2 - r['T1']['SV'] / 2 - r['T1']['S2'] / 4)) > 1e-3 * r['T1']['K'])
    # [7] half lemma
    joint('[10]', 'ED half lemma: S2 + SV <= 8 A L + 7.136 A  (engine)', [(8 * r['A'] * r['L'] + 7.136 * r['A'] - r['T1']['S2'] - r['T1']['SV']) / (8 * r['A'] * r['L'] + 7.136 * r['A']) for r in S])
    q2 = json.load(open(os.path.join(C.RES_EF, 'ef_q2_2d.json'), encoding='utf-8'))['rows']; R2 = [x for f in q2.values() for x in f]
    joint('[10]', 'ED half lemma, SV part: SV <= 8 A L + 7.136 A  (2-D closed-form families)', [(8 * x['A'] * x['ln_1_over_A'] + 7.136 * x['A'] - x['SV']) / (8 * x['A'] * x['ln_1_over_A'] + 7.136 * x['A']) for x in R2])
    # [9],[15] -> k_eff^2: (Q2 + SV/2)/A <= psi(max(L, 4)) on the 2-D families; K/A <= psi on the engine
    joint('[12,18]->[23]', 'Theorem EG-1 (iv): (Q2 + SV/2)/A <= psi(max(L,4))  (2-D families; Q2, SV from the sums)', [1 - (x['Q2'] + x['SV'] / 2) / x['A'] / psi(max(x['ln_1_over_A'], 4)) for x in R2])
    joint('[18]->[23]', 'Theorem EG-1 (v): k_eff^2 = K/A <= psi(max(L, 4))  (engine)', [1 - r['T1']['k_eff2'] / psi(max(r['L'], 4)) for r in S])
    joint('[12]', 'K <= Q2 + SV/2 is used through D = Q2 + SV - A; on the engine: K - D - A/2 + SV/2 = -S2/4 <= 0', [(r['T1']['S2'] / 4) / r['T1']['K'] for r in S], tol=1e-12, unit='(S2/4)/K >= 0')
    # [16] Cramer-Rao floor and L <= L_bar, k_eff <= k_bar
    H2 = [r for r in S if r['H2']]
    joint('[8]', 'EE-3 (states with (H2)): A >= alpha sigma^2', [r['A'] / (alpha(r['sigma']) * r['sigma'] ** 2) - 1 for r in H2], unit='A/(alpha s^2) - 1')
    joint('[8]->[28]', 'L = ln(1/A) <= L_bar = 2 ln(1/sigma) + ln(1/alpha)  ((H2))', [(2 * math.log(1 / r['sigma']) + math.log(1 / alpha(r['sigma'])) - r['L']) for r in H2], unit='L_bar - L')
    joint('[8,18]->[22,23]', 'k_eff <= k_bar(sigma)  ((H2))', [1 - math.sqrt(r['T1']['k_eff2']) / kbar(r['sigma']) for r in H2])
    out['n_H2'] = len(H2); out['n_no_H2'] = len(S) - len(H2)
    # [19] improved EE-C and [20] EF-3(b)
    joint('[20]', 'EE-C improved: C~_w >= 1 - (A + 1.291 sqrt A)/B - 8.4258 sqrt(theta/B)', [r['C_w2'] - (1 - (r['A'] + 1.291 * math.sqrt(r['A'])) / r['B'] - 8.4258 * math.sqrt(r['th'] / r['B'])) for r in S], unit='difference')
    joint('[19]->[21]', 'EF-3(b): sigma^2 C~_w^2 B^2 <= 16 theta B k_eff^2 A  (the factor sigma^2 is the joint of EF)', [1 - (r['sigma'] * r['C_w2'] * r['B']) ** 2 / (16 * r['th'] * r['B'] * r['k_eff2'] * r['A']) for r in S], tol=2e-5,
          note='paper F dictionary states are equality cases; tolerance = quadrature')
    joint('[19]', 'the same without sigma^2 (as first written in the instruction of stage EF): breaks', [1 - (r['C_w2'] * r['B']) ** 2 / (16 * r['th'] * r['B'] * r['k_eff2'] * r['A']) for r in S if r['th'] > 0][:1] and
          [min(1 - (r['C_w2'] * r['B']) ** 2 / (16 * r['th'] * r['B'] * r['k_eff2'] * r['A']) for r in S if r['th'] > 0)], unit='power check: expected to break')
    J[-1]['expected_to_break'] = True; J[-1]['holds_as_power_check'] = not J[-1]['holds']
    # [21] EF-O inside O
    O = [r for r in S if r['inO']]
    joint('[22]', 'EF-O inside O: B <= sigma max(280, 1268 k_eff)', [1 - r['B'] / (r['sigma'] * max(280, 1268 * math.sqrt(r['k_eff2']))) for r in O])
    joint('[22,8,18]->[28]', 'inside O ((H2)): A + B <= 112 sigma^2 + 1268 sigma k_bar', [1 - (r['A'] + r['B']) / (CO * r['sigma'] ** 2 + 1268 * r['sigma'] * kbar(r['sigma'])) for r in O if r['H2']])
    # [24] Theorem EF (I): regions, with the actual k_eff and with k_bar
    reg = {1: [], 2: [], 3: [], 4: []}; low_k = []; low_kb = []
    for r in S:
        if r['inO'] or r['sigma'] > 1e-4: continue
        s = r['sigma']; l = math.log(1 / s); A, B = r['A'], r['B']
        if A / (A + B) >= 2 * s: k = 1; b = 1.0
        elif r['th'] >= 2 * s * CC / math.sqrt(l): k = 2; b = CC / (2 * math.sqrt(l))
        elif B < KR * s: k = 3; b = 0.14 * (1 - 2 * s)
        else:
            k = 4; bp = 2 * s / (1 - 2 * s) + 1.291 * math.sqrt(2 / ((1 - 2 * s) * KR)) + 8.4258 * math.sqrt(2 * CC / (KR * math.sqrt(l)))
            b = (1 - 2 * s) * (1 - bp) / (4 * math.sqrt(2) * math.sqrt(r['k_eff2'])); r['b_kbar'] = (1 - 2 * s) * (1 - bp) / (4 * math.sqrt(2) * kbar(s)) if r['H2'] else None
            r['Cw_ge_1_minus_bp'] = r['C_w2'] - (1 - bp)
        reg[k].append(r); low_k.append(r['g'] / b - 1)
        if r['H2']:
            tb = min(1.0, CC / (2 * math.sqrt(l)), 0.14 * (1 - 2 * s), (1 - 2 * s) * (1 - 0.284) / (4 * math.sqrt(2) * kbar(s))); low_kb.append(r['g'] / tb - 1)
    joint('[23]', 'Theorem EF (I) region by region, with the actual k_eff: g >= region bound', low_k, unit='g/bound - 1')
    joint('[23]->[28]', 'Theorem EG (a) (outside O, sigma <= 1e-4, (H2)): g >= min{1, c/(2 sqrt l), 0.14(1-2s), (1-2s)(1-0.284)/(4 sqrt2 k_bar)}', low_kb, unit='g/bound - 1')
    joint('[20]->[23]', 'region 4: C~_w >= 1 - beta\'', [r['Cw_ge_1_minus_bp'] for r in reg[4]], unit='difference')
    out['regions'] = {k: len(v) for k, v in reg.items()}
    # [25] EG-D on the engine (record of stage EG) and the plain checks of stage EG
    d = json.load(open(os.path.join(C.RES_EG, 'eg_direct.json'), encoding='utf-8')); pt = json.load(open(os.path.join(C.RES_EG, 'eg_point.json'), encoding='utf-8'))
    joint('[24]', 'EG-D on two Gaussian blobs (plain quadrature, stage EG): 1 - |R\'|/bound', [1 - abs(x['R']) / x['bound'] for x in pt['f']])
    joint('[24]', 'EG-D on the three points, sigma = 1 .. 0.1 (stage EG): 1 - |R\'|/bound', [1 - abs(x['R']) / x['bound'] for x in pt['d']])
    joint('[24]', 'EG-D on 247 engine states: 1 - |R_corr|/(bound + engine tolerance)  (record of stage EG)', [1 - abs(x['R_corr']) / (x['bound'] + x['tolerance']) for x in d['rows']])
    joint('[24]', 'EG-D: A_w0 >= A/2 and |E2| <= (sigma/2) sqrt(U) ||W_zz||', [min(x['E2_CS_margin'] / (1e-30 + abs(x['E2']) + 1e-300), 1.0) if x['Aw0_ge_half_A'] else -1 for x in d['rows']], unit='>= 0')
    # [26, 27] -> [28]: ||W_zz||^2 <= h(J) for J <= e^-4
    rw = [x for x in d['rows'] if x['J'] <= math.exp(-4)]
    joint('[25,26]->[27]', 'EG-iso + EG-sm: ||W_zz||_2^2 <= h(A + B)  (states with A + B <= e^-4)', [1 - x['W_zz_L2'] ** 2 / h(x['J']) for x in rw])
    out['n_J_le_e4'] = len(rw)
    # [28] EG-O on states with A + B <= 2e-9 (and (H2))
    by = {(r['source'][3:], round(r['g'], 9)): r for r in S if r['source'].startswith('EG')}
    rO_all = [x for x in d['rows'] if x['J'] <= 2e-9]; rO = [x for x in rO_all if by.get((x['record'], round(x['g'], 9)), {}).get('H2')]; rOx = [x for x in rO_all if x not in rO]
    joint('[27]', 'EG-O (A + B <= 2e-9, (H2)): |r0| <= 1/4  (r0 with the engine defect taken out)', [0.25 - abs(x['r0_corr']) for x in rO], unit='1/4 - |r0|')
    joint('[27]', 'EG-O (A + B <= 2e-9, (H2)): g >= (3/4) A/(sigma (A + B))', [x['g'] * x['sigma'] * x['J'] / (0.75 * x['A_over_s2'] * x['sigma'] ** 2) - 1 for x in rO], unit='g/bound - 1')
    joint('[27]', 'outside the hypotheses ((H2) fails; %d states): |r0| <= 1/4 all the same' % len(rOx), [0.25 - abs(x['r0_corr']) for x in rOx], unit='1/4 - |r0|')
    out['n_J_le_2e9'] = len(rO_all); out['n_J_le_2e9_H2'] = len(rO); out['sigma_of_J_le_2e9'] = sorted(set(x['sigma'] for x in rO))
    # [15] in 2-D and 1-D: the circles (records of stage EG)
    for p in ('2d', 'iso', '1d'):
        c = json.load(open(os.path.join(C.RES_EG, 'eg_circ_%s.json' % p), encoding='utf-8'))
        joint('[14-18]', 'eg_circ %s: every step of the three-circle proof (stage EG record): min ln(right/left)' % p, [v for v in c['all_min_margins'].values() if v is not None], unit='ln(right/left)')
    out['joints'] = J; out['env'] = C.env()
    out['ok'] = bool(all(j['holds'] or j.get('expected_to_break') for j in J) and all(j['holds_as_power_check'] for j in J if j.get('expected_to_break')))
    C.save('eh_joints.json', out, C.RES_EH)
    print('engine states %d (inside O %d; (H2) %d); regions outside O, sigma <= 1e-4: %s; states with A + B <= e^-4: %d, <= 2e-9: %d (sigma %s); ok %s' % (
        out['n_engine'], out['n_engine_in_O'], out['n_H2'], out['regions'], out['n_J_le_e4'], out['n_J_le_2e9'], out['sigma_of_J_le_2e9'], out['ok']))


if __name__ == '__main__':
    np.seterr(all='ignore'); main()
