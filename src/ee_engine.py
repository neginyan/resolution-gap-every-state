"""Stage EE: one evaluation of a state with every factor of F_c, the true gap g, g sqrt(l), the case flags and the T1 sums.
On top of stage ED's engine (ed_engine.evaluate: paper F's exact point fields, v3_tools / at_tools of resolution-gap-two-sided v1.5, unchanged).
Added point field (paper F's point_fields2 wrapped, the release is not modified):
  'Cw2_int' = s_v^2 cos(m_z) (1 - sigma^2 (H_ww + H_vw))   -> C_w' = E[Cw2_int]/J_vv
            (the window inequality with psi = sqrt2 sin(m_z) in place of 2 sin(m_z/2): psi^2 = 2 phi (2 - phi) <= 4 phi, d_w psi = cos(m_z) (1 - sigma^2(H_ww + H_vw));
             the same Cauchy-Schwarz step gives (C_w' J_vv)^2 <= 16 E[s_v^2 phi] int (d_w f)^2, so Theorem AT2 holds with |C_w'| in place of |C_w|.)
Outputs (all floats; 'case3_op' is bd_case3.row's operational test, 'case3_text' the text's case split with c^2 = 0.1):
  g, g_sqrt_l, F_c = (a/den)(1 - eps)|C_w| sqrt((1 - Phi/4)(1 - delta_w)), its factors, AT2 = F_c/(2 k_eff), AT2/g; the same with C_w';
  the branches (a), (a') and the branch identity; theta, Q, W2, BB1; the T1 sums (sigma = 1 units) A, B, D, K, S2, SV.
Python 3.11 grammar."""
import os, sys, json, hashlib
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import ed_common as C
import ed_engine as G
import at_tools                                   # paper F (put on the path by ed_engine -> C.use_F(); bd_case3 already added 'svd2')

_pf_prev = at_tools.point_fields2


def _pf_ee(pre, sig, X):
    p, q = _pf_prev(pre, sig, X); q = dict(q)
    q['Cw2_int'] = q['sv'] ** 2 * np.cos(q['mz']) * (1 - q['aHww'] - sig * sig * q['Hvw'])
    return p, q


if 'Cw2_int' not in at_tools.KT:
    at_tools.point_fields2 = _pf_ee; at_tools.KT = at_tools.KT + ('Cw2_int',)

C2 = 0.1
FAST = dict(nv=16, n0=11, nper=12, coarse=100)
FULL = dict(nv=32, n0=21, nper=40, coarse=300)
CACHE = os.environ.get('EH_CACHE') or os.environ.get('EG_CACHE') or os.path.join(C.ROOT, 'results_EH', 'cache')     # EH: the cache of this stage        # EF: the cache of this stage
NOCACHE = os.environ.get('EE_NOCACHE') == '1' or os.environ.get('EF_NOCACHE') == '1' or os.environ.get('EG_NOCACHE') == '1' or os.environ.get('EH_NOCACHE') == '1'


def _key(sig, bg, rip, Lv, o):
    def cv(x):
        if isinstance(x, dict): return {k: cv(v) for k, v in sorted(x.items())}
        if isinstance(x, (list, tuple)): return [cv(v) for v in x]
        if isinstance(x, (float, np.floating)): return float.hex(float(x))
        if isinstance(x, np.integer): return int(x)
        return x
    far = any(abs(b[1]) > 0 or (len(b) > 3 and abs(b[3]) > 0) for b in bg)          # states with a background off the origin: computed with the fixed terms()
    return hashlib.sha1(json.dumps(cv(['EG3' if far else 'EG2', sig, bg, rip, Lv, o]), sort_keys=True).encode()).hexdigest()


def var_xz(T):
    """Var(X_z), X_z = (x_v + x_w)/sqrt2, x ~ rho, exact from the terms (as ed_engine.W2_exact)."""
    M0 = M1 = M2 = 0.0
    for lc, m, S, kap in T:
        cf = np.exp(lc) * np.exp(-0.5 * kap @ S @ kap); mu = m + 1j * (S @ kap)
        bz = mu[0] + mu[1]; vb = S[0, 0] + 2 * S[0, 1] + S[1, 1]
        M0 += cf; M1 += cf * bz; M2 += cf * (vb + bz * bz)
    M0 = M0.real; m1 = (M1 / M0).real; m2 = (M2 / M0).real
    return float((m2 - m1 * m1) / 2), float(m1 / np.sqrt(2))


def quantities(sig, bg, rip, Lv_sig, **kw):
    o = G.evaluate(sig, bg, rip, Lv_sig, **kw); E, cm = o['E'], o['cm']
    r = at_tools.analyse(E, sig); s2 = sig * sig; ln = np.log(1 / sig)
    A, B = r['sig2Jww'], r['sig2Jvv']; den = A + B
    theta = (1 - r['delta_c']) * r['Phi']; W2 = G.W2_exact(sig, o['T']); Q = 4 * theta / W2
    keff = r['k_eff']; sq = np.sqrt((1 - r['Phi'] / 4) * (1 - r['delta_w']))
    Fc = r['a_over_den'] * (1 - r['eps']) * abs(r['C_w']) * sq
    Cw2 = E['Cw2_int'] / E['svv']; Fc2 = r['a_over_den'] * (1 - r['eps']) * abs(Cw2) * sq
    bb1 = (1 - abs(r['eps'])) * theta / (2 * sig)
    br_a = 2 * s2 * E['sww_w'] / (sig * den); br_ap = 0.5 * s2 * E['phi_svv'] / (sig * den)
    c = np.sqrt(C2)
    case1 = bool(theta >= 2 * sig * c / np.sqrt(ln)); case2 = bool(Q >= C2 / ln)
    case3_op = bool(Q * ln < 0.1 and bb1 < r['bound_AT2'] and A / den < 2 * sig)
    tm = cm['tm']; tA, tB, tD, tK, tS2, tSV = (cm[k] / tm for k in ('tA', 'tB', 'tD', 'tK', 'tS2', 'tSV'))
    g = float(r['g'])
    out = {'sigma': sig, 'ln': ln, 'g': g, 'g_sqrt_l': g * np.sqrt(ln), 'F_c': float(Fc), 'a_over_den': float(r['a_over_den']), 'eps': float(r['eps']),
           'C_w': float(r['C_w']), 'Phi': float(r['Phi']), 'delta_w': float(r['delta_w']), 'sqrt_factor': float(sq), 'k_eff2': float(keff ** 2),
           'AT2': float(r['bound_AT2']), 'AT2_over_g': float(r['bound_AT2'] / g) if g > 0 else float('nan'),
           'C_w2': float(Cw2), 'F_c2': float(Fc2), 'AT2b': float(Fc2 / (2 * keff)), 'AT2b_over_g': float(Fc2 / (2 * keff) / g) if g > 0 else float('nan'),
           'theta': float(theta), 'theta_over_case1': float(theta / (2 * sig * c / np.sqrt(ln))), 'Q_over_case2': float(Q * ln / C2), 'Q': float(Q), 'W2': float(W2), 'ln_W2': float(np.log(W2)), 'Q_ln': float(Q * ln), 'BB1': float(bb1), 'A_frac': float(A / den),
           'A_over_2sigma': float(A / den / (2 * sig)), 'branch_a': float(br_a), 'branch_ap': float(br_ap),
           'branch_lower': float((1 - abs(r['eps'])) * (br_a + br_ap)), 'branch_upper': float((1 + abs(r['eps'])) * (br_a + br_ap)),
           'case1': case1, 'case2': case2, 'case3_text': bool(not case1 and not case2), 'case3_op': case3_op,
           'E_piv_d2': float(E['svd2'] / E['svv']), 'BD1_bound': float(np.pi ** 2 / 2 * theta),
           'T1': {'A': float(tA), 'B': float(tB), 'D': float(tD), 'K': float(tK), 'S2': float(tS2), 'SV': float(tSV), 'max_Var_Zv': float(cm['tVmax']),
                  'k_eff2': float(tK / tA), 'ln_1_over_A': float(np.log(1 / tA))},
           'n_points': int(o['n_points'])}
    # EG: the direct expansion around Y (Lemma EG-D):  sigma (A + B) g = 2 A_w0 + N_phi0 + E2 + R',  |R'| <= (sqrt15/6) sigma^2 sqrt(J_pp)
    dA0, dN0, dE2, dJ, dW2, dU = (cm[k] / tm for k in ('eA0', 'eN0', 'eE2', 'eJpp', 'eW2', 'eU')); lhs = sig * den * g
    out['direct'] = {'A_w0': float(dA0), 'N_phi0': float(dN0), 'E2': float(dE2), 'J_pp': float(dJ), 'W_zz_L2sq': float(dW2), 'U': float(dU), 'lhs': float(lhs),
                     'residual': float(lhs - (2 * dA0 + dN0 + dE2)), 'bound': float(np.sqrt(15) / 6 * s2 * np.sqrt(dJ)), 'A': float(A), 'B': float(B)}
    vz, mzx = var_xz(o['T']); EQ = cm['qz'] / cm['q0']
    out.update({'Var_Xz': vz, 'mean_Xz': mzx, 'E_Q_yz2': float(EQ), 'H1': bool(EQ <= np.pi ** 2 + s2 / 2), 'H2': bool(vz <= np.pi ** 2),
                'H1_ratio': float(EQ / (np.pi ** 2 + s2 / 2)), 'H2_ratio': float(vz / np.pi ** 2)})
    # Lemma EE-C (the lower bound for C_w' in sigma = 1 units; checked in ee_fc.py): C_w' >= 1 - (S2 + 9.73 sqrt(A))/B - 8.43 sqrt(theta/B)
    out['Cw2_lower'] = float(1 - (tS2 + (2 * np.sqrt(np.sqrt(315.0)) + np.sqrt(15.0) / 3) * np.sqrt(tA)) / tB - np.sqrt(2 * np.sqrt(1260.0) * max(theta, 0.0) / tB))
    return out


def finish(q):
    """theta recomputed without rounding and the case flags redone (applied when the records are read).
    paper F's analyse() returns theta = (1 - delta_c) Phi with delta_c = 1 - E[s_v^2 phi]/(J_vv Phi): when theta/Phi < 1e-16 (v-information at a straight
    point, mass at a bent one, sigma <= 1e-16) this rounds to 0.  The branch (a') = sigma^2 E[s_v^2 phi]/(2 sigma den) keeps it exactly:
    theta = 2 sigma (a')/(a/den).  Q, BB1, the case flags and the lower bound of Lemma EE-C are recomputed from it."""
    q = dict(q); sig, ln = q['sigma'], q['ln']; th = 2 * sig * q['branch_ap'] / q['a_over_den']
    q['theta_paperF_formula'] = q['theta']; q['theta'] = float(th); Q = 4 * th / q['W2']; q['Q'] = float(Q); q['Q_ln'] = float(Q * ln)
    c = np.sqrt(C2); q['theta_over_case1'] = float(th / (2 * sig * c / np.sqrt(ln))); q['Q_over_case2'] = float(Q * ln / C2)
    q['BB1'] = float((1 - abs(q['eps'])) * th / (2 * sig)); q['case1'] = bool(th >= 2 * sig * c / np.sqrt(ln)); q['case2'] = bool(Q >= C2 / ln)
    q['case3_text'] = bool(not q['case1'] and not q['case2'])
    q['case3_op'] = bool(Q * ln < 0.1 and q['BB1'] < q['AT2'] and q['A_frac'] < 2 * sig)
    q['BD1_bound'] = float(np.pi ** 2 / 2 * th)
    t = q['T1']; q['Cw2_lower'] = float(1 - (t['S2'] + (2 * np.sqrt(np.sqrt(315.0)) + np.sqrt(15.0) / 3) * np.sqrt(t['A'])) / t['B']
                                        - np.sqrt(2 * np.sqrt(1260.0) * max(th, 0.0) / t['B']))
    return q


def ev(sig, bg, rip, Lv, st=FAST, **kw):
    """one evaluation, stored at once in results_EE/cache (skipped when stored); EE_NOCACHE=1 computes afresh and stores nothing."""
    o = dict(st); o.update(kw); fn = os.path.join(CACHE, _key(sig, bg, rip, Lv, o) + '.json')
    if not NOCACHE and os.path.exists(fn): return json.load(open(fn))
    q = quantities(sig, bg, rip, Lv, **o)
    q = json.loads(json.dumps(q, default=lambda x: bool(x) if isinstance(x, np.bool_) else float(x)))
    if not NOCACHE:
        os.makedirs(CACHE, exist_ok=True); tmp = fn + '.tmp'; json.dump(q, open(tmp, 'w')); os.replace(tmp, fn)
    return q
