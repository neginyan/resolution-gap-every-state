"""Stage EH, 1.2: every numerical constant of the chain of Theorem EG, recomputed with the outward-rounded interval arithmetic of eh_iv.py.
Each check is 'claimed value  vs  enclosure of the exact expression'; a check passes when the enclosure lies on the right side of the claimed value.
Where a claim holds 'for all L >= 4' (or all sigma <= ...), the report gives the analytic monotonicity argument; here the claim is checked at the
end point with intervals and, in addition, on a logarithmic grid of L (sampled, not a proof; marked 'grid').
Output results_EH/eh_constants.json.  Python 3.11 grammar."""
import os, sys, json, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ed_common as C
from eh_iv import Iv, PI, E, sqrt, exp, ln
from fractions import Fraction as Fr

OUT = []
PRINTED = []      # rounded values printed in the reports of stages EE-EG and used there as bounds (their direction is checked; a failure is a slip)


def chk(item, what, enclosure, claim, side, note='', printed=False):
    """side '<=' : enclosure <= claim ;  '>=' : enclosure >= claim ; '<' / '>' strict."""
    e = Iv.c(enclosure); c = Iv.c(claim)
    ok = {'<=': e <= c, '>=': e >= c, '<': e < c, '>': e > c}[side]
    (PRINTED if printed else OUT).append({'item': item, 'what': what, 'enclosure': [float(e.lo), float(e.hi)], 'claim': float(c.lo), 'side': side, 'ok': bool(ok), 'note': note})
    print('%-6s %-78s %s %s %-12s %s' % (item, what, repr(e), side, '%.6g' % float(c.lo), 'ok' if ok else ('slip (printed value)' if printed else '** FAILS **')))
    return ok


def psi(L): return 289 * Iv.c(L) ** Iv(Fr(5, 2)) / ln(2 * Iv.c(L)) ** Iv(Fr(3, 2)) + 4 * Iv.c(L) + Iv('3.57')


def main():
    one = Iv(1)
    # ---------------------------------------------------------------- [E3] EE-3 / Cramer-Rao floor
    s4 = Iv('1e-4'); al4 = one / (2 * (PI * PI + s4 * s4))
    chk('E3', 'alpha = 1/(2(pi^2 + sigma^2)) at sigma = 1e-4 (alpha decreases in sigma)', al4, '0.0506', '>=')
    # ---------------------------------------------------------------- [D2] ED half lemma (r = 2, lambda = 1/4)
    chk('D2', 'C_{2,1/4} = M_2^{1/2} M_2^{1/2} (1/2)^{-5/4} = 3 * 2^{5/4}', 3 * Iv(2) ** Iv(Fr(5, 4)), '7.136', '<=')
    chk('D2', 'coefficient r\'/lambda = 2/(1/4)', Iv(8), 8, '<=')
    # ---------------------------------------------------------------- [F1] paper F G2: Lagrange remainder
    chk('F1', '(sqrt960)/6  (E delta^6 <= 2^6 * 15 sigma^6)', sqrt(960) / 6, '5.1640', '<=')
    # ---------------------------------------------------------------- [EF1] numerator
    nx = 2 * sqrt(2) + sqrt(Fr(3, 2)); c1 = sqrt(3) / 2 + sqrt(3) + sqrt(960) / 6; c2 = nx / sqrt(2); c3 = nx / 2
    chk('EF1', '||X_+-||_2 <= 2 sqrt2 + sqrt(3/2)', nx, '4.0533', '<=')
    chk('EF1', 'c1 = sqrt3/2 + sqrt3 + sqrt960/6', c1, '7.77', '<=')
    chk('EF1', 'c2 = (2 sqrt2 + sqrt(3/2))/sqrt2', c2, '2.87', '<=')
    chk('EF1', 'c3 = (2 sqrt2 + sqrt(3/2))/2', c3, '2.03', '<=')
    # ---------------------------------------------------------------- [EF2] threshold 112
    c4 = sqrt(Iv('2.03') ** 2 + 2 * Iv('2.87') ** 2)
    chk('EF2', 'c4 = sqrt(2.03^2 + 2 * 2.87^2)', c4, '4.5382', '<=')
    chk('EF2', '(1/2) 112 - c4 sqrt112 - 7.77  (> 0: the threshold M >= 112 sigma^2)', Iv(56) - c4 * sqrt(112) - Iv('7.77'), 0, '>')
    chk('EF2', 'alternative: 4 (c1 + c2^2 + c3^2)  (A >= 80.5 sigma^2 suffices)', 4 * (c1 + c2 * c2 + c3 * c3), '80.5', '<=')
    # ---------------------------------------------------------------- [EEC] improved EE-C constants
    chk('EEC', 'sqrt15/3  (|E[s_v^3 s_w]|/3 <= sqrt(15 A)/3)', sqrt(15) / 3, '1.291', '<=')
    chk('EEC', 'printed in EE 3.3 (original EE-C, not used by Theorem EG): c1 = 2*315^{1/4} + sqrt15/3 <= 9.7167', 2 * Iv(315) ** Iv(Fr(1, 4)) + sqrt(15) / 3, '9.7167', '<=', printed=True)
    chk('G1', 'printed in EG 2.4 step 7: 1 - s >= 0.9837 at L = 4', 1 - ln(Iv(8)) / 128, '0.9837', '>=', printed=True)
    chk('GI', 'printed in EG 4.2 step 4: (1 - s)^3 >= 0.952 at L = 4', (1 - ln(Iv(8)) / 128) ** 3, '0.952', '>=', printed=True)
    chk('EFO', 'printed in EF 2.4: 1.291 sqrt112 <= 13.67', Iv('1.291') * sqrt(112), '13.67', '<=', printed=True)
    chk('EFO', 'printed in EF 2.4: 8.4258 sqrt224 <= 126.11', Iv('8.4258') * sqrt(224), '126.11', '<=', printed=True)
    chk('EFO', 'printed in EF 2.4: 4 sqrt(224 * 112) <= 633.6', 4 * sqrt(224 * 112), '633.6', '<=', printed=True)
    chk('GO', 'printed in EG 4.3: alpha >= 0.0506', one / (2 * (PI * PI + Iv('1e-6'))), '0.0506', '>=', printed=True)
    chk('GO', 'printed in EG 4.3: Xi^2/alpha <= 0.0132', Iv('6.64e-4') / Iv('0.0506'), '0.0132', '<=', printed=True)
    chk('GO', 'printed in EG 4.3: Xi/(2 sqrt alpha) <= 0.0574', Iv('0.0258') / (2 * sqrt(Iv('0.0506'))), '0.0574', '<=', printed=True)
    chk('EEC', 'c_2 = sqrt(2 sqrt1260)  (E G^4 <= sqrt(105 * 1680) = 420)', sqrt(2 * sqrt(1260)), '8.4258', '<=')
    OUT.append({'item': 'EEC', 'what': '105 * 1680 = 420^2 (integers)', 'enclosure': [176400, 176400], 'claim': 176400, 'side': '==', 'ok': 105 * 1680 == 420 ** 2, 'note': ''}); print('EEC    105 * 1680 = 420^2 (integers):', 105 * 1680 == 420 ** 2)
    # ---------------------------------------------------------------- [EFO] Lemma EF-O
    s3 = Iv('1e-3')
    chk('EFO', '1.291 sqrt112 + 8.4258 sqrt224 + 112 sigma  (sigma <= 1e-3)', Iv('1.291') * sqrt(112) + Iv('8.4258') * sqrt(224) + 112 * s3, 140, '<=')
    chk('EFO', '8 sqrt(224 * 112)  (B <= 2 * 4 k sqrt(theta B A)/sigma^2 ... ; C~_w >= 1/2)', 8 * sqrt(224 * 112), 1268, '<=')
    # ---------------------------------------------------------------- [EFI] Theorem EF (I), K = 400, c^2 = 0.1, sigma <= 1e-4
    K = Iv(400); cc = sqrt(Fr(1, 10)); l4 = ln(Iv(10000))
    bp = 2 * s4 / (1 - 2 * s4) + Iv('1.291') * sqrt(2 / ((1 - 2 * s4) * K)) + Iv('8.4258') * sqrt(2 * cc / (K * sqrt(l4)))
    chk('EFI', 'beta\' at sigma = 1e-4 (decreasing as sigma decreases)', bp, '0.284', '<=')
    chk('EFI', '112/K = 112/400', Iv(112) / K, '0.28', '<=')
    # ---------------------------------------------------------------- [G1] Theorem EG-1, L >= 4 (end point L = 4)
    L = Iv(4); dl = one / (4 * L); T = one + dl; eta = ln(2 * L) / (4 * L); s = ln(2 * L) / (8 * L * L); r1 = dl * (1 - eta); r2 = dl + s; r3 = T / 2
    chk('G1', 'eta = ln(2L)/(4L) at L = 4 (decreasing)', eta, '0.13', '<')
    chk('G1', 's = ln(2L)/(8L^2) at L = 4 (decreasing)', s, '0.0163', '<')
    chk('G1', '-ln(1 - eta) - 2 eta at eta = 0.13 (<= 0: the function is increasing... see report)', -ln(1 - Iv('0.13')) - 2 * Iv('0.13'), 0, '<=')
    chk('G1', 'r2 < r3 at L = 4', r2, r3, '<')
    chk('G1', 'theta = ln(r2/r1)/ln(r3/r1) at L = 4, times L (< 1)', ln(r2 / r1) / ln(r3 / r1) * L, 1, '<')
    al = (s + dl) / (2 * (1 - s))
    chk('G1', 'alpha K^2 = 4 L alpha at L = 4 (< 1)', 4 * L * al, 1, '<')
    chk('G1', 'T/(1-s) at L = 4', T / (1 - s), '1.0801', '<=')
    cb = sqrt(6 * T / (1 - s)); chk('G1', 'sqrt(6 T/(1 - s))  (contribution of E_bad, x A)', cb, '2.55', '<=')
    c3c = 4 * E / sqrt(3); chk('G1', '4e/sqrt3  (three circles, x sqrt(L/ln 2L) A)', c3c, '6.2777', '<=')
    chk('G1', '2e * 4e/sqrt3  (outside E_bad)', 2 * E * c3c, '34.13', '<=')
    r4 = sqrt(L / ln(2 * L)); chk('G1', 'printed in EG 2.4 step 7 as a lower bound: sqrt(4/ln 8) >= 1.387', r4, '1.387', '>=', printed=True)
    chk('G1', 'sqrt(4/ln 8) >= 1.3869 (the correct lower bound)', r4, '1.3869', '>=')
    chk('G1', 'total: 2e (4e/sqrt3) + 2.55/sqrt(4/ln8)', 2 * E * c3c + Iv('2.55') / r4, '36.1', '<=')
    chk('G1', '36.1 * 8  (Q2 <= F(1-s)/s, 1/s = 8L^2/ln 2L)', Iv('36.1') * 8, 289, '<=')
    chk('G1', 'psi(4) = 289 * 32/(ln 8)^{3/2} + 19.57', psi(4), '3100', '>=')
    chk('G1', 'e^4 + 4*4 + 3.57  (A > e^-4: k_eff^2 < e^4 + 4L + 3.57)', E ** 4 + Iv('19.57'), '74.2', '<=')
    # ---------------------------------------------------------------- [GI] Lemma EG-iso, EG-sm
    chk('GI', 'outer, isotropic: T^3/(T^2 - (T/2)^2)^2 = 16/(9T) <= 16/9', T ** 3 / (T * T - T * T / 4) ** 2, Iv(16) / 9, '<=')
    ci = E * (Iv(16) / 9) ** Iv(Fr(1, 4)) * 4; chk('GI', 'three circles, isotropic: e (16/9)^{1/4} * 4', ci, '12.56', '<=')
    bi = sqrt(6) * T * T / (1 - s) ** 3; chk('GI', 'sqrt6 T^2/(1 - s)^3  (E_bad, isotropic)', bi, '2.91', '<=')
    chk('GI', 'printed in EG 4.2 step 4 as a lower bound: 4/ln 8 >= 1.924', L / ln(2 * L), '1.924', '>=', printed=True)
    chk('GI', '4/ln 8 >= 1.9235 (the correct lower bound)', L / ln(2 * L), '1.9235', '>=')
    chk('GI', 'total: 2e * 12.56 + 2.91/(4/ln 8)', 2 * E * Iv('12.56') + Iv('2.91') / (L / ln(2 * L)), '70.2', '<=')
    chk('GI', '70.2 * 8  (E[(d_e^2 p/p)^2] <= G(1-s)/s)', Iv('70.2') * 8, 562, '<=')
    # ---------------------------------------------------------------- [GD] Lemma EG-D
    chk('GD', 'sqrt15/6  (|R\'| <= (sqrt15/6) sigma^2 sqrt(J_pp))', sqrt(15) / 6, '0.6455', '<=')
    # ---------------------------------------------------------------- [GO] Theorem EG-O, J0 = 2e-9
    J0 = Iv('2e-9'); LJ = ln(one / J0); h = 562 * J0 * LJ ** 3 / ln(2 * LJ) ** 2; Xi = h.sqrt()
    alo = one / (2 * (PI * PI + Iv('1e-6')))                    # sigma^2 <= J/alpha <= 4.0e-8 is automatic; alpha >= this
    chk('GO', 'h(2e-9) = 562 J0 L^3/(ln 2L)^2', h, '6.64e-4', '<=')
    chk('GO', 'Xi = sqrt(h(2e-9))', Xi, '0.0258', '<=')
    br = h / alo + Xi / (2 * alo.sqrt()) + sqrt(15) / 6 * J0.sqrt() / alo
    chk('GO', 'bracket Xi^2/alpha + Xi/(2 sqrt alpha) + 0.6455 sqrt(J0)/alpha', br, '0.25', '<=')
    chk('GO', 'h increasing on (0, e^-4]: 1 - 3/L + 2/(L ln 2L) at L = 4', 1 - Iv(3) / 4 + 2 / (Iv(4) * ln(Iv(8))), 0, '>')
    # ---------------------------------------------------------------- [GT] Theorem EG at sigma = 1e-15 (and the monotonicity in sigma)
    rows = []
    for e10 in (15, 20, 30, 50, 100, 300):
        sg = Iv('1e-%d' % e10); l = ln(one / sg); alp = one / (2 * (PI * PI + sg * sg)); Lb = 2 * l + ln(one / alp); kb = psi(Lb).sqrt()
        Jo = 112 * sg * sg + 1268 * sg * kb; bnd = Iv(3) / 4 * alp / (alp * sg + 1268 * kb); cf = Iv('1.7e-6') * ln(2 * Lb) ** Iv(Fr(3, 4)) / Lb ** Iv(Fr(5, 4))
        rows.append({'sigma_exp': -e10, 'L_bar': Lb.f(), 'k_bar': kb.f(), 'J_in_O': Jo.f(), 'bound': bnd.f(), 'closed_form': cf.f()})
        if e10 == 15:
            chk('GT', 'L_bar = 2 ln(1e15) + ln(2(pi^2 + sigma^2))', Lb, 72, '>=')
            chk('GT', 'k_bar = sqrt(psi(L_bar)) at sigma = 1e-15', kb, '1072.4', '<=')
            chk('GT', 'J inside O <= 112 sigma^2 + 1268 sigma k_bar at sigma = 1e-15', Jo, '2e-9', '<=')
            chk('GT', '(3/4) alpha/(alpha sigma + 1268 k_bar) * k_bar', bnd * kb, '2.9e-5', '>=')
            chk('GT', 'psi(L) <= 289.1 L^{5/2}/(ln 2L)^{3/2} at L = 72: (4L + 3.57)(ln 2L)^{3/2}/L^{5/2}', (4 * Iv(72) + Iv('3.57')) * ln(Iv(144)) ** Iv(Fr(3, 2)) / Iv(72) ** Iv(Fr(5, 2)), '0.1', '<=')
            chk('GT', 'sqrt(289.1)', sqrt(Iv('289.1')), '17.01', '<=')
            chk('GT', '0.75 alpha/1268/17.01  (>= 1.7e-6)', Iv(3) / 4 * alp / 1268 / Iv('17.01'), '1.7e-6', '>=')
        chk('GT', 'bound >= 1.7e-6 (ln 2L)^{3/4}/L^{5/4} at sigma = 1e-%d' % e10, bnd, cf, '>=')
        chk('GT', 'J inside O <= 2e-9 at sigma = 1e-%d' % e10, Jo, '2e-9', '<=')
        other = [one, Iv(1) / 2 * cc / l.sqrt(), Iv('0.14') * (1 - 2 * sg), (1 - 2 * sg) * (1 - Iv('0.284')) / (4 * sqrt(2) * kb)]
        chk('GT', 'the O-term is the smallest of the five at sigma = 1e-%d' % e10, bnd, min(o.lo for o in other), '<')
    # ---------------------------------------------------------------- [EXT] section 1.4 (a separate, new claim; report EH, 5): J0 = 1.9e-8 and sigma_0 = 1.5e-14
    J1 = Iv('1.9e-8'); LJ1 = ln(one / J1); h1 = 562 * J1 * LJ1 ** 3 / ln(2 * LJ1) ** 2
    br1 = h1 / alo + h1.sqrt() / (2 * alo.sqrt()) + sqrt(15) / 6 * J1.sqrt() / alo
    chk('EXT', 'new (section 1.4): bracket of Theorem EG-O at J0 = 1.9e-8', br1, '0.25', '<=', printed=False)
    sg = Iv('1.5e-14'); l = ln(one / sg); alp = one / (2 * (PI * PI + sg * sg)); Lb = 2 * l + ln(one / alp); kb = psi(Lb).sqrt()
    chk('EXT', 'new (section 1.4): J inside O <= 112 sigma^2 + 1268 sigma k_bar at sigma = 1.5e-14', 112 * sg * sg + 1268 * sg * kb, '1.9e-8', '<=')
    chk('EXT', 'new (section 1.4): alpha >= 0.050660 for sigma^2 <= J0/alpha (sigma^2 <= 4e-7)', one / (2 * (PI * PI + Iv('4e-7'))), '0.05066', '>=')
    # ---------------------------------------------------------------- grids (sampled; the report gives the analytic arguments)
    grid = []
    Ls = [4 * 1.25 ** k for k in range(0, 70)]
    for x in Ls:
        Lx = Iv(Fr(x).limit_denominator(10 ** 6)); dl = one / (4 * Lx); T = one + dl; eta = ln(2 * Lx) / (4 * Lx); s = ln(2 * Lx) / (8 * Lx * Lx)
        r1 = dl * (1 - eta); r2 = dl + s; r3 = T / 2; th = ln(r2 / r1) / ln(r3 / r1); al = (s + dl) / (2 * (1 - s))
        grid.append({'L': float(Lx.lo), 'theta_L': (th * Lx).f(), 'alpha_4L': (4 * Lx * al).f(), 'psi_increasing': bool(psi(Lx * Iv('1.0001')) > psi(Lx))})
    gok = all(g['theta_L'][1] < 1 and g['alpha_4L'][1] < 1 and g['psi_increasing'] for g in grid)
    print('grid (L = 4 .. %.3g, %d points): theta L < 1, 4 L alpha < 1, psi increasing: %s' % (Ls[-1], len(Ls), gok))
    sgrid = []
    for e10 in range(15, 301, 5):
        sg = Iv('1e-%d' % e10); l = ln(one / sg); alp = one / (2 * (PI * PI + sg * sg)); Lb = 2 * l + ln(one / alp); sgrid.append((e10, (sg * psi(Lb).sqrt()).f()))
    sok = all(sgrid[i + 1][1][1] < sgrid[i][1][0] for i in range(len(sgrid) - 1))
    print('grid sigma = 1e-15 .. 1e-300: sigma k_bar(sigma) decreasing as sigma decreases: %s' % sok)
    out = {'checks': OUT, 'printed': PRINTED, 'printed_slips': [p['what'] for p in PRINTED if not p['ok']], 'n': len(OUT), 'n_ok': sum(c['ok'] for c in OUT), 'theorem_rows': rows, 'grid_L': grid, 'grid_L_ok': gok, 'grid_sigma_ok': sok, 'env': C.env(),
           'arithmetic': 'decimal, 60 digits, outward rounding (eh_iv.py)'}
    out['ok'] = bool(out['n_ok'] == out['n'] and gok and sok); C.save('eh_constants.json', out, C.RES_EH)
    print('%d/%d constant checks hold; printed values: %d checked, slips: %s; ok %s' % (out['n_ok'], out['n'], len(PRINTED), out['printed_slips'], out['ok']))


if __name__ == '__main__':
    main()
