"""Stage EI (paper K): the constants that the manuscript adds to those of stage EH (results_EH/eh_constants.json), by the same interval arithmetic (eh_iv.py):
  - Theorem EG-O with J0 = 1.9e-8 (the range sigma <= 1.5e-14 of Theorem 1 (b)); the sigma bound sigma^2 <= 2 pi^2 J/(1 - 2J) written out (review point 5);
  - the five terms of Theorem 1 at sigma = 1.5e-14 (the O-term is the smallest) and the value of the bound there;
  - the closed form at sigma = 1e-15 (as in stage EH);
  - the alternative choice delta = 0.4/L of the review (omega K_b^2 < 1 at L = 4);
  - the monotonicity of sigma k_bar(sigma): psi'/psi < 2.5/L (written as inequalities, review point 6), checked at L = 4 and on a grid.
Output results_EI/ei_constants.json.  Python 3.11 grammar.  Model: Claude Opus 5.5 (claude-opus-5-5)."""
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ed_common as C
from eh_iv import Iv, PI, E, sqrt, exp, ln
from fractions import Fraction as Fr

OUT = []
RES_EI = os.path.join(C.ROOT, 'results_EI')


def chk(key, what, enc, claim, side):
    e = Iv.c(enc); c = Iv.c(claim); ok = {'<=': e <= c, '>=': e >= c, '<': e < c, '>': e > c}[side]
    OUT.append({'key': key, 'what': what, 'lo': str(e.lo), 'hi': str(e.hi), 'claim': str(c.lo), 'side': side, 'ok': bool(ok)})
    print('%-14s %-74s [%.10g, %.10g] %s %s  %s' % (key, what, float(e.lo), float(e.hi), side, str(c.lo), 'ok' if ok else '** FAILS **'))


def psi(L): L = Iv.c(L); return 289 * L ** Iv(Fr(5, 2)) / ln(2 * L) ** Iv(Fr(3, 2)) + 4 * L + Iv('3.57')
def dpsi(L):
    L = Iv.c(L); a = 289 * L ** Iv(Fr(5, 2)) / ln(2 * L) ** Iv(Fr(3, 2))
    return a * (Iv(Fr(5, 2)) / L - Iv(Fr(3, 2)) / (L * ln(2 * L))) + 4


def main():
    one = Iv(1)
    # Theorem EG-O with J0 = 1.9e-8
    J0 = Iv('1.9e-8'); s2max = 2 * PI * PI * J0 / (1 - 2 * J0); al = one / (2 * (PI * PI + s2max))
    chk('O.s2', 'sigma^2 <= 2 pi^2 J0/(1 - 2 J0)', s2max, '3.751e-7', '<=')
    chk('O.alpha', 'alpha >= 1/(2(pi^2 + sigma^2_max))', al, '0.0506605', '>=')
    LJ = ln(one / J0); h = 562 * J0 * LJ ** 3 / ln(2 * LJ) ** 2; Xi = h.sqrt()
    chk('O.LJ', 'L_J = ln(1/J0) >= 4 (h_O increasing there)', LJ, 4, '>=')
    chk('O.h', 'h_O(J0) = 562 J0 L^3/(ln 2L)^2', h, '5.105e-3', '<=')
    chk('O.Xi', 'Xi = sqrt(h_O(J0))', Xi, '0.07145', '<=')
    t1 = h / al; t2 = Xi / (2 * al.sqrt()); t3 = sqrt(15) / 6 * J0.sqrt() / al
    chk('O.t1', 'Xi^2/alpha', t1, '0.1008', '<=')
    chk('O.t2', 'Xi/(2 sqrt alpha)', t2, '0.1588', '<=')
    chk('O.t3', '(sqrt15/6) sqrt(J0)/alpha', t3, '0.00176', '<=')
    chk('O.br', 'the bracket Xi^2/alpha + Xi/(2 sqrt alpha) + (sqrt15/6) sqrt(J0)/alpha', t1 + t2 + t3, '0.2471', '<=')
    chk('O.br4', 'the bracket <= 1/4', t1 + t2 + t3, Fr(1, 4), '<=')
    # Theorem 1 (b) at sigma = 1.5e-14 and (c) at 1e-15
    for tag, sg_s in (('b', '1.5e-14'), ('c', '1e-15')):
        sg = Iv(sg_s); l = ln(one / sg); alp = one / (2 * (PI * PI + sg * sg)); Lb = 2 * l + ln(one / alp); kb = psi(Lb).sqrt()
        J = 112 * sg * sg + 1268 * sg * kb; bnd = Iv(3) / 4 * alp / (alp * sg + 1268 * kb); cc = sqrt(Fr(1, 10)); K = Iv(400)
        bp = 2 * sg / (1 - 2 * sg) + Iv('1.291') * sqrt(2 / ((1 - 2 * sg) * K)) + Iv('8.4258') * sqrt(2 * cc / (K * sqrt(l)))
        others = [one, cc / (2 * l.sqrt()), Iv('0.14') * (1 - 2 * sg), (1 - 2 * sg) * (1 - Iv('0.284')) / (4 * sqrt(2) * kb)]
        chk(tag + '.Lbar', 'L_bar at sigma = %s' % sg_s, Lb, '66.6' if tag == 'b' else '72.06', '>=')
        chk(tag + '.kbar', 'k_bar at sigma = %s' % sg_s, kb, '984.16' if tag == 'b' else '1072.3', '<=')
        chk(tag + '.J', 'A + B inside O <= 112 sigma^2 + 1268 sigma k_bar at sigma = %s' % sg_s, J, '1.872e-8' if tag == 'b' else '1.3598e-9', '<=')
        chk(tag + '.JleJ0', 'that is <= J0 = 1.9e-8', J, J0, '<=')
        chk(tag + '.bp', 'beta\' <= 0.284 at sigma = %s' % sg_s, bp, '0.284', '<=')
        chk(tag + '.min', 'the O-term (3/4) alpha/(alpha sigma + 1268 k_bar) is the smallest of the five', bnd, min(o.lo for o in others), '<')
        chk(tag + '.bound', 'the O-term at sigma = %s' % sg_s, bnd, '3.044e-8' if tag == 'b' else '2.794e-8', '>=')
        chk(tag + '.bk', 'the O-term times k_bar', bnd * kb, '2.996e-5', '>=')
        if tag == 'c':
            chk('c.cf', '1.7e-6 (ln 2 L_bar)^{3/4}/L_bar^{5/4} <= the O-term at sigma = 1e-15', Iv('1.7e-6') * ln(2 * Lb) ** Iv(Fr(3, 4)) / Lb ** Iv(Fr(5, 4)), bnd, '<=')
    # review point 7: delta = 0.4/L also passes the return step at L = 4
    L = Iv(4); tau = ln(2 * L) / (8 * L * L); dl = Iv('0.4') / L; om = (tau + dl) / (2 * (1 - tau))
    chk('d04', 'delta = 0.4/L at L = 4: omega K_b^2 = 4 L omega = (0.8 + ln(2L)/(4L))/(1 - tau) < 1', 4 * L * om, '0.9455', '<=')
    # review point 6: psi'/psi < 2.5/L <= 1 for L >= 4  (so d ln(sigma k_bar)/d ell < -1 + psi'/psi < 0)
    chk('psi.ratio4', 'psi\'(L)/psi(L) * L/2.5 at L = 4 (< 1)', dpsi(4) / psi(4) * Iv(4) / Iv(Fr(5, 2)), 1, '<')
    grid_ok = True
    for k in range(0, 60):
        x = Iv(Fr(4 * 13 ** k, 10 ** k).limit_denominator(10 ** 9)) if k else Iv(4)
        grid_ok = grid_ok and bool(dpsi(x) / psi(x) * x / Iv(Fr(5, 2)) < 1)
    print('grid L = 4 .. 4 * 1.3^59: psi\'/psi < 2.5/L:', grid_ok)
    out = {'checks': OUT, 'n': len(OUT), 'n_ok': sum(c['ok'] for c in OUT), 'grid_psi_ratio_ok': grid_ok, 'env': C.env(), 'arithmetic': 'decimal, 60 digits, outward rounding (eh_iv.py)'}
    out['ok'] = bool(out['n_ok'] == out['n'] and grid_ok); os.makedirs(RES_EI, exist_ok=True); C.save('ei_constants.json', out, RES_EI)
    print('%d/%d hold; ok %s' % (out['n_ok'], out['n'], out['ok']))


if __name__ == '__main__':
    main()
