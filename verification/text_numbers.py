r"""Paper K -- the numbers printed in the text, checked against the stored records and the interval arithmetic.  Python 3.11 grammar.
    python verification/text_numbers.py [--paper DIR]        (DIR holds sections/*.tex; default: paper/ of this repository)
Two passes.
 (1) Every entry of ENTRIES names a file of the paper, a piece of its text (which must occur there verbatim), the number printed in it, how the
     printed number relates to the true one, and where the true one comes from:
       ub     the text uses the printed number as an upper bound: the upper end of the true value (interval) must not exceed it;
       lb     a lower bound: the lower end of the true value must not be below it;
       lt     strict upper bound;
       round  the printed number is the true value rounded to the printed digits (|true - printed| <= half a unit of the last printed digit);
       trunc  the printed number ends with \ldots: the true value starts with the printed digits;
       oom    an order of magnitude (|log10 true - log10 printed| <= 1/2);
       eq     exact (integers, counts, parameters of a recorded state, exact rationals);
       cond   the number is a hypothesis or a range end, and the check is the condition that justifies it (a list of interval checks, or a
              property of a record).
     True values come from results_EH/eh_constants.json and results_EI/ei_constants.json (outward-rounded intervals, 60 digits), from the records
     of the recorded states (results_EF, results_EG, results_EH, results_EI) and, for one entry, from verification/common.py.
 (2) An inventory: every number in the running text (decimal numbers, a x 10^k, 10^k and integers >= 10; the generated tables C_table, C_printed,
     D_table and the references are not part of it) must be covered by an entry of (1) or by an entry of ALLOW (a definition, a label of a
     scenario, a layout width, an exact constant inside a formula), with its reason.
Output results_EI/ei_text_numbers.json (ok = every entry holds, every context is found, every number is covered)."""
import json, math, os, re, sys
from decimal import Decimal as D, getcontext

getcontext().prec = 60
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PAPER = os.path.join(ROOT, 'paper')
if '--paper' in sys.argv: PAPER = os.path.abspath(sys.argv[sys.argv.index('--paper') + 1])
SKIP = ('C_table.tex', 'C_printed.tex', 'D_table.tex', '11_references.tex')


def rec(path): return json.load(open(os.path.join(ROOT, path), encoding='utf-8'))


EHC = rec('results_EH/eh_constants.json'); EIC = rec('results_EI/ei_constants.json'); RNG = rec('results_EH/eh_range.json')
JO = rec('results_EH/eh_joints.json'); EGC = {p: rec('results_EG/eg_circ_%s.json' % p) for p in ('2d', 'iso', '1d', '1ds')}
EGP = rec('results_EG/eg_point.json'); EGD = rec('results_EG/eg_direct.json'); EF2 = rec('results_EF/ef_q2_2d.json'); EIP = rec('results_EI/ei_point.json')


def ehc(what):
    """interval (lo, hi) of the check of eh_constants.json with this description (exactly one)."""
    m = [c for c in EHC['checks'] if c['what'] == what]
    assert len(m) == 1, what
    return D(repr(m[0]['enclosure'][0])), D(repr(m[0]['enclosure'][1]))


def ehc_ok(what): return [c for c in EHC['checks'] if c['what'] == what][0]['ok']


def eic(key):
    m = [c for c in EIC['checks'] if c['key'] == key]
    assert len(m) == 1, key
    return D(m[0]['lo']), D(m[0]['hi'])


def eic_ok(*keys): return all([c for c in EIC['checks'] if c['key'] == k][0]['ok'] for k in keys)


def joint(J): return [j for j in JO['joints'] if j['J'] == J][0]


def rng(prefix): return D(repr([r for r in RNG['rows'] if r['scenario'].startswith(prefix)][0]['sigma_0']))


def pt(x): x = D(repr(x)) if not isinstance(x, D) else x; return (x, x)


def egc_min(key): return pt(min(EGC[p]['all_min_margins'][key] for p in EGC))


def common_rtol():
    sys.path.insert(0, os.path.join(ROOT, 'verification'))
    import common
    return pt(common.FLOAT_RTOL)


C224 = ehc('8 sqrt(224 * 112)  (B <= 2 * 4 k sqrt(theta B A)/sigma^2 ... ; C~_w >= 1/2)')
C21 = ehc('C_{2,1/4} = M_2^{1/2} M_2^{1/2} (1/2)^{-5/4} = 3 * 2^{5/4}')
THR = '(1/2) 112 - c4 sqrt112 - 7.77  (> 0: the threshold M >= 112 sigma^2)'
B_ALL = ('b.Lbar', 'b.kbar', 'b.J', 'b.JleJ0', 'b.bp', 'b.min', 'b.bound', 'b.bk')
C_ALL = ('c.Lbar', 'c.kbar', 'c.J', 'c.JleJ0', 'c.bp', 'c.min', 'c.bound', 'c.bk', 'c.cf')
N_STATES = JO['n_engine'] + EF2['n_states'] + EGC['1d']['n_states']
E4 = D(4).exp()
TV = lambda: (D(1) - 2 * D('1.5e-14')) * (D(1) - D('0.284')) / (4 * D(2).sqrt())      # the fourth term of (a) times k_bar, at sigma = 1.5e-14

# (file, context, printed, kind, true value or condition, what)
ENTRIES = [
    # abstract
    ('00_abstract', r'every $\sigma\le1.5\times10^{-14}$', r'1.5\times10^{-14}', 'cond', lambda: eic_ok(*B_ALL), 'Theorem 1(b) at its range end (interval)'),
    ('00_abstract', r'\alpha\sigma+1268\,\kb}', '1268', 'ub', lambda: C224, 'C_O >= 8 sqrt(224 112)'),
    ('00_abstract', r'\kb^2=289\,\frac', '289', 'ub', lambda: ehc('36.1 * 8  (Q2 <= F(1-s)/s, 1/s = 8L^2/ln 2L)'), 'C_Q >= 36.1 * 8'),
    ('00_abstract', r'4\Lb+3.57,', '3.57', 'ub', lambda: (C21[0] / 2, C21[1] / 2), '3.57 >= (3 * 2^{5/4})/2'),
    ('00_abstract', r'for $\sigma\le10^{-15}$ this gives', r'10^{-15}', 'cond', lambda: eic_ok(*C_ALL), 'Theorem 1(c) at sigma = 1e-15 (interval)'),
    ('00_abstract', r'$g\ge1.7\times10^{-6}\,(\ln2\Lb)', r'1.7\times10^{-6}', 'lb', lambda: ehc('0.75 alpha/1268/17.01  (>= 1.7e-6)'), '0.75 alpha/(1268 * 17.01)'),
    ('00_abstract', r'on $6062$ recorded states', '6062', 'eq', lambda: pt(N_STATES), '692 + 2449 + 2921 states of the records'),
    # introduction
    ('01_introduction', r'It holds for $\sigma\le1.5\times10^{-14}$', r'1.5\times10^{-14}', 'cond', lambda: eic_ok(*B_ALL), 'range end'),
    ('01_introduction', r'\psi(L)=289L^{5/2}/(\ln2L)^{3/2}+4L+3.57', '289', 'ub', lambda: ehc('36.1 * 8  (Q2 <= F(1-s)/s, 1/s = 8L^2/ln 2L)'), 'C_Q'),
    ('01_introduction', r'\psi(L)=289L^{5/2}/(\ln2L)^{3/2}+4L+3.57', '3.57', 'ub', lambda: (C21[0] / 2, C21[1] / 2), '3.57'),
    ('01_introduction', r'together exceed $112\sigma^2$', '112', 'cond', lambda: ehc_ok(THR) and ehc(THR)[0] > 0, 'the threshold: 56 - c4 sqrt112 - 7.77 > 0'),
    ('01_introduction', r'on $6062$ recorded states', '6062', 'eq', lambda: pt(N_STATES), 'number of recorded states'),
    ('01_introduction', r'fails by a factor of $9\times10^{71}$', r'9\times10^{71}', 'round', lambda: pt(-joint('J18')['min_margin']), 'J18: 1 - wrong/true'),
    ('01_introduction', r'remainder $(\sqrt{960}/6)\sigma^2$', '960', 'eq', lambda: pt(D(2) ** 6 * 15), '2^6 * 15'),
    # setting and Theorem 1
    ('02_main_theorem', r'\Oset=\{M<112\,\sigma^2\}', '112', 'cond', lambda: ehc_ok(THR) and ehc(THR)[0] > 0, 'threshold'),
    ('02_main_theorem', r'\psi(L)=289\,\frac{L^{5/2}}{(\ln2L)^{3/2}}+4L+3.57', '289', 'ub', lambda: ehc('36.1 * 8  (Q2 <= F(1-s)/s, 1/s = 8L^2/ln 2L)'), 'C_Q'),
    ('02_main_theorem', r'\psi(L)=289\,\frac{L^{5/2}}{(\ln2L)^{3/2}}+4L+3.57', '3.57', 'ub', lambda: (C21[0] / 2, C21[1] / 2), '3.57'),
    ('02_main_theorem', r'If $0<\sigma\le10^{-4}$ and $M\ge112\sigma^2$', r'10^{-4}', 'cond',
     lambda: ehc_ok("beta' at sigma = 1e-4 (decreasing as sigma decreases)") and ehc_ok('alpha = 1/(2(pi^2 + sigma^2)) at sigma = 1e-4 (alpha decreases in sigma)'), "beta' and alpha at 1e-4"),
    ('02_main_theorem', r'If $0<\sigma\le10^{-4}$ and $M\ge112\sigma^2$', '112', 'cond', lambda: ehc(THR)[0] > 0, 'threshold'),
    ('02_main_theorem', r'0.14(1-2\sigma)', '0.14', 'eq', lambda: (ehc('112/K = 112/400')[0] / 2, ehc('112/K = 112/400')[1] / 2), 'half of 112/400'),
    ('02_main_theorem', r'(1-2\sigma)(1-0.284)', '0.284', 'ub', lambda: ehc("beta' at sigma = 1e-4 (decreasing as sigma decreases)"), "beta' at 1e-4"),
    ('02_main_theorem', r'If $0<\sigma\le1.5\times10^{-14}$, then for every state', r'1.5\times10^{-14}', 'cond', lambda: eic_ok(*B_ALL), 'range end'),
    ('02_main_theorem', r'\frac{\alpha}{\alpha\sigma+1268\,\kb}', '1268', 'ub', lambda: C224, 'C_O'),
    ('02_main_theorem', r'If $0<\sigma\le10^{-15}$, then for every state', r'10^{-15}', 'cond', lambda: eic_ok(*C_ALL), '(c)'),
    ('02_main_theorem', r'\frac{2.9\times10^{-5}}{\kb}', r'2.9\times10^{-5}', 'lb', lambda: eic('c.bk'), 'O-term times k_bar at 1e-15'),
    ('02_main_theorem', r'\ge\ 1.7\times10^{-6}\,\frac', r'1.7\times10^{-6}', 'lb', lambda: ehc('0.75 alpha/1268/17.01  (>= 1.7e-6)'), '(c) closed form'),
    ('02_main_theorem', r'one has $\Lb=72.06$', '72.06', 'round', lambda: eic('c.Lbar'), 'L_bar at 1e-15'),
    ('02_main_theorem', r'and $\kb=1072.3$', '1072.3', 'round', lambda: eic('c.kbar'), 'k_bar at 1e-15'),
    ('02_main_theorem', r'the bound of (b) is $2.79\times10^{-8}$', r'2.79\times10^{-8}', 'round', lambda: eic('c.bound'), 'bound at 1e-15'),
    ('02_main_theorem', r'$\kb\sim17\,\Lb^{5/4}', '17', 'eq', lambda: pt(D(289).sqrt()), 'sqrt 289'),
    ('02_main_theorem', r'the range $\sigma\le1.5\times10^{-14}$ is narrow', r'1.5\times10^{-14}', 'cond', lambda: eic_ok(*B_ALL), 'range end'),
    ('02_main_theorem', r'For $10^{-14}\lesssim\sigma\le10^{-4}$', r'10^{-14}', 'oom', lambda: pt(D('1.5e-14')), 'the range end, in order of magnitude'),
    # three circles
    ('04_three_circles', r'$\mathrm{SV}\le8AL+7.136A$', '7.136', 'ub', lambda: C21, '3 * 2^{5/4}'),
    ('04_three_circles', r'(ii) $F(1-\tau)\le36.1\sqrt{L/\ln2L}\;A$', '36.1', 'ub', lambda: ehc('total: 2e (4e/sqrt3) + 2.55/sqrt(4/ln8)'), 'total coefficient'),
    ('04_three_circles', r'(iii) $Q_2\le289\,A\,L^{5/2}', '289', 'ub', lambda: ehc('36.1 * 8  (Q2 <= F(1-s)/s, 1/s = 8L^2/ln 2L)'), '36.1 * 8'),
    ('04_three_circles', r'$\keff^2\le74.2\le\psi(4)$', '74.2', 'ub', lambda: ehc('e^4 + 4*4 + 3.57  (A > e^-4: k_eff^2 < e^4 + 4L + 3.57)'), 'e^4 + 19.57'),
    ('04_three_circles', r'$\eta\le0.12997$', '0.12997', 'ub', lambda: ehc('eta = ln(2L)/(4L) at L = 4 (decreasing)'), 'eta at L = 4'),
    ('04_three_circles', r'$\tau\le0.016246$', '0.016246', 'ub', lambda: ehc('s = ln(2L)/(8L^2) at L = 4 (decreasing)'), 'tau at L = 4'),
    ('04_three_circles', r'(for $\eta\le0.13$, $-\ln(1-\eta)\le2\eta$)', '0.13', 'cond',
     lambda: ehc('eta = ln(2L)/(4L) at L = 4 (decreasing)')[1] <= D('0.13') and ehc('-ln(1 - eta) - 2 eta at eta = 0.13 (<= 0: the function is increasing... see report)')[1] <= 0, 'eta <= 0.13 and the inequality there'),
    ('04_three_circles', r'A\le2.55A$', '2.55', 'ub', lambda: ehc('sqrt(6 T/(1 - s))  (contribution of E_bad, x A)'), 'E_bad'),
    ('04_three_circles', r'\le34.13\sqrt{L/\ln2L}\,A$', '34.13', 'ub', lambda: ehc('2e * 4e/sqrt3  (outside E_bad)'), 'outside E_bad'),
    ('04_three_circles', r'is $1.38693\ldots$ at $L=4$', r'1.38693\ldots', 'trunc', lambda: ehc('sqrt(4/ln 8) >= 1.3869 (the correct lower bound)'), 'sqrt(4/ln 8)'),
    ('04_three_circles', r'the sum is at most $36.1\sqrt{L/\ln2L}\,A$', '36.1', 'ub', lambda: ehc('total: 2e (4e/sqrt3) + 2.55/sqrt(4/ln8)'), 'total'),
    ('04_three_circles', r'and $36.1\cdot8\le289$', '289', 'ub', lambda: ehc('36.1 * 8  (Q2 <= F(1-s)/s, 1/s = 8L^2/ln 2L)'), '36.1 * 8'),
    ('04_three_circles', r'$\keff^2=K/A\le Q_2/A+4L+3.568$', '3.568', 'ub', lambda: (C21[0] / 2, C21[1] / 2), '(3 * 2^{5/4})/2'),
    ('04_three_circles', r'$\frac12\mathrm{SV}/A<19.57$ for $L<4$', '19.57', 'ub', lambda: ((32 + C21[0]) / 2, (32 + C21[1]) / 2), '(8 * 4 + 7.136)/2'),
    ('04_three_circles', r'$e^4+19.57<74.17\le\psi(4)=3103.7$', '74.17', 'lt', lambda: pt(E4 + D('19.57')), 'e^4 + 19.57'),
    ('04_three_circles', r'$e^4+19.57<74.17\le\psi(4)=3103.7$', '3103.7', 'round', lambda: ehc('psi(4) = 289 * 32/(ln 8)^{3/2} + 19.57'), 'psi(4)'),
    ('04_three_circles', r'for instance $\delta=0.4/L$ gives', '0.4', 'cond', lambda: eic_ok('d04'), 'the alternative delta (interval)'),
    ('04_three_circles', r'$\omega K_b^2=(0.8+\ln(2L)/(4L))/(1-\tau)=0.9453$', '0.8', 'eq', lambda: pt(D(2) * D('0.4')), '4L * 0.4/(4L) * 2'),
    ('04_three_circles', r'$\omega K_b^2=(0.8+\ln(2L)/(4L))/(1-\tau)=0.9453$', '0.9453', 'round', lambda: eic('d04'), 'omega K_b^2 at L = 4'),
    ('04_three_circles', r'$1441$ two-dimensional states', '1441', 'eq', lambda: pt(EGC['iso']['n_states']), 'eg_circ: complex-time states'),
    ('04_three_circles', r'($889$ of them with $L\ge4$)', '889', 'eq', lambda: pt(EGC['iso']['n_L_ge_4']), 'of them with L >= 4'),
    ('04_three_circles', r'(ii)--(iv) on $2449$ two-dimensional states', '2449', 'eq', lambda: pt(EF2['n_states']), 'two-dimensional states'),
    ('04_three_circles', r'and $2921$ one-variable states', '2921', 'eq', lambda: pt(EGC['1d']['n_states']), 'one-variable states'),
    ('04_three_circles', r'are $0.702$ (inner circle)', '0.702', 'lb', lambda: egc_min('inner_margin'), 'smallest inner-circle margin'),
    ('04_three_circles', r'and $1.498$ (three circles)', '1.498', 'lb', lambda: egc_min('three_margin'), 'smallest three-circle margin'),
    ('04_three_circles', r'and $\vartheta L\le0.75$', '0.75', 'ub', lambda: pt(max(EGC[p]['theta_L_max'] for p in EGC)), 'largest theta L'),
    ('04_three_circles', r'$F(1-\tau)/A$ is at most $1.06$', '1.06', 'ub',
     lambda: pt(max(f['F_over_A_max'] for p in ('2d', '1d', '1ds') for f in EGC[p]['families'].values() if f.get('F_over_A_max') is not None)), 'largest F(1 - tau)/A'),
    # outside O
    ('05_outside', r'\frac{\sqrt{960}}6\le7.77$', '7.77', 'ub', lambda: ehc('c1 = sqrt3/2 + sqrt3 + sqrt960/6'), 'c1'),
    ('05_outside', r'\frac{\sqrt{960}}6\le7.77$', '960', 'eq', lambda: pt(D(2) ** 6 * 15), '2^6 * 15'),
    ('05_outside', r'/\sqrt2\le2.87$', '2.87', 'ub', lambda: ehc('c2 = (2 sqrt2 + sqrt(3/2))/sqrt2'), 'c2'),
    ('05_outside', r'/2\le2.03$', '2.03', 'ub', lambda: ehc('c3 = (2 sqrt2 + sqrt(3/2))/2'), 'c3'),
    ('05_outside', r'$c_4=\sqrt{2.03^2+2\cdot2.87^2}\le4.5382$', '4.5382', 'ub', lambda: ehc('c4 = sqrt(2.03^2 + 2 * 2.87^2)'), 'c4'),
    ('05_outside', r'equals $0.2028\ldots>0$ at $x=\sqrt{112}$', r'0.2028\ldots', 'trunc', lambda: ehc(THR), 'the threshold margin'),
    ('05_outside', r'equals $0.2028\ldots>0$ at $x=\sqrt{112}$', '112', 'cond', lambda: ehc(THR)[0] > 0, 'threshold'),
    ('05_outside', r'$\tilde C_w\ge1-(A+1.291\sqrt A)/B-8.4258\sqrt{\theta/B}$', '1.291', 'ub', lambda: ehc('sqrt15/3  (|E[s_v^3 s_w]|/3 <= sqrt(15 A)/3)'), 'sqrt15/3'),
    ('05_outside', r'$\tilde C_w\ge1-(A+1.291\sqrt A)/B-8.4258\sqrt{\theta/B}$', '8.4258', 'ub', lambda: ehc('c_2 = sqrt(2 sqrt1260)  (E G^4 <= sqrt(105 * 1680) = 420)'), 'sqrt(2 sqrt1260)'),
    ('05_outside', r'Let $\sigma\le10^{-3}$. If $M<112\sigma^2$', r'10^{-3}', 'cond', lambda: ehc_ok('1.291 sqrt112 + 8.4258 sqrt224 + 112 sigma  (sigma <= 1e-3)'), '140 at sigma = 1e-3'),
    ('05_outside', r'$B\le\sigma\max(280,\,1268\,\keff)$', '1268', 'ub', lambda: C224, 'C_O'),
    ('05_outside', r'$B\le\sigma\max(280,\,1268\,\keff)$', '280', 'eq', lambda: pt(D(2) * 140), '2 * 140'),
    ('05_outside', r'and $\theta B<224\sigma^2$', '224', 'eq', lambda: pt(D(2) * 112), '2 * 112'),
    ('05_outside', r'\ge1-140\sigma/B$', '140', 'ub', lambda: ehc('1.291 sqrt112 + 8.4258 sqrt224 + 112 sigma  (sigma <= 1e-3)'), '1.291 sqrt112 + 8.4258 sqrt224 + 112 sigma'),
    ('05_outside', r'$B\le8\sqrt{224\cdot112}\,\keff\sigma\le1267.2\,\keff\sigma$', '1267.2', 'ub', lambda: C224, '8 sqrt(224 * 112)'),
    ('05_outside', r'Let \textup{(M1)}, \textup{(H2)}, $\sigma\le10^{-4}$', r'10^{-4}', 'cond', lambda: ehc_ok("beta' at sigma = 1e-4 (decreasing as sigma decreases)"), "beta' at 1e-4"),
    ('05_outside', r'$K_r=400$', '400', 'cond', lambda: ehc_ok('112/K = 112/400') and ehc_ok("beta' at sigma = 1e-4 (decreasing as sigma decreases)"), 'K_r = 400 is a choice; it enters 0.28 and beta\''),
    ('05_outside', r'$M/(\sigma(A+B))>112(1-2\sigma)/K_r=0.28(1-2\sigma)$', '0.28', 'eq', lambda: ehc('112/K = 112/400'), '112/400'),
    ('05_outside', r'\le0.284$ for $\sigma\le10^{-4}$', '0.284', 'ub', lambda: ehc("beta' at sigma = 1e-4 (decreasing as sigma decreases)"), "beta'"),
    ('05_outside', r'On $492$ recorded engine states outside', '492', 'eq', lambda: pt(joint('J3')['n']), 'J3 states'),
    ('05_outside', r'a margin of at least $0.99998$ in units of $M/2$', '0.99998', 'lb', lambda: pt(joint('J3')['min_margin']), 'J3'),
    ('05_outside', r'on the $481$ of them', '481', 'eq', lambda: pt(joint('J22')['n']), 'J22 states'),
    ('05_outside', r'by a factor of at least $1809$', '1809', 'lb', lambda: pt(1 + joint('J22')['min_margin']), 'J22: g/bound'),
    # inside O
    ('06_inside', r'its first term $7.77\sigma^2$', '7.77', 'ub', lambda: ehc('c1 = sqrt3/2 + sqrt3 + sqrt960/6'), 'c1'),
    ('06_inside', r'\le0.6455\,\sigma^2\sqrt{J} .', '0.6455', 'ub', lambda: ehc("sqrt15/6  (|R'| <= (sqrt15/6) sigma^2 sqrt(J_pp))"), 'sqrt15/6'),
    ('06_inside', r'$I_e(1-\tau)\le70.2\,(L_e/\ln2L_e)\,I_e(1)$', '70.2', 'ub', lambda: ehc('total: 2e * 12.56 + 2.91/(4/ln 8)'), 'isotropic total'),
    ('06_inside', r'\le12.56(L_e/\ln2L_e)I_e(1)$', '12.56', 'ub', lambda: ehc('three circles, isotropic: e (16/9)^{1/4} * 4'), 'e (16/9)^{1/4} 4'),
    ('06_inside', r'(1-\tau)^{-3}I_e(1)\le2.91\,I_e(1)$', '2.91', 'ub', lambda: ehc('sqrt6 T^2/(1 - s)^3  (E_bad, isotropic)'), 'E_bad isotropic'),
    ('06_inside', r'$L_e/\ln2L_e\ge4/\ln8=1.92359\ldots$', r'1.92359\ldots', 'trunc', lambda: ehc('4/ln 8 >= 1.9235 (the correct lower bound)'), '4/ln 8'),
    ('06_inside', r'the total coefficient is at most $69.80\le70.2$', '69.80', 'ub', lambda: ehc('total: 2e * 12.56 + 2.91/(4/ln 8)'), 'isotropic total'),
    ('06_inside', r'the total coefficient is at most $69.80\le70.2$', '70.2', 'ub', lambda: ehc('total: 2e * 12.56 + 2.91/(4/ln 8)'), 'isotropic total'),
    ('06_inside', r'\hO(x)=562\,x\,(\ln\frac1x)^3', '562', 'ub', lambda: ehc('70.2 * 8  (E[(d_e^2 p/p)^2] <= G(1-s)/s)'), '70.2 * 8'),
    ('06_inside', r'$70.2\cdot8\le562$', '70.2', 'ub', lambda: ehc('total: 2e * 12.56 + 2.91/(4/ln 8)'), 'isotropic total'),
    ('06_inside', r'$J=A+B\le J_0=1.9\times10^{-8}$', r'1.9\times10^{-8}', 'cond', lambda: eic_ok('O.s2', 'O.alpha', 'O.LJ', 'O.h', 'O.Xi', 'O.t1', 'O.t2', 'O.t3', 'O.br', 'O.br4'), 'Theorem 19 at J0 (interval)'),
    ('06_inside', r'\le3.7505\times10^{-7}$', r'3.7505\times10^{-7}', 'ub', lambda: eic('O.s2'), 'sigma^2 bound'),
    ('06_inside', r'and $\alpha\ge0.0506605$', '0.0506605', 'lb', lambda: eic('O.alpha'), 'alpha'),
    ('06_inside', r'\le\hO(J_0)\le4.706\times10^{-3}$', r'4.706\times10^{-3}', 'ub', lambda: eic('O.h'), 'h_O(J0)'),
    ('06_inside', r'and $\Xi\le0.06860$', '0.06860', 'ub', lambda: eic('O.Xi'), 'Xi'),
    ('06_inside', r"and $|R'|\le0.6455\,\sigma^2\sqrt J", '0.6455', 'ub', lambda: ehc("sqrt15/6  (|R'| <= (sqrt15/6) sigma^2 sqrt(J_pp))"), 'sqrt15/6'),
    ('06_inside', r'The bracket is $0.09288+', '0.09288', 'ub', lambda: eic('O.t1'), 'Xi^2/alpha'),
    ('06_inside', r'0.09288+0.15238+0.00176', '0.15238', 'ub', lambda: eic('O.t2'), 'Xi/(2 sqrt alpha)'),
    ('06_inside', r'0.09288+0.15238+0.00176', '0.00176', 'ub', lambda: eic('O.t3'), '0.6455 sqrt(J0)/alpha'),
    ('06_inside', r'0.15238+0.00176\le0.24702\le\frac14', '0.24702', 'ub', lambda: eic('O.br'), 'the bracket'),
    ('06_inside', r'0.15238+0.00176\le0.24702\le\frac14', '0.24702', 'eq', lambda: pt(D('0.09288') + D('0.15238') + D('0.00176')), 'the sum of the three printed terms'),
    ('06_inside', r'three-point mixture ($\sigma=1$ to $0.1$)', '0.1', 'eq', lambda: pt(min(x['sigma'] for x in EGP['d'])), 'smallest sigma of the check'),
    ('06_inside', r'($\sigma=10^{-2}$, $10^{-4}$)', r'10^{-2}', 'eq', lambda: pt(EGP['f'][0]['sigma']), 'blobs'),
    ('06_inside', r'($\sigma=10^{-2}$, $10^{-4}$)', r'10^{-4}', 'eq', lambda: pt(EGP['f'][1]['sigma']), 'blobs'),
    ('06_inside', r'stays below $3\%$ of its bound', '3', 'ub', lambda: pt(100 * max(1 - joint('J24')['min_margin'], 1 - joint('J25')['min_margin'])), "largest |R'|/bound in per cent"),
    ('06_inside', r'with points $(0,0)$, $(1.3,0.45)$, $(-0.7,1.2)$', '1.3', 'eq', lambda: pt(EIP['state']['points_zp'][1][0]), 'state'),
    ('06_inside', r'with points $(0,0)$, $(1.3,0.45)$, $(-0.7,1.2)$', '0.45', 'eq', lambda: pt(EIP['state']['points_zp'][1][1]), 'state'),
    ('06_inside', r'with points $(0,0)$, $(1.3,0.45)$, $(-0.7,1.2)$', '0.7', 'eq', lambda: pt(-EIP['state']['points_zp'][2][0]), 'state'),
    ('06_inside', r'with points $(0,0)$, $(1.3,0.45)$, $(-0.7,1.2)$', '1.2', 'eq', lambda: pt(EIP['state']['points_zp'][2][1]), 'state'),
    ('06_inside', r'weights $0.2$, $0.5$, $0.3$', '0.2', 'eq', lambda: pt(EIP['state']['weights'][0]), 'state'),
    ('06_inside', r'weights $0.2$, $0.5$, $0.3$', '0.5', 'eq', lambda: pt(EIP['state']['weights'][1]), 'state'),
    ('06_inside', r'weights $0.2$, $0.5$, $0.3$', '0.3', 'eq', lambda: pt(EIP['state']['weights'][2]), 'state'),
    ('06_inside', r'and $\sigma=0.55$, and found', '0.55', 'eq', lambda: pt(EIP['state']['sigma']), 'state'),
    ('06_inside', r"and found $|R'|=0.00988$", '0.00988', 'round', lambda: pt(abs(EIP['grid_b']['R_prime'])), "our R' (the reviewer's value agrees)"),
    ('06_inside', r'against the bound $0.6455\,\sigma^2\sqrt J=0.2286$', '0.2286', 'round', lambda: pt(EIP['grid_b']['bound_0.6455_sqrtJ']), 'our bound'),
    ('06_inside', r'a quadrature converged to $10^{-12}$', r'10^{-12}', 'ub', lambda: pt(EIP['max_grid_difference']), 'difference of two grids'),
    ('06_inside', r'gives $2A+N=1.2496652$', '1.2496652', 'round', lambda: pt(EIP['grid_b']['2A+N']), 'our 2A + N'),
    ('06_inside', r"$R'=0.00988$ and the bound $0.2286$", '0.00988', 'round', lambda: pt(abs(EIP['grid_b']['R_prime'])), "our R'"),
    ('06_inside', r"$R'=0.00988$ and the bound $0.2286$", '0.2286', 'round', lambda: pt(EIP['grid_b']['bound_0.6455_sqrtJ']), 'our bound'),
    ('06_inside', r'On the $36$ recorded states', '36', 'eq', lambda: pt(joint('J29')['n']), 'J29 states'),
    ('06_inside', r'with $J\le2\times10^{-9}$ and (H2)', r'2\times10^{-9}', 'cond', lambda: JO['n_J_le_2e9_H2'] == joint('J29')['n'], 'the selection of J29'),
    ('06_inside', r'is at most $1.2\times10^{-6}$', r'1.2\times10^{-6}', 'ub', lambda: pt(D('0.25') - D(repr(joint('J29')['min_margin']))), 'J29: largest |r0|'),
    ('06_inside', r'on all $200$ engine states inside', '200', 'eq', lambda: pt(EGD['n_in_O']), 'engine states inside O'),
    ('06_inside', r'at most $2.9\times10^{-5}$ (the engine', r'2.9\times10^{-5}', 'ub', lambda: pt(EGD['max_abs_r0_corr_in_O']), 'largest |r0| inside O'),
    # assembly
    ('07_assembly', r'(a) applies ($1.5\times10^{-14}\le10^{-4}$)', r'1.5\times10^{-14}', 'cond', lambda: eic_ok(*B_ALL), 'range end'),
    ('07_assembly', r'(a) applies ($1.5\times10^{-14}\le10^{-4}$)', r'10^{-4}', 'cond', lambda: D('1.5e-14') <= D('1e-4'), 'inside the range of (a)'),
    ('07_assembly', r'$B\le\sigma\max(280,1268\kb)=1268\sigma\kb$', '280', 'eq', lambda: pt(D(2) * 140), '2 * 140'),
    ('07_assembly', r'$B\le\sigma\max(280,1268\kb)=1268\sigma\kb$', '1268', 'ub', lambda: C224, 'C_O'),
    ('07_assembly', r'and $A<112\sigma^2$', '112', 'cond', lambda: ehc(THR)[0] > 0, 'threshold'),
    ('07_assembly', r"because $\psi'/\psi<2.5/L\le1$", '2.5', 'cond', lambda: eic_ok('psi.ratio4') and EIC['grid_psi_ratio_ok'], "psi'/psi < 2.5/L (interval, L = 4 and grid)"),
    ('07_assembly', r'$\Lb=66.644$', '66.644', 'round', lambda: eic('b.Lbar'), 'L_bar at 1.5e-14'),
    ('07_assembly', r'$\kb\le984.16$', '984.16', 'ub', lambda: eic('b.kbar'), 'k_bar at 1.5e-14'),
    ('07_assembly', r'$J\le1.8719\times10^{-8}\le J_0$', r'1.8719\times10^{-8}', 'ub', lambda: eic('b.J'), 'J inside O at 1.5e-14'),
    ('07_assembly', r'for $\kb\ge984$ it is below', '984', 'lb', lambda: eic('b.kbar'), 'k_bar at 1.5e-14 (k_bar decreases with sigma)'),
    ('07_assembly', r'it is below $3.1\times10^{-5}/\kb$', r'3.1\times10^{-5}', 'ub', lambda: eic('b.bk'), 'O-term times k_bar'),
    ('07_assembly', r'is above $0.12/\kb$', '0.12', 'lb', lambda: pt(TV()), 'fourth term of (a) times k_bar'),
    ('07_assembly', r'(c) At $\sigma=10^{-15}$, $\Lb=72.060$', '72.060', 'round', lambda: eic('c.Lbar'), 'L_bar at 1e-15'),
    ('07_assembly', r'and $\kb\le1072.3$;', '1072.3', 'ub', lambda: eic('c.kbar'), 'k_bar at 1e-15'),
    ('07_assembly', r'\ge2.996\times10^{-5}/\kb$', r'2.996\times10^{-5}', 'lb', lambda: eic('c.bk'), 'O-term times k_bar'),
    ('07_assembly', r'For $\Lb\ge72$,', '72', 'lb', lambda: eic('c.Lbar'), 'L_bar at 1e-15 (L_bar increases as sigma decreases)'),
    ('07_assembly', r'(4\Lb+3.57)(\ln2\Lb)^{3/2}/\Lb^{5/2}\le0.0735', '0.0735', 'ub', lambda: ehc('psi(L) <= 289.1 L^{5/2}/(ln 2L)^{3/2} at L = 72: (4L + 3.57)(ln 2L)^{3/2}/L^{5/2}'), 'at L = 72'),
    ('07_assembly', r'\kb\le\sqrt{289.1}\,\Lb^{5/4}', '289.1', 'ub', lambda: pt(D(289) + D('0.0735')), '289 + 0.0735'),
    ('07_assembly', r'\le17.01\,\Lb^{5/4}/(\ln2\Lb)^{3/4}$', '17.01', 'ub', lambda: ehc('sqrt(289.1)'), 'sqrt 289.1'),
    ('07_assembly', r'$0.75\,\alpha/(1268\cdot17.01)\ge1.76\times10^{-6}$', '0.75', 'eq', lambda: pt(D(3) / 4), '3/4'),
    ('07_assembly', r'$0.75\,\alpha/(1268\cdot17.01)\ge1.76\times10^{-6}$', r'1.76\times10^{-6}', 'lb', lambda: ehc('0.75 alpha/1268/17.01  (>= 1.7e-6)'), '0.75 alpha/(1268 17.01)'),
    ('07_assembly', r'$101$ claims in total', '101', 'eq', lambda: pt(EHC['n'] + EIC['n']), 'eh_constants + ei_constants'),
    ('07_assembly', r'interval arithmetic at $60$ digits', '60', 'eq', lambda: pt(int(re.search(r'(\d+) digits', EHC['arithmetic']).group(1))), 'precision of eh_iv'),
    ('07_assembly', r'where $C_O=1268$ comes', '1268', 'ub', lambda: C224, 'C_O'),
    ('07_assembly', r'contains $C_Q=289$', '289', 'ub', lambda: ehc('36.1 * 8  (Q2 <= F(1-s)/s, 1/s = 8L^2/ln 2L)'), 'C_Q'),
    ('07_assembly', r'fixed by $C_h=562$', '562', 'ub', lambda: ehc('70.2 * 8  (E[(d_e^2 p/p)^2] <= G(1-s)/s)'), 'C_h'),
    ('07_assembly', r'(best $J_0=1.94\times10^{-8}$)', r'1.94\times10^{-8}', 'round', lambda: pt([r for r in RNG['rows'] if r['scenario'].startswith('as proved')][0]['J_0']), 'best J0'),
    ('07_assembly', r'& $1.6\times10^{-14}$', r'1.6\times10^{-14}', 'round', lambda: pt(rng('as proved')), 'range'),
    ('07_assembly', r'$C_O/10$ & $1.7\times10^{-13}$', r'1.7\times10^{-13}', 'round', lambda: pt(rng('C_O / 10')), 'range'),
    ('07_assembly', r'$C_Q/10$ & $5.1\times10^{-14}$', r'5.1\times10^{-14}', 'round', lambda: pt(rng('C_Q / 10')), 'range'),
    ('07_assembly', r'$C_h/10$ & $2.5\times10^{-13}$', r'2.5\times10^{-13}', 'round', lambda: pt(rng('C_h / 10')), 'range'),
    ('07_assembly', r'$C_O/10$ and $C_h/10$ & $2.7\times10^{-12}$', r'2.7\times10^{-12}', 'round', lambda: pt(rng('C_O / 10 and C_h / 10')), 'range'),
    ('07_assembly', r'sharp $T1$ & $2.0\times10^{-12}$', r'2.0\times10^{-12}', 'round', lambda: pt(rng('sharp T1: ')), 'range'),
    ('07_assembly', r'sharp $T1$, $C_O/10$, $C_h/10$ & $3.1\times10^{-10}$', r'3.1\times10^{-10}', 'round', lambda: pt(rng('sharp T1, C_O / 10, C_h / 10')), 'range'),
    ('07_assembly', r'sharp $T1$, $C_O/100$, $C_h/100$ & $4.9\times10^{-8}$', r'4.9\times10^{-8}', 'round', lambda: pt(rng('sharp T1, C_O / 100, C_h / 100')), 'range'),
    ('07_assembly', r'the actual $J$ is about $10^{-4}$ of', r'10^{-4}', 'oom', lambda: pt(1 - joint('J20')['min_margin']), 'J20: largest J/(112 s^2 + 1268 s k_bar)'),
    ('07_assembly', r'$(4\Lb+3.57)(\ln2\Lb)^{3/2}', '3.57', 'ub', lambda: (C21[0] / 2, C21[1] / 2), '3.57'),
    ('06_inside', r'|R'+"'"+r'|\le\tfrac{\sqrt{15}}6\sigma^2', '15', 'eq', lambda: pt(D(15)), 'E Z^6 = 15; sqrt15/6 is checked as 0.6455'),
    ('04_three_circles', r'$\delta\le1/16$', '16', 'eq', lambda: pt(D(4) * 4), 'delta = 1/(4L) at L = 4'),
    ('05_outside', r'\sigma^2\tilde C_w^2B^2\le16\,\theta B\,\keff^2A', '16', 'eq', lambda: pt(D(4) ** 2), 'the factor 4 of Paper F, squared'),
    ('06_inside', r'e(16/9)^{1/4}(4L_e/\ln2L_e)', '16', 'eq', lambda: pt(D(16)) if ehc_ok('outer, isotropic: T^3/(T^2 - (T/2)^2)^2 = 16/(9T) <= 16/9') else None, '16/(9T) <= 16/9 (interval check of App. C)'),
    # outlook
    ('08_outlook', r'at most $2.700$ on $2921$ one-variable states', '2.700', 'ub', lambda: pt(max(EGC[p]['families'][f]['rho_star_max'] for p in ('1d', '1ds') for f in EGC[p]['families'])), 'largest rho*, one variable'),
    ('08_outlook', r'at most $2.700$ on $2921$ one-variable states', '2921', 'eq', lambda: pt(EGC['1d']['n_states']), 'one-variable states'),
    ('08_outlook', r'and $2.648$ on $2449$ two-dimensional ones', '2.648', 'ub', lambda: pt(EF2['all']['rho_star_max']), 'largest rho*, two dimensions'),
    ('08_outlook', r'and $2.648$ on $2449$ two-dimensional ones', '2449', 'eq', lambda: pt(EF2['n_states']), 'two-dimensional states'),
    ('08_outlook', r'and with the factor $36.1\sqrt{L/\ln2L}$', '36.1', 'ub', lambda: ehc('total: 2e (4e/sqrt3) + 2.55/sqrt(4/ln8)'), 'total'),
    ('08_outlook', r'the state-independent $(\sqrt{960}/6)\sigma^2$', '960', 'eq', lambda: pt(D(2) ** 6 * 15), '2^6 * 15'),
    ('08_outlook', r'so the threshold $112\sigma^2$', '112', 'cond', lambda: ehc(THR)[0] > 0, 'threshold'),
    ('08_outlook', r'the constant $1268$ of Lemma', '1268', 'ub', lambda: C224, 'C_O'),
    # appendix B
    ('B_borrowed', r'$|R_L|\le(\sqrt{960}/6)\sigma^2$', '960', 'eq', lambda: pt(D(2) ** 6 * 15), '2^6 * 15'),
    ('B_borrowed', r'$\E\zeta^6\le2^6\cdot15\sigma^6$', '15', 'eq', lambda: pt(D(15)), 'E Z^6 = 15 (Gaussian sixth moment)'),
    ('B_borrowed', r'=3\cdot2^{5/4}A\le7.136A$', '7.136', 'ub', lambda: C21, '3 * 2^{5/4}'),
    ('B_borrowed', r'$\E[\st_v^2\Gamma]\ge B-A-1.291\sqrt A$', '1.291', 'ub', lambda: ehc('sqrt15/3  (|E[s_v^3 s_w]|/3 <= sqrt(15 A)/3)'), 'sqrt15/3'),
    ('B_borrowed', r'$\E\Gamma^4\le(105\cdot1680)^{1/2}=420$', '105', 'eq', lambda: pt(D(105)), 'E Z^8 = 105'),
    ('B_borrowed', r'$\E\Gamma^4\le(105\cdot1680)^{1/2}=420$', '1680', 'eq', lambda: pt(D(105) * 16), 'E (Z_v + Z_w)^8 = 105 * 2^4'),
    ('B_borrowed', r'$\E\Gamma^4\le(105\cdot1680)^{1/2}=420$', '420', 'eq', lambda: pt((D(105) * 1680).sqrt()), 'sqrt(105 * 1680)'),
    ('B_borrowed', r'(2\sqrt{1260})^{1/2}\sqrt{\theta B}\le8.4258', '1260', 'eq', lambda: pt(D(3) * 420), '3 * 420'),
    ('B_borrowed', r'(2\sqrt{1260})^{1/2}\sqrt{\theta B}\le8.4258', '8.4258', 'ub', lambda: ehc('c_2 = sqrt(2 sqrt1260)  (E G^4 <= sqrt(105 * 1680) = 420)'), 'sqrt(2 sqrt1260)'),
    # appendix C
    ('C_constants', r'(\code{src/eh\_iv.py}): $60$ digits', '60', 'eq', lambda: pt(int(re.search(r'(\d+) digits', EHC['arithmetic']).group(1))), 'precision'),
    ('C_constants', r'($L=4$ to $2\times10^7$', r'2\times10^7', 'round', lambda: pt(EHC['grid_L'][-1]['L']), 'largest L of the grid'),
    ('C_constants', r'$\sigma=10^{-15}$ to $10^{-300}$', r'10^{-15}', 'eq', lambda: pt(D(10) ** EHC['theorem_rows'][0]['sigma_exp']), 'grid'),
    ('C_constants', r'$\sigma=10^{-15}$ to $10^{-300}$', r'10^{-300}', 'eq', lambda: pt(D(10) ** EHC['theorem_rows'][-1]['sigma_exp']), 'grid'),
    ('C_constants', r'All $101$ claims', '101', 'eq', lambda: pt(EHC['n'] + EIC['n']), 'number of claims'),
    ('C_constants', r'($\sqrt{4/\ln8}=1.38693\ldots$', r'1.38693\ldots', 'trunc', lambda: ehc('sqrt(4/ln 8) >= 1.3869 (the correct lower bound)'), 'sqrt(4/ln 8)'),
    ('C_constants', r'printed as $1.387$ and used', '1.387', 'eq', lambda: pt([c for c in EHC['printed'] if 'sqrt(4/ln 8) >= 1.387' in c['what']][0]['claim']), 'the slip as printed in the record'),
    ('C_constants', r'($4/\ln8=1.92359\ldots$', r'1.92359\ldots', 'trunc', lambda: ehc('4/ln 8 >= 1.9235 (the correct lower bound)'), '4/ln 8'),
    ('C_constants', r'printed as $1.924$)', '1.924', 'eq', lambda: pt([c for c in EHC['printed'] if '4/ln 8 >= 1.924' in c['what']][0]['claim']), 'the slip as printed'),
    ('C_constants', r'the constants $36.1$ and $70.2$', '36.1', 'ub', lambda: ehc('total: 2e (4e/sqrt3) + 2.55/sqrt(4/ln8)'), 'total'),
    ('C_constants', r'the constants $36.1$ and $70.2$', '70.2', 'ub', lambda: ehc('total: 2e * 12.56 + 2.91/(4/ln 8)'), 'isotropic total'),
    # appendix D
    ('D_joints', r'\item $692$ states of Paper', '692', 'eq', lambda: pt(JO['n_engine']), 'engine states'),
    ('D_joints', r'of which $200$ lie inside', '200', 'eq', lambda: pt(JO['n_engine_in_O']), 'engine states inside O'),
    ('D_joints', r'\item $2449$ two-dimensional states', '2449', 'eq', lambda: pt(EF2['n_states']), 'two-dimensional states'),
    ('D_joints', r'\item $2921$ one-variable states', '2921', 'eq', lambda: pt(EGC['1d']['n_states']), 'one-variable states'),
    ('D_joints', r'which is up to $0.5\%$ of the main terms', '0.5', 'cond', lambda: joint('J1')['min_margin'] >= 0, 'J1: the defect is within the accuracy filter 0.5 %'),
    ('D_joints', r'At $\sigma=10^{-24}$ that defect', r'10^{-24}', 'cond', lambda: 1e-24 in set(r['sigma'] for r in rec('results_EH/eh_states.json')['rows']), 'an engine sigma'),
    ('D_joints', r'it fails by a factor of $9\times10^{71}$', r'9\times10^{71}', 'round', lambda: pt(-joint('J18')['min_margin']), 'J18'),
    ('D_joints', r'all have $\sigma\ge10^{-6}$', r'10^{-6}', 'eq', lambda: pt(min(x['sigma'] for x in EGD['rows'] if x.get('in_O', x.get('inO')))), 'smallest sigma inside O'),
    ('D_joints', r'$\sigma=10^{-12},10^{-24},10^{-36}$ lie outside', r'10^{-12}', 'cond', lambda: SIG_OUT(1e-12), 'engine states at this sigma, none inside O'),
    ('D_joints', r'$\sigma=10^{-12},10^{-24},10^{-36}$ lie outside', r'10^{-24}', 'cond', lambda: SIG_OUT(1e-24), 'engine states at this sigma, none inside O'),
    ('D_joints', r'$\sigma=10^{-12},10^{-24},10^{-36}$ lie outside', r'10^{-36}', 'cond', lambda: SIG_OUT(1e-36), 'engine states at this sigma, none inside O'),
    # appendix E
    ('E_verification', r'up to a relative $10^{-9}$', r'10^{-9}', 'eq', common_rtol, 'FLOAT_RTOL of verification/common.py'),
    # chain
    ('03_chain', r'fails there by a factor of $9\times10^{71}$', r'9\times10^{71}', 'round', lambda: pt(-joint('J18')['min_margin']), 'J18'),
]


# wrong versions: the check must reject each of them (the slips of the working record and nearby wrong roundings)
WRONG = [
    ('2.647', 'ub', lambda: pt(EF2['all']['rho_star_max']), 'largest rho* in two dimensions, rounded the wrong way (first draft of Sec. 8)'),
    ('1.387', 'lb', lambda: ehc('sqrt(4/ln 8) >= 1.3869 (the correct lower bound)'), 'sqrt(4/ln 8) as printed in the working record'),
    ('1.924', 'lb', lambda: ehc('4/ln 8 >= 1.9235 (the correct lower bound)'), '4/ln 8 as printed in the working record'),
    ('9.7167', 'ub', lambda: ehc_printed('9.7167'), 'the third slip of the working record (EE-C, not used)'),
    ('0.99999', 'lb', lambda: pt(joint('J3')['min_margin']), 'J3 margin rounded up'),
    ('984.15', 'ub', lambda: eic('b.kbar'), 'k_bar at 1.5e-14 rounded down'),
    ('1810', 'lb', lambda: pt(1 + joint('J22')['min_margin']), 'J22 factor rounded up'),
    ('1.2496656113', 'round', lambda: pt(EIP['grid_b']['2A+N']), "the reviewer's 2A + N, which our quadrature does not reproduce beyond 6 digits"),
    ('1266', 'ub', lambda: C224, 'C_O below 8 sqrt(224 112)'),
    (r'1.6\times10^{-14}', 'cond', lambda: eic_ok(*B_ALL) and rng('as proved') >= D('1.6e-14'), 'a range end beyond the solved sigma_0'),
]


def ehc_printed(claim):
    return [(D(repr(c['enclosure'][0])), D(repr(c['enclosure'][1]))) for c in EHC['printed'] if repr(c['claim']) == claim][0]


def SIG_OUT(s):
    rows = rec('results_EH/eh_states.json')['rows']; at = [r for r in rows if r['sigma'] == s]
    return len(at) > 0 and not any(r['T1']['A'] + max(r['theta'], 0.0) * r['T1']['B'] / 2 < 112 * s * s for r in at)


# numbers of the running text that are not results: (file, printed, reason)
ALLOW = [
    ('02_main_theorem', '0.1', 'c_r = sqrt(0.1): a choice of the proof (region boundaries of Proposition 15)'),
    ('05_outside', '0.1', 'c_r = sqrt(0.1): a choice'),
    ('E_verification', '0.1', 'c_r = sqrt(0.1), quoted as an example of a definition'),
    ('07_assembly', '10', 'labels of the scenarios C/10 of Table 2'),
    ('07_assembly', '100', 'labels of the scenarios C/100 of Table 2'),
    ('07_assembly', r'10^{-15}', 'the point at which (c) is evaluated (checked as cond in 02 and by c.* of ei_constants)'),
]


def parse(s):
    s = s.replace(r'\ldots', '')
    m = re.fullmatch(r'(\d+(?:\.\d+)?)\\times10\^\{(-?\d+)\}', s) or re.fullmatch(r'(\d+(?:\.\d+)?)\\times10\^(\d)', s)
    if m: return D(m.group(1)) * D(10) ** int(m.group(2)), D(1) * D(10) ** (int(m.group(2)) - (len(m.group(1).split('.')[1]) if '.' in m.group(1) else 0))
    m = re.fullmatch(r'10\^\{(-?\d+)\}', s)
    if m: return D(10) ** int(m.group(1)), D(10) ** int(m.group(1))
    return D(s), D(10) ** (-(len(s.split('.')[1]) if '.' in s else 0))


def judge(kind, P, ulp, v):
    if kind == 'cond': return bool(v), ''
    lo, hi = v
    if kind == 'ub': return hi <= P, 'true <= %s' % hi
    if kind == 'lt': return hi < P, 'true < %s' % hi
    if kind == 'lb': return lo >= P, 'true >= %s' % lo
    if kind == 'eq': return lo == P == hi, 'true = %s' % lo
    if kind == 'round': return abs(lo - P) <= ulp / 2 * (1 + D('1e-12')) and abs(hi - P) <= ulp / 2 * (1 + D('1e-12')), 'true = %.6g' % lo
    if kind == 'trunc': return P <= lo and hi < P + ulp, 'true = %s' % str(lo)[:14]
    if kind == 'oom': return abs(D(math.log10(lo)) - D(math.log10(P))) <= D('0.5'), 'log10 true = %.3f' % math.log10(lo)
    raise ValueError(kind)


TOK = re.compile(r'(?<![\w.])(\d+\.\d+(?:\\ldots)?(?:\\times10\^\{-?\d+\})?|\d+\\times10\^\{?-?\d+\}?|10\^\{-?\d+\}|\d{2,})(?![\w.])')


def tokens(tex):
    t = re.sub(r'(?m)(?<!\\)%.*$', '', tex)
    t = re.sub(r'\\(?:ref|eqref|label|cite|input|includegraphics|code|url|href)\*?(?:\[[^\]]*\])?\{[^}]*\}', ' ', t)
    t = re.sub(r'\\[dt]?frac\d\d', ' ', t)                         # \frac12, \tfrac34: fractions of single digits
    t = t.replace(r'\times', '\x00').replace(r'\ldots', '\x01')
    t = re.sub(r'\\[A-Za-z]+\*?', ' ', t)
    t = t.replace('\x00', r'\times').replace('\x01', r'\ldots')
    return TOK.findall(t)


def main():
    texts = {os.path.basename(f)[:-4]: open(f, encoding='utf-8').read() for f in sorted(os.listdir(os.path.join(PAPER, 'sections')) and
             [os.path.join(PAPER, 'sections', x) for x in os.listdir(os.path.join(PAPER, 'sections')) if x.endswith('.tex')])}
    rows = []; ok = True
    for f, ctx, printed, kind, src, what in ENTRIES:
        found = f in texts and ctx in texts[f] and printed in ctx
        try:
            P, ulp = parse(printed) if kind != 'cond' else (None, None)
            good, info = judge(kind, P, ulp, src())
        except Exception as e:  # a missing record entry is a failure, not a crash
            good, info = False, 'error: %r' % (e,)
        rows.append({'file': f, 'printed': printed, 'kind': kind, 'what': what, 'context_found': bool(found), 'holds': bool(good), 'info': info})
        ok = ok and found and good
    wrong = []
    for printed, kind, src, what in WRONG:
        P, ulp = parse(printed) if kind != 'cond' else (None, None)
        try: good, info = judge(kind, P, ulp, src())
        except Exception as e: good, info = True, 'error: %r' % (e,)     # an error does not count as caught
        wrong.append({'printed': printed, 'kind': kind, 'what': what, 'caught': not good, 'info': info})
    # inventory
    covered = set((f, p) for f, c, p, k, s, w in ENTRIES) | set((f, p) for f, p, r in ALLOW)
    inv = []; n_tok = 0
    for f, t in texts.items():
        if f + '.tex' in SKIP: continue
        for tok in sorted(set(tokens(t))):
            n_tok += 1
            if (f, tok) not in covered: inv.append({'file': f, 'number': tok})
    allow_unused = [{'file': f, 'number': p} for f, p, r in ALLOW if p not in set(tokens(texts.get(f, '')))]
    out = {'paper': os.path.relpath(PAPER, ROOT) if PAPER.startswith(ROOT) else 'external', 'n_entries': len(rows), 'n_hold': sum(r['holds'] and r['context_found'] for r in rows),
           'by_kind': {k: sum(1 for r in rows if r['kind'] == k) for k in sorted(set(r['kind'] for r in rows))},
           'n_numbers_in_text': n_tok, 'uncovered': inv, 'allow': [{'file': f, 'number': p, 'reason': r} for f, p, r in ALLOW], 'allow_unused': allow_unused,
           'failures': [r for r in rows if not (r['holds'] and r['context_found'])], 'wrong_versions': wrong,
           'wrong_versions_caught': sum(x['caught'] for x in wrong), 'rows': rows}
    out['ok'] = bool(ok and not inv and not allow_unused and all(x['caught'] for x in wrong))
    os.makedirs(os.path.join(ROOT, 'results_EI'), exist_ok=True)
    json.dump(out, open(os.path.join(ROOT, 'results_EI', 'ei_text_numbers.json'), 'w', encoding='utf-8'), indent=1)
    for r in out['failures']: print('FAIL', r)
    for r in inv: print('UNCOVERED', r)
    for r in allow_unused: print('ALLOW UNUSED', r)
    for x in wrong:
        if not x['caught']: print('WRONG VERSION NOT CAUGHT', x)
    print('%d/%d printed numbers hold; %d distinct numbers in the text, %d uncovered; %d/%d wrong versions caught; ok %s'
          % (out['n_hold'], out['n_entries'], n_tok, len(inv), out['wrong_versions_caught'], len(wrong), out['ok']))
    return 0 if out['ok'] else 1


if __name__ == '__main__':
    sys.exit(main())
