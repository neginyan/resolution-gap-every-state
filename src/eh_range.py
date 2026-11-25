"""Stage EH, 1.4 (derivation, numerical): how the range sigma <= sigma_0 = 1e-15 of Theorem EG (b) depends on the constants of the chain.
The O-part of the proof needs  J_max(sigma) := 112 sigma^2 + C_O sigma k_bar(sigma) <= J_0,  where
  k_bar(sigma)^2 = psi(L_bar) = C_Q L^{5/2}/(ln 2L)^{3/2} + 4L + 3.57  (C_Q = 289, Theorem EG-1; or a sharp T1: c0 L + C0),  L_bar = 2 ln(1/sigma) + ln(1/alpha),
  J_0 = the largest J with  h(J)/alpha + sqrt(h(J))/(2 sqrt alpha) + 0.6455 sqrt(J)/alpha <= 1/4,  h(J) = C_h J L_J^3/(ln 2 L_J)^2  (C_h = 562; Theorem EG-O),
  C_O = 1268 (Lemma EF-O).
sigma_0 is found by bisection in ln sigma (the left side decreases with sigma).  Output results_EH/eh_range.json.  Python 3.11 grammar."""
import os, sys, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ed_common as C

AL = 1 / (2 * math.pi ** 2)


def J0(Ch):
    def br(J):
        L = math.log(1 / J); h = Ch * J * L ** 3 / math.log(2 * L) ** 2
        return h / AL + math.sqrt(h) / (2 * math.sqrt(AL)) + 0.6455 * math.sqrt(J) / AL
    lo, hi = math.log(1e-30), math.log(math.exp(-4))
    for _ in range(200):
        m = (lo + hi) / 2; (lo, hi) = (m, hi) if br(math.exp(m)) <= 0.25 else (lo, m)
    return math.exp(lo)


def sigma0(CO=1268.0, CQ=289.0, Ch=562.0, sharp=None):
    j0 = J0(Ch)
    def kb(s):
        Lb = 2 * math.log(1 / s) + math.log(2 * (math.pi ** 2 + s * s))
        return math.sqrt(sharp[0] * Lb + sharp[1]) if sharp else math.sqrt(CQ * Lb ** 2.5 / math.log(2 * Lb) ** 1.5 + 4 * Lb + 3.57)
    f = lambda ls: 112 * math.exp(2 * ls) + CO * math.exp(ls) * kb(math.exp(ls)) - j0
    lo, hi = math.log(1e-300), math.log(1e-1)
    if f(hi) <= 0: return 1e-1, j0
    for _ in range(200):
        m = (lo + hi) / 2; (lo, hi) = (m, hi) if f(m) <= 0 else (lo, m)
    return math.exp(lo), j0


def main():
    rows = []
    for name, kw in (('as proved (C_O 1268, C_Q 289, C_h 562)', {}), ('C_O / 10', dict(CO=126.8)), ('C_Q / 10', dict(CQ=28.9)), ('C_h / 10', dict(Ch=56.2)),
                     ('C_O / 10 and C_h / 10', dict(CO=126.8, Ch=56.2)), ('sharp T1: k_eff^2 <= L + 4 (c0 = 1)', dict(sharp=(1.0, 4.0))),
                     ('sharp T1 and C_O / 10', dict(sharp=(1.0, 4.0), CO=126.8)), ('sharp T1, C_O / 10, C_h / 10', dict(sharp=(1.0, 4.0), CO=126.8, Ch=56.2)),
                     ('sharp T1, C_O / 100, C_h / 100', dict(sharp=(1.0, 4.0), CO=12.68, Ch=5.62))):
        s0, j0 = sigma0(**kw); rows.append({'scenario': name, 'sigma_0': s0, 'log10_sigma_0': math.log10(s0), 'J_0': j0})
        print('%-40s sigma_0 = %.3e (log10 %.2f)   J_0 = %.3e' % (name, s0, math.log10(s0), j0))
    out = {'rows': rows, 'env': C.env(), 'ok': True, 'proved_sigma_0_at_least_1e-15': bool(rows[0]['sigma_0'] >= 1e-15)}
    C.save('eh_range.json', out, C.RES_EH)


if __name__ == '__main__':
    main()
