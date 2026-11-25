"""Stage EH: a small, rigorous interval arithmetic on Python's decimal module (mpmath is not installed in this environment).
Every operation rounds outward: +, -, *, / with the context rounding ROUND_FLOOR for the lower end and ROUND_CEILING for the upper end;
sqrt, exp, ln are correctly rounded by decimal (ROUND_HALF_EVEN, error <= 1/2 ulp) and the result is then widened by one ulp on each side.
pi and e enter as enclosures from their decimal expansions (pi to 60 digits; e = exp(1) as above).  Python 3.11 grammar."""
from decimal import Decimal as Dc, localcontext, ROUND_FLOOR, ROUND_CEILING, ROUND_HALF_EVEN, getcontext
from fractions import Fraction

PREC = 60
PI_S = '3.14159265358979323846264338327950288419716939937510582097494'      # 60 digits after the point: |pi - PI_S| < 1e-59


def _d(x):
    if isinstance(x, Dc): return x
    if isinstance(x, Fraction):
        with localcontext() as c: c.prec = PREC; return Dc(x.numerator) / Dc(x.denominator)
    return Dc(str(x)) if isinstance(x, float) else Dc(x)


def _op(f, a, b, mode):
    with localcontext() as c:
        c.prec = PREC; c.rounding = mode; return f(a, b)


class Iv:
    __slots__ = ('lo', 'hi')

    def __init__(self, lo, hi=None):
        if isinstance(lo, Iv): self.lo, self.hi = lo.lo, lo.hi; return
        if isinstance(lo, Fraction) and hi is None:                 # an exact rational: enclose its decimal value
            self.lo = _op(lambda a, b: a / b, Dc(lo.numerator), Dc(lo.denominator), ROUND_FLOOR); self.hi = _op(lambda a, b: a / b, Dc(lo.numerator), Dc(lo.denominator), ROUND_CEILING); return
        self.lo = _d(lo); self.hi = _d(lo if hi is None else hi)
        assert self.lo <= self.hi

    @staticmethod
    def c(x): return x if isinstance(x, Iv) else Iv(Fraction(x) if isinstance(x, int) else x)

    def __add__(s, o): o = Iv.c(o); return Iv(_op(Dc.__add__, s.lo, o.lo, ROUND_FLOOR), _op(Dc.__add__, s.hi, o.hi, ROUND_CEILING))
    __radd__ = __add__
    def __neg__(s): return Iv(-s.hi, -s.lo)
    def __sub__(s, o): return s + (-Iv.c(o))
    def __rsub__(s, o): return Iv.c(o) - s

    def __mul__(s, o):
        o = Iv.c(o); ps = [(a, b) for a in (s.lo, s.hi) for b in (o.lo, o.hi)]
        return Iv(min(_op(Dc.__mul__, a, b, ROUND_FLOOR) for a, b in ps), max(_op(Dc.__mul__, a, b, ROUND_CEILING) for a, b in ps))
    __rmul__ = __mul__

    def __truediv__(s, o):
        o = Iv.c(o); assert o.lo > 0 or o.hi < 0, 'division by an interval containing 0'
        ps = [(a, b) for a in (s.lo, s.hi) for b in (o.lo, o.hi)]
        return Iv(min(_op(Dc.__truediv__, a, b, ROUND_FLOOR) for a, b in ps), max(_op(Dc.__truediv__, a, b, ROUND_CEILING) for a, b in ps))
    def __rtruediv__(s, o): return Iv.c(o) / s

    def _mono(s, f, inc=True):
        with localcontext() as c:
            c.prec = PREC; c.rounding = ROUND_HALF_EVEN; a, b = f(s.lo), f(s.hi)
            if not inc: a, b = b, a
            return Iv(a.next_minus(c), b.next_plus(c))
    def sqrt(s): assert s.lo >= 0; return s._mono(lambda x: x.sqrt())
    def exp(s): return s._mono(lambda x: x.exp())
    def ln(s): assert s.lo > 0; return s._mono(lambda x: x.ln())
    def __pow__(s, a):
        if isinstance(a, int) and a >= 0:
            r = Iv(1)
            for _ in range(a): r = r * s
            return r
        return (Iv.c(a) * s.ln()).exp()

    def __le__(s, o): o = Iv.c(o); return s.hi <= o.lo          # certainly <=
    def __lt__(s, o): o = Iv.c(o); return s.hi < o.lo
    def __ge__(s, o): return Iv.c(o) <= s
    def __gt__(s, o): return Iv.c(o) < s
    def f(s): return (float(s.lo), float(s.hi))
    def __repr__(s): return '[%.12g, %.12g]' % (float(s.lo), float(s.hi))


PI = Iv(Dc(PI_S) - Dc('1e-59'), Dc(PI_S) + Dc('1e-59'))
E = Iv(1).exp()


def sqrt(x): return Iv.c(x).sqrt()
def exp(x): return Iv.c(x).exp()
def ln(x): return Iv.c(x).ln()
