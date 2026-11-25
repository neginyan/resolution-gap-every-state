"""Stage EI (paper K): the chain figure (Fig. 1): which parts Theorem 1 rests on.  usage: python ei_figure.py OUTDIR -> fig_chain_K.pdf, fig_chain_K.png.
Python 3.11 grammar.  Model: Claude Opus 5.5 (claude-opus-5-5)."""
import os, sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

NEW = '#f2c14e'; OWN = '#cfe3f3'; BOR = '#e3e3e3'
X = dict(A=0.09, B=0.29, C=0.49, D=0.69, E=0.89); Yr = {1: 0.84, 2: 0.67, 3: 0.50, 4: 0.33, 5: 0.14}
def P(c): return X[c[0]], Yr[int(c[1])]
NODES = {  # key: (cell, text, colour)
    'K': ('A1', 'identities,\nhalf lemma', OWN),
    'sm': ('B1', 'smoothing\n$Q_2\\leq F(1-u)/u$', OWN),
    'gap': ('C1', 'gap formula,\nbranch identity (E)', BOR),
    'G2': ('D1', 'remainder of F\n(Lemma G2)', BOR),
    'win': ('E1', 'window inequality\n(F, AT2)', BOR),
    'circ': ('A2', 'three circles\n$k_{\\rm eff}^2\\leq\\psi(L)$', NEW),
    'EGD': ('C2', 'expansion around\nthe observed point', NEW),
    'EF1': ('D2', 'numerator', OWN),
    'EEC': ('E2', '$\\tilde C_w$ (EE)', OWN),
    'iso': ('A3', 'isotropic circles,\ndiagonal smoothing', NEW),
    'kb': ('C3', '$k_{\\rm eff}\\leq\\bar k$', OWN),
    'EF2': ('D3', 'threshold\n$M\\geq112\\sigma^2$', NEW),
    'EGO': ('A4', 'small information\n$J\\leq1.9\\times10^{-8}$', NEW),
    'out': ('D4', 'outside $\\mathcal{O}$:\nfour regions', OWN),
    'EFO': ('E4', 'inside $\\mathcal{O}$:\n$B\\leq1268\\,\\sigma k_{\\rm eff}$', OWN),
    'main': ('C5', 'Theorem 1:\nevery state', NEW),
}
EDGES = [('gap', 'G2'), ('G2', 'EF1'), ('EF1', 'EF2'), ('K', 'circ'), ('sm', 'circ'), ('win', 'EEC'), ('gap', 'EGD'), ('circ', 'iso'), ('circ', 'kb'),
         ('EGD', 'EGO'), ('iso', 'EGO'), ('EF2', 'out'), ('EEC', 'out'), ('kb', 'out'), ('EEC', 'EFO'), ('kb', 'EFO'),
         ('EGO', 'main'), ('out', 'main'), ('EFO', 'main')]
W, H = 0.17, 0.10
ROUTE = {('kb', 'out'): (0.425, 0.49), ('kb', 'EFO'): (0.425, 0.49), ('EEC', 'out'): (0.405, 0.84)}


def border(x0, y0, x1, y1):
    """the point where the segment from the centre (x0, y0) towards (x1, y1) leaves the box around (x0, y0)"""
    dx, dy = x1 - x0, y1 - y0; t = min((W / 2 + 0.006) / abs(dx) if dx else 9, (H / 2 + 0.006) / abs(dy) if dy else 9)
    return x0 + t * dx, y0 + t * dy


def main(out):
    fig, ax = plt.subplots(figsize=(7.0, 5.0)); ax.set_xlim(-0.01, 0.99); ax.set_ylim(0.02, 0.97); ax.axis('off')
    ax.add_patch(FancyBboxPatch((0.005, 0.915), 0.97, 0.04, boxstyle='round,pad=0.004,rounding_size=0.01', fc=BOR, ec='#333333', lw=0.8))
    ax.text(0.49, 0.935, 'used throughout: integration by parts for any $\\rho$ with $E|X|<\\infty$; Tweedie\'s formulas; Cram\u00e9r-Rao floor $A\\geq\\alpha\\sigma^2$ from (H2)', ha='center', va='center', fontsize=7.0)
    for a, b in EDGES:
        x0, y0 = P(NODES[a][0]); x1, y1 = P(NODES[b][0])
        xs, ys = border(x0, y0, x1, y1); xe, ye = border(x1, y1, x0, y0)
        if (a, b) in ROUTE:                       # elbow routes through the channel between rows 3 and 4
            yc, xa_ = ROUTE[(a, b)]; xt = x1 + (0.03 if (a, b) == ('EEC', 'out') else -0.03); pts = [(xa_, y0 - H / 2 - 0.006), (xa_, yc), (xt, yc), (xt, y1 + H / 2 + 0.006)]
            ax.plot([p[0] for p in pts[:-1]], [p[1] for p in pts[:-1]], color='#444444', lw=0.8)
            ax.annotate('', xy=pts[-1], xytext=pts[-2], arrowprops=dict(arrowstyle='-|>', color='#444444', lw=0.8, mutation_scale=8, shrinkA=0, shrinkB=0)); continue
        ax.annotate('', xy=(xe, ye), xytext=(xs, ys), arrowprops=dict(arrowstyle='-|>', color='#444444', lw=0.8, mutation_scale=8, shrinkA=0, shrinkB=0))
    for k, (cell, t, c) in NODES.items():
        x, y = P(cell)
        ax.add_patch(FancyBboxPatch((x - W / 2, y - H / 2), W, H, boxstyle='round,pad=0.004,rounding_size=0.012', fc=c, ec='#333333', lw=0.8, zorder=3))
        ax.text(x, y, t, ha='center', va='center', fontsize=7.0, zorder=4)
    for i, (c, t) in enumerate(((NEW, 'new in this paper'), (OWN, 'lab stages behind this paper'), (BOR, 'borrowed (Papers E, F; standard)'))):
        ax.add_patch(FancyBboxPatch((0.64, 0.16 - 0.04 * i), 0.025, 0.022, boxstyle='square,pad=0', fc=c, ec='#333333', lw=0.6)); ax.text(0.675, 0.171 - 0.04 * i, t, fontsize=6.6, va='center')
    os.makedirs(out, exist_ok=True); fig.tight_layout()
    fig.savefig(os.path.join(out, 'fig_chain_K.pdf')); fig.savefig(os.path.join(out, 'fig_chain_K.png'), dpi=200); print('figure written')


if __name__ == '__main__':
    main(sys.argv[1])
