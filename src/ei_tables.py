"""Stage EI (paper K): write the tables of App. C (constants by interval arithmetic) and App. D (joints on recorded states) of the manuscript from the
records results_EH/eh_constants.json, results_EI/ei_constants.json and results_EH/eh_joints.json.  usage: python ei_tables.py OUTDIR (the sections folder).
Python 3.11 grammar.  Model: Claude Opus 5.5 (claude-opus-5-5)."""
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ed_common as C

LAB2PAPER = {'E3': 'L.\\ref{lem:CR}', 'D2': 'L.\\ref{lem:K}', 'F1': 'L.\\ref{lem:G2}', 'EF1': 'L.\\ref{lem:EF1}', 'EF2': 'L.\\ref{lem:EF2}', 'EEC': 'L.\\ref{lem:EEC}',
             'EFO': 'L.\\ref{lem:EFO}', 'EFI': 'P.\\ref{prop:outside}', 'G1': 'T.\\ref{thm:circles}', 'GI': 'L.\\ref{lem:iso},\\ref{lem:sm}', 'GD': 'L.\\ref{lem:EGD}', 'GO': 'T.\\ref{thm:EGO}',
             'GT': 'T.\\ref{thm:main}', 'EXT': 'T.\\ref{thm:EGO},\\ref{thm:main}'}


# wording of the record that reads badly in a table of the paper (the record itself is not changed)
WHAT = {'-ln(1 - eta) - 2 eta at eta = 0.13 (<= 0: the function is increasing... see report)':
        '-ln(1 - eta) - 2 eta at eta = 0.13 (<= 0; it is <= 0 for every eta in [0, 1/2], since it vanishes at 0 and its derivative 1/(1 - eta) - 2 is negative there)'}

# the rounded values printed in the working record (eh_constants.json, 'printed'): record string -> (place in the working record, LaTeX formula).
# The record is not changed; s of the working record is tau of the paper.
PRINTED = {
    'printed in EE 3.3 (original EE-C, not used by Theorem EG): c1 = 2*315^{1/4} + sqrt15/3 <= 9.7167':
        ('EE 3.3$^\\dagger$', '$c_1=2\\cdot315^{1/4}+\\sqrt{15}/3\\le9.7167$'),
    'printed in EG 2.4 step 7: 1 - s >= 0.9837 at L = 4': ('EG 2.4 step 7', '$1-\\tau\\ge0.9837$ at $L=4$'),
    'printed in EG 4.2 step 4: (1 - s)^3 >= 0.952 at L = 4': ('EG 4.2 step 4', '$(1-\\tau)^3\\ge0.952$ at $L=4$'),
    'printed in EF 2.4: 1.291 sqrt112 <= 13.67': ('EF 2.4', '$1.291\\sqrt{112}\\le13.67$'),
    'printed in EF 2.4: 8.4258 sqrt224 <= 126.11': ('EF 2.4', '$8.4258\\sqrt{224}\\le126.11$'),
    'printed in EF 2.4: 4 sqrt(224 * 112) <= 633.6': ('EF 2.4', '$4\\sqrt{224\\cdot112}\\le633.6$'),
    'printed in EG 4.3: alpha >= 0.0506': ('EG 4.3', '$\\alpha\\ge0.0506$'),
    'printed in EG 4.3: Xi^2/alpha <= 0.0132': ('EG 4.3', '$\\Xi^2/\\alpha\\le0.0132$'),
    'printed in EG 4.3: Xi/(2 sqrt alpha) <= 0.0574': ('EG 4.3', '$\\Xi/(2\\sqrt\\alpha)\\le0.0574$'),
    'printed in EG 2.4 step 7 as a lower bound: sqrt(4/ln 8) >= 1.387': ('EG 2.4 step 7', '$\\sqrt{4/\\ln8}\\ge1.387$'),
    'printed in EG 4.2 step 4 as a lower bound: 4/ln 8 >= 1.924': ('EG 4.2 step 4', '$4/\\ln8\\ge1.924$'),
}


C_GROUPS = 3          # Table 4 (101 rows) in three table* floats of at most 34 rows


def clean(w):
    w = WHAT.get(w, w)
    for p in ('new (section 1.4): ', 'new (section 1.4) '):
        if w.startswith(p): w = w[len(p):]
    return w

def tex(s):
    out = []
    i = 0
    while i < len(s):
        two = s[i:i + 2]
        if two in ('<=', '>='): out.append('$\\le$' if two == '<=' else '$\\ge$'); i += 2; continue
        if two == '->': out.append('$\\to$'); i += 2; continue
        ch = s[i]
        out.append({'\\': '\\textbackslash{}', '_': '\\_', '%': '\\%', '&': '\\&', '#': '\\#', '{': '\\{', '}': '\\}', '^': '\\textasciicircum{}', '|': '\\textbar{}',
                    '<': '$<$', '>': '$>$', '~': '\\textasciitilde{}'}.get(ch, ch)); i += 1
    return ''.join(out)


def fmt(x):
    x = float(x); return ('%.6g' % x).replace('e-0', 'e-').replace('e+0', 'e+')


def main(out):
    eh = json.load(open(os.path.join(C.RES_EH, 'eh_constants.json'), encoding='utf-8')); ei = json.load(open(os.path.join(C.ROOT, 'results_EI', 'ei_constants.json'), encoding='utf-8'))
    jo = json.load(open(os.path.join(C.RES_EH, 'eh_joints.json'), encoding='utf-8'))
    # Table 4 as plain table* floats (no longtable: REVTeX 4.2 with a recent longtable stops on longtable*), split into C_GROUPS groups of rows; the first
    # carries the caption and the label; the others have no caption (so no number of their own) and a heading line "Table IV (continued)".
    side = {'<=': '$\\le$', '>=': '$\\ge$', '<': '$<$', '>': '$>$', '==': '$=$'}
    rows = ['%s & \\texttt{%s} & %s & %s & %s\\\\' % (LAB2PAPER.get(c['item'], c['item']), tex(clean(c['what'])), fmt(c['enclosure'][0]), fmt(c['claim']), side[c['side']]) for c in eh['checks']]
    rows += ['%s & \\texttt{%s} & %s & %s & %s\\\\' % ('T.\\ref{thm:EGO}' if c['key'].startswith('O.') else ('R.\\ref{rem:delta}' if c['key'] == 'd04' else 'T.\\ref{thm:main}'), tex(clean(c['what'])),
                                                     fmt(c['lo']), fmt(c['claim']), side[c['side']]) for c in ei['checks']]
    n = len(rows); k = -(-n // C_GROUPS); groups = [rows[i:i + k] for i in range(0, n, k)]
    head = 'part & claim (the expression, as evaluated by the code) & lower end of the enclosure & claimed & side\\\\'
    L = ['% generated by src/ei_tables.py from results_EH/eh_constants.json and results_EI/ei_constants.json -- do not edit by hand',
         '%% Table 4 in %d groups of at most %d rows (plain table* floats)' % (len(groups), k)]
    for g, grp in enumerate(groups):
        L += ['\\begin{table*}[p]', '\\scriptsize']
        if g == 0:
            L.append('\\caption{Every constant of the chain, recomputed with outward-rounded interval arithmetic (60 digits; \\code{src/eh\\_constants.py}, \\code{src/ei\\_constants.py}). '
                     'The expressions are printed as evaluated by the code; two explanatory remarks of the record are reworded. '
                     'The enclosure is the interval of the exact expression; ``claimed\'\' is the value used in the text, and ``side\'\' says on which side the enclosure must lie. '
                     'All %d hold. The table continues in %d more parts.}\\label{tab:constants}' % (n, len(groups) - 1))
        else:
            L += ['{\\small\\centering Table~\\ref{tab:constants} (continued)\\par}', '\\smallskip']     # no caption: no new number, no counter changes
        L += ['\\begin{ruledtabular}', '\\begin{tabular}{@{}lp{0.50\\textwidth}lll@{}}', head, '\\hline'] + grp + ['\\end{tabular}', '\\end{ruledtabular}', '\\end{table*}', '']
    open(os.path.join(out, 'C_table.tex'), 'w', encoding='utf-8').write('\n'.join(L) + '\n')
    P = ['% generated by src/ei_tables.py from results_EH/eh_constants.json (the strings of the record set as formulas by PRINTED) -- do not edit by hand',
         '\\begin{table}[t]', '\\caption{Rounded values printed in the working record and used there as one-sided bounds, checked for their direction against the enclosure of the exact value. '
         'The three marked ``wrong\'\' are corrected in this paper; the constants derived from them hold when recomputed (Table~\\ref{tab:constants}). '
         'The last two were used as lower bounds. $^\\dagger$The original form of Lemma EE-C, not used by Theorem~\\ref{thm:main}.}',
         '\\label{tab:printed}', '\\begin{ruledtabular}', '\\footnotesize\\setlength{\\tabcolsep}{2pt}', '\\begin{tabular}{@{}llll@{}}', 'record & printed value & enclosure & check\\\\', '\\hline']
    for c in eh['printed']:
        place, val = PRINTED[c['what']]          # a record string without a formula in PRINTED stops the script (KeyError)
        P.append('%s & %s & %s & %s\\\\' % (place, val, fmt(c['enclosure'][0]), 'right' if c['ok'] else '\\textbf{wrong}'))
    P += ['\\end{tabular}', '\\end{ruledtabular}', '\\end{table}']
    open(os.path.join(out, 'C_printed.tex'), 'w', encoding='utf-8').write('\n'.join(P) + '\n')
    JL = {'[4]': 'P.\\ref{prop:branch}', '[5]->[6]': 'L.\\ref{lem:G2}$\\to$L.\\ref{lem:EF1}', '[6]->[7]': 'L.\\ref{lem:EF1}$\\to$L.\\ref{lem:EF2}', '[7]->[23]': 'L.\\ref{lem:EF2}$\\to$P.\\ref{prop:outside}',
          '[4](b)': 'P.\\ref{prop:branch}', '[9]': 'L.\\ref{lem:K}', '[11]': 'L.\\ref{lem:K}', '[10]': 'L.\\ref{lem:K}', '[12,18]->[23]': 'L.\\ref{lem:K},T.\\ref{thm:circles}$\\to$P.\\ref{prop:outside}',
          '[18]->[23]': 'T.\\ref{thm:circles}$\\to$P.\\ref{prop:outside}', '[12]': 'L.\\ref{lem:K}', '[8]': 'L.\\ref{lem:CR}', '[8]->[28]': 'L.\\ref{lem:CR}$\\to$T.\\ref{thm:main}',
          '[8,18]->[22,23]': 'L.\\ref{lem:CR},T.\\ref{thm:circles}$\\to$L.\\ref{lem:EFO},P.\\ref{prop:outside}', '[20]': 'L.\\ref{lem:EEC}', '[19]->[21]': 'L.\\ref{lem:window}', '[19]': 'L.\\ref{lem:window} (wrong version)',
          '[22]': 'L.\\ref{lem:EFO}', '[22,8,18]->[28]': 'L.\\ref{lem:EFO}$\\to$T.\\ref{thm:main}', '[23]': 'P.\\ref{prop:outside}', '[23]->[28]': 'P.\\ref{prop:outside}$\\to$T.\\ref{thm:main}',
          '[20]->[23]': 'L.\\ref{lem:EEC}$\\to$P.\\ref{prop:outside}', '[24]': 'L.\\ref{lem:EGD}', '[25,26]->[27]': 'L.\\ref{lem:iso},\\ref{lem:sm}$\\to$T.\\ref{thm:EGO}', '[27]': 'T.\\ref{thm:EGO}',
          '[14-18]': 'L.\\ref{lem:kernel}--T.\\ref{thm:circles}'}
    D = ['% generated by src/ei_tables.py from results_EH/eh_joints.json -- do not edit by hand (one table* float)', '\\begin{table*}[p]', '\\scriptsize',
         '\\caption{The joints of the chain on recorded states (\\code{src/eh\\_joints.py}). Each row inserts the conclusion of one part into the hypothesis of the next, or checks one part, on the states listed; the margin is in the unit stated in the row (non-negative means the inequality holds). J18 is a deliberately wrong version (the window inequality without $\\sigma^2$) and must fail. Two margins sit at the edge: J16 is $0$ on a state with $\\sigma=10^{-36}$ where both sides equal $1$ in double precision, and J17 is $-1.2\\times10^{-7}$, inside the quadrature tolerance $2\\times10^{-5}$ of the engine (Paper~F\'s dictionary states are equality cases of the window inequality). The expressions are printed as evaluated by the code (\\code{src/eh\\_joints.py}).}\\label{tab:joints}',
         '\\begin{ruledtabular}', '\\begin{tabular}{@{}llp{0.52\\textwidth}rl@{}}', '& parts & what is checked & states & smallest margin\\\\', '\\hline']
    for j in jo['joints']:
        res = 'breaks (expected)' if j.get('expected_to_break') else fmt(j['min_margin'])
        D.append('%s & %s & \\texttt{%s} & %d & %s\\\\' % (j['J'], JL[j['items']], tex(j['what']), j['n'], res))
    D += ['\\end{tabular}', '\\end{ruledtabular}', '\\end{table*}']
    open(os.path.join(out, 'D_table.tex'), 'w', encoding='utf-8').write('\n'.join(D) + '\n')
    print('tables written:', len(eh['checks']) + len(ei['checks']), 'constants,', len(eh['printed']), 'printed values,', len(jo['joints']), 'joints')


if __name__ == '__main__':
    main(sys.argv[1])
