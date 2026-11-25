"""Stage EB: a general state engine on top of paper F's exact point fields (v3_tools / at_tools of resolution-gap-two-sided v1.5, unchanged).
A state is
    rho(v, w) = sum_b beta_b N(v; 0, Lv^2) N(w; wc_b, Lw_b^2) * [ 1 + sum_j c_j E_j(v, w) cos((kv_j v + kw_j w)/sigma + ph_j) ],
    E_j = exp(-(v - v0_j)^2/(2 Rv_j^2) - (w - w0_j)^2/(2 ell_j^2))   (Rv_j = None: no v-envelope),
lengths in v in units of sigma, in w in units of sqrt(sigma) (as paper F's terms_of); sum_j c_j <= 1 gives rho >= 0.
Each (background x ripple) product is a Gaussian times a plane wave, so paper F's v3_tools evaluates p, the score, the posterior mean and the
curvature term exactly (blur in closed form).  Paper F's own windows (terms_of) are the special case kw = 0, Rv = None, kv = m kb.
Quantities (one pass over a grid):
  paper F: g, theta = (1 - delta_c) Phi, C_w, k_eff^2, eps, the bounds of Lemma BB1 and Theorem AT2, A = sigma^2 J_ww, B = sigma^2 J_vv;
           W2 exact from the moments of rho (and paper F's background-only W2_of for comparison); Q = 4 theta/W2; case 3 as bd_case3.row.
  CM-1 (paper H, Sec. VIII A), with p = p_bg (1 + r), p_bg the blurred background alone, v-averages over the focus (the observed v-marginal):
           on the line w* where the ripple energy is largest: E0 = E_v[r^2], x E1 = sigma^2 E_v[(d_v r)^2], mu = x E1/E0, E0' (focus twice as wide),
           b = sigma^2/V (V the observed v-variance), B0, the right side of CM-1;
           over the whole state: E0hat = E_pbg[r^2], muhat; the ripple's share of the w-information u = J_rip/J_bg (J_rip = int p_bg (d_w r)^2,
           J_bg = int (d_w p_bg)^2/p_bg), its w-weighted fineness mu_w = sigma^2 int p_bg (d_v d_w r)^2 / J_rip, the scale ratio Lam_s = u/E0.
Python 3.11 grammar.
Stage EE copy (changes, marked "EE:"; the ED records are unaffected, since none of their states has a background off v = 0 or a density that underflows):
  a background may carry a v-centre: bg entries (beta, wc, Lw) or (beta, wc, Lw, vc), vc in sigma units; the v-nodes then cover every v-centre;
  p and p_bg are floored at 1e-300 in the extra sums (far-apart blobs: the density underflows between them);
  extra sums q0 = int p^2, qz = int p^2 y_z^2 (paper E's (H1): E_Q[y_z^2] with Q = p^2/||p||^2);
  the log-normalisation of a (background x ripple window) product is written without cancellation, (mb - mj).Sbi S P.(mb - mj) in place of
  mb.Sbi.mb + mj.P.mj - m.S^-1.m (the two agree to rounding; the old form loses all digits when the window sits far from w = 0, e.g. at w = 2 pi sqrt2
  for sigma = 1e-24, and about 0.07 in the log at sigma = 1e-12; for windows near w = 0, as in stages EB-ED, nothing changes beyond rounding)."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import ed_common as C

C.use_F()
import at_tools, v3_tools                     # noqa: E402  (paper F)
from u_cgm import trap_weights               # noqa: E402

KT = ('svv', 'sww', 'Nibp', 'N1', 'phi_svv', 'sww_w', 'ibp_rest', 'phi', 'Cw_int', 'X_int')


# ----------------------------------------------------------------------------------------------------------------- terms
def terms(sig, bg, rip, Lv_sig):
    """bg: list of (beta, wc, Lw) [wc, Lw in sqrt(sigma)]; rip: list of dicts c, v0 [sigma], w0 [sqrt sigma], Rv [sigma or None], ell [sqrt sigma],
    kv, kw [1/sigma], ph.  Returns (T_all, T_bg) in v3_tools' native form (log_coef, centre, S, kappa)."""
    s = np.sqrt(sig); Lv = Lv_sig * sig; bsum = sum(b[0] for b in bg); Tb = []
    for b_ in bg:                                                       # EE: optional v-centre
        beta, wc, Lw = b_[:3]; vc = b_[3] if len(b_) > 3 else 0.0
        Tb.append((np.log(beta / bsum), np.array([vc * sig, wc * s]), np.diag([Lv ** 2, (Lw * s) ** 2]), np.zeros(2)))
    T = list(Tb)
    for b_ in bg:
        beta, wc, Lw = b_[:3]; vc = b_[3] if len(b_) > 3 else 0.0
        Sb = np.diag([Lv ** 2, (Lw * s) ** 2]); mb = np.array([vc * sig, wc * s]); Sbi = np.linalg.inv(Sb)
        for r in rip:
            if r['c'] <= 0: continue
            Rv = r.get('Rv'); P = np.diag([0.0 if Rv is None else 1 / (Rv * sig) ** 2, 1 / (r['ell'] * s) ** 2])
            mj = np.array([r.get('v0', 0.0) * sig, r.get('w0', 0.0) * s])
            S = np.linalg.inv(Sbi + P); m = S @ (Sbi @ mb + P @ mj)
            dmj = mb - mj                                                  # EE: cancellation-free form of the same exponent,
            la = 0.5 * np.log(np.linalg.det(S) / np.linalg.det(Sb)) - 0.5 * (dmj @ Sbi @ S @ P @ dmj)   # mb.Sbi.mb + mj.P.mj - m.S^-1.m = (mb - mj).Sbi S P.(mb - mj)
            kap = np.array([r['kv'], r.get('kw', 0.0)]) / sig; ph = r.get('ph', 0.0) + kap @ m
            for sg in (1, -1):
                T.append((np.log(beta / bsum * r['c'] / 2) + la + 1j * sg * ph, m, S, sg * kap))
    return T, Tb


# ----------------------------------------------------------------------------------------------------------------- p and derivatives
def fields(pre, X):
    """p, grad p, Hessian of p (real) at points X (n x 2) from v3_tools' pre-processed terms (plain sum, used for r = p/p_bg - 1)."""
    p = np.zeros(len(X)); g = np.zeros((len(X), 2)); H = np.zeros((len(X), 2, 2))
    for c0, mu, Ti, K, varz, Kk, Ck, kCk, Sref, S in pre:
        d = X - mu; ex = c0 - 0.5 * np.einsum('ni,ij,nj->n', d, Ti, d) + 1j * (d @ Kk) - 0.5 * kCk
        f = np.exp(ex); gr = -(d @ Ti.T) + 1j * Kk[None, :]
        p += f.real; g += (f[:, None] * gr).real
        H += (f[:, None, None] * (gr[:, :, None] * gr[:, None, :] - Ti[None, :, :])).real
    return p, g, H


# ----------------------------------------------------------------------------------------------------------------- moments of rho (W2)
def W2_exact(sig, T):
    """W2 = E[y_w^2] + (E[x_w^2 z2^2] + sigma^2 E[z2^2])/8, z2 = x_v + x_w, x ~ rho (exact: each term is a Gaussian times a plane wave,
    int N(x; m, S) e^{i k.(x - m)} f(x) dx = e^{-k'Sk/2} E_{N(m + iSk, S)}[f])."""
    M0 = M2w = M2z = M4 = 0.0
    for lc, m, S, kap in T:
        cf = np.exp(lc) * np.exp(-0.5 * kap @ S @ kap)
        mu = m + 1j * (S @ kap)
        a = mu[1]; bz = mu[0] + mu[1]; va = S[1, 1]; vb = S[0, 0] + 2 * S[0, 1] + S[1, 1]; cab = S[0, 1] + S[1, 1]
        M0 += cf; M2w += cf * (va + a * a); M2z += cf * (vb + bz * bz)
        M4 += cf * ((va + a * a) * (vb + bz * bz) + 2 * cab ** 2 + 4 * a * bz * cab)
    M0 = M0.real; Ew2 = (M2w / M0).real; Ez2 = (M2z / M0).real; E4 = (M4 / M0).real
    return Ew2 + sig ** 2 + (E4 + sig ** 2 * Ez2) / 8


# ----------------------------------------------------------------------------------------------------------------- grids
def w_grid(sig, bg, rip, nper=40, coarse=300, kw_res=16):
    s = np.sqrt(sig)
    lo = min(b_[1] - 6.5 * b_[2] for b_ in bg); hi = max(b_[1] + 6.5 * b_[2] for b_ in bg)            # EE: entries may carry a v-centre
    pts = [np.linspace(lo, hi, 4 * coarse + 1)]
    for b_ in bg: pts.append(np.arange(b_[1] - 6.5 * b_[2], b_[1] + 6.5 * b_[2] + 1e-12, b_[2] / coarse))
    for r in rip:
        e = r['ell']; w0 = r.get('w0', 0.0); a, b = max(w0 - 10 * e, lo), min(w0 + 10 * e, hi)
        if b <= a: continue
        h = e / nper
        if abs(r.get('kw', 0.0)) > 0: h = min(h, 2 * np.pi * s / (kw_res * abs(r['kw'])))     # resolve the w-oscillation (w in sqrt(sigma) units)
        pts.append(np.arange(a, b + 1e-12, h))
    return np.unique(np.concatenate(pts)) * s


def v_grid(sig, Lv_sig, rip, per_res=16, ext=7.0, vcs=(0.0,)):
    """v-grid over +-ext observed v-spreads around every v-centre of the background; spacing <= period/per_res of the finest kv where a ripple lives
    (inside its v-envelope, +-7 Rv, or everywhere for a ripple without v-envelope), Lo/40 elsewhere.
    EG: 'vcs' added.  Before stage EG the grid covered only the centre 0, so a background blob with a v-centre (stage EF's family W, second blob) was cut off."""
    Lo = np.sqrt(Lv_sig ** 2 + 1.0); coarse = Lo / 40; pts = []
    for vc in sorted(set(vcs)):
        lo, hi = vc - ext * Lo, vc + ext * Lo; pts.append(np.arange(lo, hi + 1e-12, coarse))
        for r in rip:
            if abs(r['kv']) == 0 and r.get('Rv') is None: continue
            h = min(2 * np.pi / (per_res * max(abs(r['kv']), 1e-9)), (r['Rv'] / 10 if r.get('Rv') else coarse))
            a, b = (lo, hi) if r.get('Rv') is None else (max(lo, r.get('v0', 0.0) - 7 * r['Rv']), min(hi, r.get('v0', 0.0) + 7 * r['Rv']))
            if b > a: pts.append(np.arange(a, b + 1e-12, h))
    return np.unique(np.concatenate(pts)) * sig


# ----------------------------------------------------------------------------------------------------------------- one pass
def evaluate(sig, bg, rip, Lv_sig, periodic_kb=None, nv=32, n0=21, nper=40, coarse=300, per_res=16, max_entries=2.5e6, wline=None, TT=None,
             local_kmin=None, nb=8, kw_res=16):
    """periodic_kb: base wavenumber of a periodic state (all kv integer multiples of it, no v-envelopes): the v-integral is paper F's
    sums_envelope (one-period averages along the v-envelope); otherwise a plain grid over v.  Returns the dict of sums."""
    if TT is None:
        T, Tb = terms(sig, bg, rip, Lv_sig); wn = w_grid(sig, bg, rip, nper, coarse, kw_res)
    else:
        T, Tb, wn = TT                                              # given terms (paper F's hidden rows), background terms and w-nodes
    pre = v3_tools._pre(T, sig); preb = v3_tools._pre(Tb, sig); ww = trap_weights(wn)
    Lo = np.sqrt((Lv_sig * sig) ** 2 + sig ** 2)
    if periodic_kb is not None:
        per = 2 * np.pi * sig / periodic_kb
        vcs = sorted(set(float(b_[3]) if len(b_) > 3 else 0.0 for b_ in bg))      # EE: one envelope grid per v-centre
        v0s = np.concatenate([vc * sig + np.linspace(-6.5 * Lo, 6.5 * Lo, n0) for vc in vcs]); w0 = np.concatenate([trap_weights(np.linspace(-6.5 * Lo, 6.5 * Lo, n0))] * len(vcs))
        blocks = [(v0 + (np.arange(nv) + 0.5) * per / nv - per / 2, cw * np.full(nv, 1.0 / nv)) for v0, cw in zip(v0s, w0)]
    elif local_kmin is not None:
        # block-local averages (not periodic: chirps, quasi-periodic and self-similar ripples whose v-envelope is much wider than a period):
        # at each envelope node v0 the local mean is taken over a block of nb periods of the lowest wavenumber, resolved at per_res points per
        # period of the highest one.  Exact for a periodic state; for the others the block length is the accuracy parameter (P3: nb doubled).
        kmax = max(abs(r['kv']) for r in rip); Lb = nb * 2 * np.pi * sig / local_kmin; nbp = int(np.ceil(nb * per_res * kmax / local_kmin))
        v0s = np.linspace(-6.5 * Lo, 6.5 * Lo, n0); w0 = trap_weights(v0s)
        blocks = [(v0 + (np.arange(nbp) + 0.5) * Lb / nbp - Lb / 2, cw * np.full(nbp, 1.0 / nbp)) for v0, cw in zip(v0s, w0)]
    else:
        vn = v_grid(sig, Lv_sig, rip, per_res, vcs=[float(b_[3]) if len(b_) > 3 else 0.0 for b_ in bg]); wv = trap_weights(vn); rows = max(1, int(max_entries // max(1, len(pre) * len(wn))))     # EG: every v-centre covered
        blocks = [(vn[a:a + rows], wv[a:a + rows]) for a in range(0, len(vn), rows)]
    acc = {k: 0.0 for k in at_tools.KT}; m = 0.0
    cm = {k: 0.0 for k in ('mb', 'r2', 'dvr2', 'Jrip', 'Jbg', 'dvwr2', 'keff_num', 'keff_den', 'fww_num', 'svv_d', 'svvphi_d', 'tA', 'tB', 'tD', 'tK', 'tS2', 'tSV', 'tm', 'tVmax', 'q0', 'qz', 'eA0', 'eN0', 'eE2', 'eJpp', 'eW2', 'eU')}     # EG: the sums of the direct expansion around Y (report EG, Lemma EG-D)
    nW = len(wn); line = np.zeros((4, nW))          # per w: sum pb, pb r^2, pb (dv r)^2, (and pb at 2x focus is handled below)
    npts = max(256, int(max_entries // max(1, len(pre))))
    for vb, wvb in blocks:
        V, W = np.meshgrid(vb, wn, indexing='ij'); X = np.stack([V.ravel(), W.ravel()], 1); dA = np.outer(wvb, ww).ravel()
        for i0 in range(0, len(X), npts):
            Xc = X[i0:i0 + npts]; dc = dA[i0:i0 + npts]
            m += at_tools._acc(pre, sig, Xc, dc, acc)
            p, gp, Hp = fields(pre, Xc); pb, gb, Hb = fields(preb, Xc)
            p = np.maximum(p, 1e-300); pb = np.maximum(pb, 1e-300)                                    # EE: floor (far-apart blobs)
            yz = (Xc[:, 0] + Xc[:, 1]) / np.sqrt(2); cm['q0'] += np.sum(p * p * dc); cm['qz'] += np.sum(p * p * yz * yz * dc)   # EE: (H1)
            q = p / pb; r = q - 1
            # derivatives of q = p/pb
            sv, sw = gp[:, 0] / p, gp[:, 1] / p; svb, swb = gb[:, 0] / pb, gb[:, 1] / pb
            qv = q * (sv - svb); qw = q * (sw - swb)
            Hvw = Hp[:, 0, 1] / p - sv * sw; Hvwb = Hb[:, 0, 1] / pb - svb * swb      # d_v d_w ln p
            qvw = q * ((sv - svb) * (sw - swb) + (Hvw - Hvwb))
            w_ = pb * dc
            cm['mb'] += w_.sum(); cm['r2'] += np.sum(w_ * r * r); cm['dvr2'] += sig ** 2 * np.sum(w_ * qv * qv)
            cm['Jrip'] += np.sum(w_ * qw * qw); cm['Jbg'] += np.sum(dc * gb[:, 1] ** 2 / pb); cm['dvwr2'] += sig ** 2 * np.sum(w_ * qvw * qvw)
            # k_eff^2 directly from sqrt p (cross-check of paper F's X_int)
            sq = np.sqrt(p); dw_sq = gp[:, 1] / (2 * sq); dvw_sq = Hp[:, 0, 1] / (2 * sq) - gp[:, 0] * gp[:, 1] / (4 * p * sq)
            cm['keff_num'] += sig ** 2 * np.sum(dc * dvw_sq ** 2); cm['keff_den'] += np.sum(dc * dw_sq ** 2)
            mz = ((Xc[:, 0] + sig ** 2 * sv) + (Xc[:, 1] + sig ** 2 * sw)) / np.sqrt(2)      # Tweedie posterior mean, z = (v + w)/sqrt 2
            cm['svv_d'] += np.sum(dc * p * sv * sv); cm['svvphi_d'] += np.sum(dc * p * sv * sv * 2 * np.sin(mz / 2) ** 2)
            # stage ED: the dimensionless sums of T1 (sigma = 1 units): A, B, D = E[(d_v s_w)^2], K = E[(d_v s_w + s_v s_w/2)^2], S2, SV
            pd = p * dc; s2_ = sig * sig; Hvv_ = Hp[:, 0, 0] / p - sv * sv
            cm['tm'] += pd.sum(); cm['tA'] += s2_ * np.sum(pd * sw * sw); cm['tB'] += s2_ * np.sum(pd * sv * sv)
            cm['tD'] += s2_ ** 2 * np.sum(pd * Hvw ** 2); cm['tK'] += s2_ ** 2 * np.sum(pd * (Hvw + 0.5 * sv * sw) ** 2)
            cm['tS2'] += s2_ ** 2 * np.sum(pd * sv * sv * sw * sw); cm['tSV'] += s2_ * np.sum(pd * sw * sw * (1 + s2_ * Hvv_))
            cm['tVmax'] = max(cm['tVmax'], float(np.max(1 + s2_ * Hvv_)))
            # EG: phi0 = 1 - cos y_z (the weight at the observed point, not at the posterior mean), s_p = (s_v - s_w)/sqrt2, W_zz = sigma^2 d_z^2 p/p = E[Z_z^2 - 1 | y]
            ph0 = 2 * np.sin(yz / 2) ** 2; sp_ = (sv - sw) / np.sqrt(2); Wzz = s2_ * (Hp[:, 0, 0] + 2 * Hp[:, 0, 1] + Hp[:, 1, 1]) / (2 * p)
            cm['eA0'] += s2_ * np.sum(pd * sw * sw * (1 - ph0 / 4)); cm['eN0'] += 0.5 * s2_ * np.sum(pd * sv * sv * ph0)
            cm['eE2'] += 0.5 * s2_ * np.sum(pd * sp_ * np.sin(yz) * Wzz); cm['eJpp'] += s2_ * np.sum(pd * sp_ * sp_)
            cm['eW2'] += np.sum(pd * Wzz * Wzz); cm['eU'] += s2_ * np.sum(pd * sp_ * sp_ * np.sin(yz) ** 2)
            dww_sq = Hp[:, 1, 1] / (2 * sq) - gp[:, 1] ** 2 / (4 * p * sq); cm['fww_num'] += np.sum(dc * dww_sq ** 2)   # (d_w^2 sqrt p)^2
            iw = (i0 + np.arange(len(Xc))) % nW
            np.add.at(line[0], iw, pb * dc); np.add.at(line[1], iw, pb * dc * r * r); np.add.at(line[2], iw, pb * dc * sig ** 2 * qv * qv)
    E = {k: v / m for k, v in acc.items()}
    return {'E': E, 'cm': cm, 'line': line, 'wn': wn, 'T': T, 'n_points': int(sum(len(b[0]) for b in blocks) * nW)}


def cm1_rhs(E0, E0p, b):
    B0 = (1 + np.sqrt(E0) + (E0 / 2) ** (1 / (2 + b))) ** 2
    return (1 + b / 2) * E0 * np.log(2 * B0 / E0) + np.sqrt(2) / np.e * b * E0p, B0


def quantities(sig, bg, rip, Lv_sig, label='', **kw):
    """all quantities of one state; kw go to evaluate (periodic_kb, nv, n0, nper, coarse, per_res, TT)."""
    o = evaluate(sig, bg, rip, Lv_sig, **kw); E, cm, line = o['E'], o['cm'], o['line']
    r = at_tools.analyse(E, sig); s2 = sig * sig; A, B = r['sig2Jww'], r['sig2Jvv']
    theta = (1 - r['delta_c']) * r['Phi']; W2 = W2_exact(sig, o['T']); Q = 4 * theta / W2; ln = np.log(1 / sig); keff2 = r['k_eff'] ** 2
    import bc_cases
    W2F = bc_cases.W2_of(sig, bg[0][2] * np.sqrt(sig), Lv_sig * sig, bg[0][1] * np.sqrt(sig)) if len(bg) == 1 else None
    bb1 = (1 - abs(r['eps'])) * theta / (2 * sig)
    case3 = bool(Q * ln < 0.1 and bb1 < r['bound_AT2'] and A / (A + B) < 2 * sig)
    # CM-1 on the line w* of largest ripple energy
    E0w = np.where(line[0] > 0, line[1] / np.maximum(line[0], 1e-300), 0.0); iw = int(np.argmax(E0w))
    E0 = float(E0w[iw]); xE1 = float(line[2][iw] / line[0][iw]); mu = xE1 / E0 if E0 > 0 else float('nan')
    V = (Lv_sig * sig) ** 2 + sig ** 2; b = s2 / V
    E0p = E0                                                    # focus twice as wide: equal for a v-uniform ripple (checked by the 'E0p_rel' entry where available)
    rhs, B0 = cm1_rhs(E0, E0p, b) if E0 > 0 else (float('nan'), float('nan'))
    E0hat = cm['r2'] / cm['mb']; muhat = cm['dvr2'] / cm['r2'] if cm['r2'] > 0 else float('nan')
    u = cm['Jrip'] / cm['Jbg']; mu_w = cm['dvwr2'] / cm['Jrip'] if cm['Jrip'] > 0 else float('nan')
    keff2_direct = cm['keff_num'] / cm['keff_den']
    Fww = cm['fww_num'] / cm['keff_den']
    theta_direct = cm['svvphi_d'] / cm['svv_d']                 # theta from the plain sums of p (independent of paper F's point_fields)                        # int (d_w^2 sqrt p)^2 / int (d_w sqrt p)^2: the squared w-fineness of the w-gradient
    out = {'label': label, 'sigma': sig, 'ln': ln, 'g': float(r['g']), 'cand1b_margin': float(np.log(W2 / s2) - keff2), 'theta': float(theta), 'W2': float(W2), 'W2_paperF_bg': W2F, 'Q': float(Q),
           'k_eff2': float(keff2), 'k_eff2_direct': float(keff2_direct), 'half_ln_1_over_Q': 0.5 * np.log(1 / Q), 'VI_excess': float(keff2 - 0.5 * np.log(1 / Q)),
           'VI_quarter_excess': float(keff2 - 0.25 * np.log(1 / Q)), 'keff2_over_ln_1_over_Q': float(keff2 / np.log(1 / Q)),
           'A_frac': float(A / (A + B)), 'eps': float(r['eps']), 'AT2': float(r['bound_AT2']), 'BB1': float(bb1), 'C_w': float(r['C_w']), 'Q_ln': float(Q * ln),
           'case3': case3,
           'T1': {'A': float(cm['tA'] / cm['tm']), 'B': float(cm['tB'] / cm['tm']), 'D': float(cm['tD'] / cm['tm']), 'K': float(cm['tK'] / cm['tm']),
                  'S2': float(cm['tS2'] / cm['tm']), 'SV': float(cm['tSV'] / cm['tm']), 'max_Var_Zv': float(cm['tVmax']),
                  'k_eff2': float(cm['tK'] / cm['tA']), 'ln_1_over_A': float(np.log(cm['tm'] / cm['tA'])),
                  'T1_margin': float(np.log(cm['tm'] / cm['tA']) - cm['tK'] / cm['tA']), 'T1_half_margin': float(0.5 * np.log(cm['tm'] / cm['tA']) - cm['tK'] / cm['tA'])},
           'cand1_margin': float(0.5 * np.log(W2 / s2) - keff2),        # candidate 1: k_eff^2 <= (1/2) ln(W2/sigma^2) + C  (Conjecture 18 with Q at its floor)
           'cand2_margin': float(0.5 * np.log(W2 * Fww) - keff2),       # candidate 2: k_eff^2 <= (1/2) ln(W2 F_ww) + C  (the w-scale of the w-gradient)
           'W2_Fww': float(W2 * Fww), 'theta_direct': float(theta_direct), 'Q_direct': float(4 * theta_direct / W2),
           'CM1': {'w_star_sqrt_sigma': float(o['wn'][iw] / np.sqrt(sig)), 'E0': E0, 'xE1': xE1, 'mu': mu, 'b': b, 'B0': B0, 'rhs': rhs, 'holds': bool(xE1 <= rhs),
                   'E0hat': float(E0hat), 'muhat': float(muhat), 'u_share': float(u), 'mu_w': float(mu_w), 'Lam_s': float(u / E0) if E0 > 0 else float('nan'),
                   'keff2_pred': float(mu_w * u / (1 + u)), 'keff2_over_pred': float(keff2 / (mu_w * u / (1 + u))) if u > 0 else float('nan')},
           'n_points': o['n_points']}
    c = out['CM1']
    c['Q_over_E0'] = Q / E0 if E0 > 0 else float('nan'); c['keff2_over_mu'] = keff2 / mu if mu > 0 else float('nan')
    c['lnQ_minus_lnE0'] = float(np.log(1 / Q) - np.log(1 / E0)) if E0 > 0 else float('nan')
    c['lnLam_over_ln_1_over_Q'] = float(np.log(c['Lam_s']) / np.log(1 / Q)) if c['Lam_s'] > 0 else float('nan')
    # the CM-1 x share reading: k_eff^2 <~ max over E of [(1+b/2) ln(2 B0(E)/E)] * E Lam/(1 + E Lam)  (mu bounded by CM-1 at brightness E)
    vals = []
    for Ee in np.logspace(-300, 0, 3001):
        rr, BB = cm1_rhs(Ee, Ee, b); vals.append(rr / Ee * Ee * c['Lam_s'] / (1 + Ee * c['Lam_s']))
    c['cm1_share_bound'] = float(max(vals)); c['cm1_share_bound_over_ln_1_over_Q'] = float(max(vals) / np.log(1 / Q))
    return out
