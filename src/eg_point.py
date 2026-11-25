"""Stage EG: independent checks on the main side's three point masses (points (0, 0), (1.5, 1), (-1, 2.2) in (v, w), weights 0.5, 0.3, 0.2).
Everything from the atoms directly (complex Gaussian kernels, posterior weights); no code shared with eg_circ.py.
 (a) the main side's numbers (sigma = 1): F(1 - s)/s = 3.4805, 1.4732, 1.1138 at s = 0.2, 0.5, s* = 0.7223; rho* = 1.1816; Q2 = 0.525690
 (b) the circles (any T > 1; here T = 1.25, delta = 0.25, eta = 0.3, r1 = delta (1 - eta), r3 = T/2, 16 angles):
       N(z) = int |q_z|^2/p_T <= 1/sqrt(1 - (r3/T)^2) on the outer circle, <= A/sqrt(1 - (r1/delta)^2) on the inner circle,
       N(1 - s) <= (sup inner)^(1 - theta) (sup outer)^theta (three circles, with the sampled sups and with the bounds)
     and the isotropic form (both variances z): outer <= T^3/(T^2 - r3^2)^2, inner <= A/(1 - (r1/delta)^2)
 (c) p_{1-s} = sqrt(T/(1-s)) p_T W and q_{1-s} = -sqrt(T/(1-s)) p_T E[Z_w e^{-alpha Z_v^2} | Y_T], with W and the conditional mean from the posterior weights at noise T
 (d) Lemma EG-D (sigma = 1, 0.5, 0.25, 0.1; z = (v + w)/sqrt2, p = (v - w)/sqrt2, g(x) = x - sin x, Z = (Y - X)/sigma):
       N_def = -sigma^-2 E[Cov(X_p, g(X_z) | Y)]  =  N_dir = E[s_p(Y) (g(X_z) - g(Y_z))]
       N = E[s~_p s~_z phi0] + (sigma/2) E[s~_p sin(Y_z) W_zz] + R',  |R'| <= (sqrt15/6) sigma^2 sqrt(J_pp),  R' ~ -(sigma^2/6) E[s~_p cos(Y_z) E[Z_z^3 | Y]]
       2A + N = 2 A_w0 + N_phi0 + E2 + R'
 (e) E[W_zz^2] <= G_z(1 - s)/s  (smoothing in a diagonal component; G_z(t) = E[(d_z ln p_t)^2] with isotropic noise t)
 (f) Lemma EG-D on two wide Gaussian blobs (the state on which the cut v-grid of stage EF's engine was found), sigma = 1e-2 and 1e-4, by a plain
     quadrature over the whole plane (closed forms for the posterior of a Gaussian mixture): the remainder R' against its bound
Output results_EG/eg_point.json.  Python 3.11 grammar."""
import os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import ed_common as C

PTS = np.array([(0.0, 0.0), (1.5, 1.0), (-1.0, 2.2)]); WT = np.array([0.5, 0.3, 0.2])
v = np.linspace(-14.0, 14.5, 1141); w = np.linspace(-11.0, 13.2, 969); V, W = np.meshgrid(v, w, indexing='ij'); DA = (v[1] - v[0]) * (w[1] - w[0])


def ker(u, z): return np.exp(-u * u / (2 * z)) / np.sqrt(2 * np.pi * z)


def pq(z, zw=1.0):
    """p_z and q_z = d_w p_z with v-variance z and w-variance zw (complex allowed)"""
    p = 0; q = 0
    for a, (xv, xw) in zip(WT, PTS):
        k = a * ker(V - xv, z) * ker(W - xw, zw); p = p + k; q = q - (W - xw) / zw * k
    return p, q


def main():
    t0 = time.time(); out = {'env': C.env()}; I = lambda x: float(np.sum(x).real * DA)
    p1, q1 = pq(1.0); p1 = p1.real; q1 = q1.real; A = I(q1 ** 2 / p1); F = lambda t: I(pq(t)[1].real ** 2 / pq(t)[0].real); out['A'] = A; out['mass'] = I(p1)
    # (a)
    ss = 1 / (1 + np.log(1 / A)); out['s_star'] = ss; out['a'] = {'F_over_s': [F(1 - s) / s for s in (0.2, 0.5, ss)], 'main_side': [3.4805, 1.4732, 1.1138], 'rho_star': F(1 - ss) / A, 'rho_star_main_side': 1.1816}
    out['a']['agrees'] = bool(max(abs(x - y) for x, y in zip(out['a']['F_over_s'], out['a']['main_side'])) < 6e-5 and abs(out['a']['rho_star'] - 1.1816) < 6e-5)
    # (b)
    T, dl, eta, s = 1.25, 0.25, 0.3, 0.1; r1, r3, r2 = dl * (1 - eta), T / 2, dl + s; ang = np.arange(16) * 2 * np.pi / 16; th = np.log(r2 / r1) / np.log(r3 / r1); B = {}
    for iso in (False, True):
        zw = (lambda z: z) if iso else (lambda z: 1.0); pT = pq(T, zw(T))[0].real; N = lambda z: I(np.abs(pq(z, zw(z))[1]) ** 2 / pT)
        Aiso = I(pq(1.0, 1.0)[1].real ** 2 / p1)
        o = [N(T + r3 * np.exp(1j * a)) for a in ang]; i_ = [N(T + r1 * np.exp(1j * a)) for a in ang]; Ns = N(1 - s)
        bo = T ** 3 / (T ** 2 - r3 ** 2) ** 2 if iso else 1 / np.sqrt(1 - (r3 / T) ** 2); bi = Aiso / (1 - (r1 / dl) ** 2) if iso else A / np.sqrt(1 - (r1 / dl) ** 2)
        B['iso' if iso else 'v'] = {'outer_sup': max(o), 'outer_bound': bo, 'inner_sup': max(i_), 'inner_bound': bi, 'N_1_minus_s': Ns, 'theta': th,
                                    'three_sampled': max(i_) ** (1 - th) * max(o) ** th, 'three_bounds': bi ** (1 - th) * bo ** th,
                                    'ok': bool(max(o) <= bo and max(i_) <= bi and Ns <= max(i_) ** (1 - th) * max(o) ** th * (1 + 1e-9))}
    out['b'] = B
    # (c)
    al = (s + dl) / (2 * (1 - s)); pT = pq(T)[0].real; ps, qs = pq(1 - s); c = np.sqrt(T / (1 - s))
    Lw = np.stack([np.log(a) - (V - xv) ** 2 / (2 * T) - (W - xw) ** 2 / 2 for a, (xv, xw) in zip(WT, PTS)]); pi = np.exp(Lw - Lw.max(0)); pi /= pi.sum(0)
    Wc = sum(pi[i] * np.exp(-al * (V - PTS[i, 0]) ** 2 / T) for i in range(3)); EZ = sum(pi[i] * (W - PTS[i, 1]) * np.exp(-al * (V - PTS[i, 0]) ** 2 / T) for i in range(3))
    ok = pT > 1e-200
    out['c'] = {'alpha': al, 'max_rel_diff_p': float(np.max(np.abs(ps.real[ok] / (c * pT[ok] * Wc[ok]) - 1))), 'max_abs_diff_q_over_p': float(np.max(np.abs(qs.real[ok] / (c * pT[ok]) + EZ[ok])))}
    # (d)
    D = []
    for sig in (1.0, 0.5, 0.25, 0.1):
        vv = np.linspace(-1 - 9.5 * sig, 1.5 + 9.5 * sig, 1201); wg = np.linspace(-9.5 * sig, 2.2 + 9.5 * sig, 1201); Vg, Wg = np.meshgrid(vv, wg, indexing='ij'); dA = (vv[1] - vv[0]) * (wg[1] - wg[0])
        Lg = np.stack([np.log(a) - ((Vg - xv) ** 2 + (Wg - xw) ** 2) / (2 * sig ** 2) for a, (xv, xw) in zip(WT, PTS)]); mx = Lg.max(0); pi = np.exp(Lg - mx); pp = pi.sum(0) * np.exp(mx) / (2 * np.pi * sig ** 2); pi /= pi.sum(0)
        E = lambda x: float(np.sum(pp * x) * dA); cm = lambda f: sum(pi[i] * f(i) for i in range(3))            # conditional mean given y
        xz = (PTS[:, 0] + PTS[:, 1]) / np.sqrt(2); xp = (PTS[:, 0] - PTS[:, 1]) / np.sqrt(2); yz = (Vg + Wg) / np.sqrt(2); yp = (Vg - Wg) / np.sqrt(2); g = lambda x: x - np.sin(x)
        Zv = lambda i: (Vg - PTS[i, 0]) / sig; Zw = lambda i: (Wg - PTS[i, 1]) / sig; Zz = lambda i: (yz - xz[i]) / sig; Zp = lambda i: (yp - xp[i]) / sig
        sv, sw, sz, sp = -cm(Zv), -cm(Zw), -cm(Zz), -cm(Zp)                                                    # dimensionless scores
        mp = cm(lambda i: xp[i] + 0 * Vg); mg = cm(lambda i: g(xz[i]) + 0 * Vg)
        N_def = -E(cm(lambda i: (xp[i] - mp) * (g(xz[i]) - mg))) / sig ** 2; N_dir = E(sp * cm(lambda i: g(xz[i]) - g(yz))) / sig
        ph0 = 1 - np.cos(yz); Wzz = cm(lambda i: Zz(i) ** 2 - 1); t1 = E(sp * sz * ph0); t2 = sig / 2 * E(sp * np.sin(yz) * Wzz); R = N_dir - t1 - t2
        Jpp = E(sp * sp); bound = np.sqrt(15) / 6 * sig ** 2 * np.sqrt(Jpp); pred = -sig ** 2 / 6 * E(sp * np.cos(yz) * cm(lambda i: Zz(i) ** 3)); A_ = E(sw * sw); B_ = E(sv * sv)
        A0 = E(sw * sw * (1 - ph0 / 4)); N0 = 0.5 * E(sv * sv * ph0)
        D.append({'sigma': sig, 'N_def': N_def, 'N_dir': N_dir, 'term1': t1, 'term2_E2': t2, 'R': R, 'bound': bound, 'R_leading': pred, 'fourth_order_bound': sig ** 3 / 24 * np.sqrt(105 * Jpp),
                  'A': A_, 'B': B_, 'two_A_plus_N': 2 * A_ + N_dir, 'two_Aw0_plus_Nphi0_plus_E2_plus_R': 2 * A0 + N0 + t2 + R, 'gap_g': (2 * A_ + N_dir) / (sig * (A_ + B_)), 'mass': E(1.0),
                  'ok': bool(abs(N_def - N_dir) < 1e-9 and abs(R) <= bound and abs(R - pred) <= sig ** 3 / 24 * np.sqrt(105 * Jpp) + 1e-10 and abs(2 * A_ + N_dir - (2 * A0 + N0 + t2 + R)) < 1e-10)})
    out['d'] = D
    # (e) sigma = 1, isotropic noise: W_zz = d_z^2 p/p, G_z(t) = E[(d_z ln p_t)^2]
    def fz(t):
        p = 0; pz = 0; pzz = 0
        for a, (xv, xw) in zip(WT, PTS):
            k = a * ker(V - xv, t) * ker(W - xw, t); u = ((V - xv) + (W - xw)) / np.sqrt(2); p = p + k; pz = pz - u / t * k; pzz = pzz + (u * u / t ** 2 - 1 / t) * k
        return p, pz, pzz
    p, pz, pzz = fz(1.0); W2 = I(pzz ** 2 / p); Ee = []
    for s_ in (0.1, 0.3, 0.6):
        pa, pza, _ = fz(1 - s_); Ee.append({'s': s_, 'W_zz_L2sq': W2, 'G_over_s': I(pza ** 2 / pa) / s_})
    out['e'] = Ee
    # (f) two Gaussian blobs: (weight, centre (v, w), variances (v, w))
    BL = [(0.5255039170606415, (0.0, 0.84940575), (0.26482592, 2.22267838)), (0.4744960829393585, (-3.0, -0.10616798), (0.26482592, 0.0143413))]; Ff = []
    for sig in (1e-2, 1e-4):
        s2 = sig * sig; n = 2001; vv = np.linspace(-3 - 8 * np.sqrt(0.2649), 8 * np.sqrt(0.2649), n); wg = np.linspace(0.8494 - 8 * np.sqrt(2.2228), 0.8494 + 8 * np.sqrt(2.2228), n)
        Vg, Wg = np.meshgrid(vv, wg, indexing='ij'); dA = (vv[1] - vv[0]) * (wg[1] - wg[0]); Pp = 0; Pv = 0; Pw = 0; Dn = 0; Pzz = 0
        for be, (mv_, mw_), (Sv, Sw) in BL:
            Tv, Tw = Sv + s2, Sw + s2; u_, x_ = Vg - mv_, Wg - mw_; k = be * np.exp(-0.5 * (u_ * u_ / Tv + x_ * x_ / Tw)) / (2 * np.pi * np.sqrt(Tv * Tw))
            mz = (mv_ + u_ * Sv / Tv + mw_ + x_ * Sw / Tw) / np.sqrt(2); vz = (Sv * s2 / Tv + Sw * s2 / Tw) / 2
            Pp = Pp + k; Pv = Pv - u_ / Tv * k; Pw = Pw - x_ / Tw * k; Dn = Dn + k * (mz - np.sin(mz) * np.exp(-vz / 2)); Pzz = Pzz + k * ((u_ / Tv + x_ / Tw) ** 2 / 2 - (1 / Tv + 1 / Tw) / 2)
        E_ = lambda f_: float(np.sum(Pp * f_) * dA); sv, sw = Pv / Pp, Pw / Pp; sp = (sv - sw) / np.sqrt(2); yz = (Vg + Wg) / np.sqrt(2); ph0 = 1 - np.cos(yz)
        A_ = s2 * E_(sw * sw); B_ = s2 * E_(sv * sv); N_ = E_(sp * (Dn / Pp - (yz - np.sin(yz)))); A0 = s2 * E_(sw * sw * (1 - ph0 / 4)); N0 = 0.5 * s2 * E_(sv * sv * ph0)
        E2 = 0.5 * s2 * E_(sp * np.sin(yz) * s2 * Pzz / Pp); bound = np.sqrt(15) / 6 * s2 * np.sqrt(s2 * E_(sp * sp)); R_ = 2 * A_ + N_ - 2 * A0 - N0 - E2
        Ff.append({'sigma': sig, 'A': A_, 'B': B_, 'sigma_AB_g': 2 * A_ + N_, 'main': 2 * A0 + N0, 'E2': E2, 'R': R_, 'bound': bound, 'mass': E_(1.0), 'ok': bool(abs(R_) <= bound and abs(E_(1.0) - 1) < 1e-9)})
    out['f'] = Ff
    out['ok'] = bool(all(x['ok'] for x in Ff) and out['a']['agrees'] and all(x['ok'] for x in B.values()) and out['c']['max_rel_diff_p'] < 1e-9 and out['c']['max_abs_diff_q_over_p'] < 1e-9 and all(x['ok'] for x in D)
                     and all(x['W_zz_L2sq'] <= x['G_over_s'] for x in Ee) and abs(out['mass'] - 1) < 1e-9)
    out['seconds'] = time.time() - t0; C.save('eg_point.json', out, C.RES_EG)
    print('(a) F(1-s)/s:', ['%.4f' % x for x in out['a']['F_over_s']], 'main side', out['a']['main_side'], ' rho* %.4f (1.1816)  agrees %s' % (out['a']['rho_star'], out['a']['agrees']))
    for k, x in B.items(): print('(b) %-3s outer sup %.4f <= %.4f; inner sup %.4f <= %.4f; N(1-s) %.4f <= %.4f (sampled sups) <= %.4f (bounds); theta %.3f  ok %s' % (
        k, x['outer_sup'], x['outer_bound'], x['inner_sup'], x['inner_bound'], x['N_1_minus_s'], x['three_sampled'], x['three_bounds'], x['theta'], x['ok']))
    print('(c) identity p: max rel diff %.1e; q: max abs diff %.1e' % (out['c']['max_rel_diff_p'], out['c']['max_abs_diff_q_over_p']))
    for x in D: print('(d) sigma %.2f: N_def %.8f N_dir %.8f | term1 %.6f E2 %.6f R %.3e (bound %.3e, leading %.3e, 4th-order %.1e) | g %.5f  ok %s' % (
        x['sigma'], x['N_def'], x['N_dir'], x['term1'], x['term2_E2'], x['R'], x['bound'], x['R_leading'], x['fourth_order_bound'], x['gap_g'], x['ok']))
    for x in Ee: print('(e) s %.1f: E[W_zz^2] %.5f <= G_z(1-s)/s %.5f' % (x['s'], x['W_zz_L2sq'], x['G_over_s']))
    for x in Ff: print('(f) two Gaussian blobs, sigma %g: sigma (A+B) g = %.9e, 2 A_w0 + N_phi0 = %.9e, E2 %.3e, R\' %.3e (bound %.3e)  ok %s' % (x['sigma'], x['sigma_AB_g'], x['main'], x['E2'], x['R'], x['bound'], x['ok']))
    print('ok', out['ok'])


if __name__ == '__main__':
    main()
