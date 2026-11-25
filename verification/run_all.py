"""Paper K -- the full verification.  Python 3.11 grammar.
    python verification/run_all.py [--e-repo PATH] [--f-repo PATH] [--full]
  (I)   verification/quick.py: the interval arithmetic of App. C, the joints of App. D, Table 2, the independent point checks, the generated tables,
        the printed numbers of the text (no other repository needed; under a minute).
  (II)  the three circles on the recorded states (Theorem 8, Lemmas 17 and 18): src/eg_circ.py 2d, iso, 1d, 1ds, each rerun from the stored states
        of stages ED and EE and compared with its record (about half an hour in all, most of it the 2921 one-variable states); and the real-time conclusions (ii)-(iv) on the 2449 two-dimensional
        states (src/ef_q2.py 2d).
  (III) with the published repository of paper F (--f-repo PATH or the environment variable F_REPO; resolution-gap-two-sided, v1.5.0, read only):
        engine states re-evaluated afresh, without cache (src/eh_fresh.py: four states of the open set O; src/eg_recheck.py: five states and the
        theorem rows of stage EG), each compared with the stored record.
        With --full also the 692 engine states of src/eh_states.py afresh (about an hour) and their record.
  (IV)  with the published repository of paper E (--e-repo PATH or E_REPO; resolution-gap-log-law, read only): the gap formula, the covariance form
        of N, the direct bound and the branch identity quoted in Propositions 22 and 23 are in its text.
Every rerun is made in a temporary copy and compared with the stored record by verification/common.py (the comparison of the published common.py of
papers I and J); the stored records are never overwritten.  Parts (III) and (IV) are skipped, with a note, when the repository is not given."""
import os, subprocess, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import ROOT, check, finish, workdir, rerun_and_compare, f_repo, e_repo
import shutil


def main():
    t0 = time.time()
    q = subprocess.run([sys.executable, os.path.join(ROOT, 'verification', 'quick.py')], capture_output=True, text=True)
    print(q.stdout.rstrip())
    check('(I) verification/quick.py passes', q.returncode == 0, q.stderr[-400:] if q.returncode else '')
    os.environ['EG_NOCACHE'] = '1'; os.environ['EF_FORCE'] = '1'; os.environ['PYTHONHASHSEED'] = '0'
    w = workdir()
    try:
        for part, what in (('2d', 'two-dimensional states with closed-form fields in complex time (1441 with the circles, 2449 in real time)'),
                           ('iso', 'the isotropic flow (Lemmas 17, 18) on the 1441 two-dimensional states'),
                           ('1d', 'the 2921 one-variable states q = (m mu) * N'),
                           ('1ds', 'a finer quadrature on 294 of the one-variable states')):
            rerun_and_compare(w, 'eg_circ.py', 'results_EG/eg_circ_%s.json' % part, '(II) Theorem 8 step by step: %s' % what, args=[part])
        rerun_and_compare(w, 'ef_q2.py', 'results_EF/ef_q2_2d.json', '(II) Theorem 8 (ii)-(iv) and rho* on the 2449 two-dimensional states', args=['2d'])
        fr = f_repo()
        if fr:
            os.environ['F_REPO'] = fr
            rerun_and_compare(w, 'eh_fresh.py', 'results_EH/eh_fresh.json', '(III) four engine states of the open set re-evaluated afresh with paper F\'s engine (no cache)', frepo=fr)
            rerun_and_compare(w, 'eg_recheck.py', 'results_EG/eg_recheck.json', '(III) five engine states and the theorem rows of stage EG afresh (no cache)', frepo=fr)
            if '--full' in sys.argv:
                rerun_and_compare(w, 'eh_states.py', 'results_EH/eh_states.json', '(III) the 692 engine states afresh (about an hour)', frepo=fr)
        else:
            print('[NOTE] (III) skipped: give the published repository of paper F with --f-repo PATH or F_REPO')
        er = e_repo()
        if er:
            try:
                t2 = open(os.path.join(er, 'paper', 'sections', '02_setting.tex'), encoding='utf-8').read()
                t3 = open(os.path.join(er, 'paper', 'sections', '03_direct.tex'), encoding='utf-8').read()
                okE = (all(x in t2 for x in ('Gap formula', r'\sigma g = \frac{2\sigma^2 J_{ww} + N}{d}', r'\Cov(X_p,\, X_z - \sin X_z \mid Y)'))
                       and all(x in t3 for x in ('Direct bound', r'\label{eq:identity}')))
            except OSError:
                okE = False
            check("(IV) paper E's gap formula, the covariance form of N, the direct bound and the branch identity are in its text (Propositions 22, 23)", okE)
        else:
            print('[NOTE] (IV) skipped: give the published repository of paper E with --e-repo PATH or E_REPO')
    finally:
        shutil.rmtree(w, ignore_errors=True)
    print('run_all.py: %.0f s' % (time.time() - t0))
    finish()


if __name__ == '__main__':
    main()
