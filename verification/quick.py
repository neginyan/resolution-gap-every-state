"""Paper K -- the short check (under a minute; it is the check of the GitHub Actions workflow; no other repository needed).  Python 3.11 grammar.
Reruns, in a temporary copy, and compares each new record with the stored one (verification/common.py, the comparison of the published common.py of
papers I and J):
  App. C   every constant of the chain by outward-rounded interval arithmetic (src/eh_constants.py, src/ei_constants.py; 101 claims);
  App. D   the 34 joints of the chain, recomputed from the stored records of the recorded states (src/eh_joints.py);
  Table 2  the range of sigma as a function of the constants (src/eh_range.py);
  Lemmas 4, 6, 16, 18  the independent checks on three point masses and two Gaussian blobs (src/eg_point.py);
  Lemma 16 our recomputation of the three-point check of the review of the manuscript (src/ei_point.py);
  App. C, D the tables of the paper are what the records say (src/ei_tables.py, regenerated and compared with paper/sections);
  the printed numbers of the text against the records and the interval arithmetic (verification/text_numbers.py).
The full verification is verification/run_all.py."""
import filecmp, os, shutil, subprocess, sys, tempfile, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import ROOT, check, finish, workdir, run, rerun_and_compare


def main():
    t0 = time.time(); w = workdir()
    try:
        rerun_and_compare(w, 'eh_constants.py', 'results_EH/eh_constants.json', 'App. C: the constants of the chain by interval arithmetic (72 claims, grids, the printed slips)')
        rerun_and_compare(w, 'ei_constants.py', 'results_EI/ei_constants.json', 'App. C: Theorem 19 at J0 = 1.9e-8, Theorem 1(b) at sigma = 1.5e-14, (c) at 1e-15, delta = 0.4/L (29 claims)')
        rerun_and_compare(w, 'eh_joints.py', 'results_EH/eh_joints.json', 'App. D: the 34 joints of the chain on the recorded states, from the stored records')
        rerun_and_compare(w, 'eh_range.py', 'results_EH/eh_range.json', 'Table 2: the range of sigma as a function of the constants')
        rerun_and_compare(w, 'eg_point.py', 'results_EG/eg_point.json', 'Lemmas 4, 6, 16, 18: three point masses and two Gaussian blobs, from the atoms directly')
        rerun_and_compare(w, 'ei_point.py', 'results_EI/ei_point.json', "Lemma 16: the three-point state of the review of the manuscript (sigma = 0.55), recomputed")
        # the generated tables of the paper
        tdir = tempfile.mkdtemp(prefix='paperK_tables_')
        rc, out, sec = run(w, 'ei_tables.py', [tdir])
        same = rc == 0 and all(filecmp.cmp(os.path.join(tdir, f), os.path.join(ROOT, 'paper', 'sections', f), shallow=False) for f in ('C_table.tex', 'C_printed.tex', 'D_table.tex'))
        check('App. C, D: the tables C_table.tex, C_printed.tex, D_table.tex of the paper equal the ones generated from the records (ei_tables.py, %.0f s)' % sec, same, '' if rc == 0 else out[-400:])
        shutil.rmtree(tdir, ignore_errors=True)
        rerun_and_compare(w, os.path.join('verification', 'text_numbers.py'), 'results_EI/ei_text_numbers.json', 'the printed numbers of the text (every entry, the inventory of the numbers of the running text, the wrong versions)')
    finally:
        shutil.rmtree(w, ignore_errors=True)
    print('quick.py: %.0f s' % (time.time() - t0))
    finish()


if __name__ == '__main__':
    main()
