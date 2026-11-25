"""Shared helpers of verification/quick.py and verification/run_all.py (paper K).  Python 3.11 grammar.  The comparison of records is that of the
published common.py of papers I and J, unchanged (strip, diff, FLOAT_RTOL, FLOAT_ATOL).
Every rerun is made in a temporary work folder that holds a copy of the repository (src/, results_*/, verification/, paper/sections; not .git); the
record it writes there is compared with the stored record, entry by entry, ignoring only timing and environment entries.  The stored records are never
overwritten.  Booleans, integers and strings (among them every interval-arithmetic bound kept as a decimal string, and every 'ok' entry) are compared
exactly.  Floats from double-precision numerical steps may differ in the last digits between machines (CPU, libm, BLAS); they are compared with a
relative tolerance FLOAT_RTOL, and two floats that are both at the rounding level (|x| <= FLOAT_ATOL) are taken as equal."""
import json, os, shutil, subprocess, sys, tempfile, time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
IGNORE = ('second', 'python', 'numpy', 'platform', 'elapsed', 'time')
RES = []
RESULTS = ('results_ED', 'results_EE', 'results_EF', 'results_EG', 'results_EH', 'results_EI')
FLOAT_RTOL = 1e-9
FLOAT_ATOL = 1e-6


def check(name, ok, info=''):
    RES.append(bool(ok)); print(('[PASS] ' if ok else '[FAIL] ') + name + (('  ' + info) if info else ''), flush=True)


def finish():
    n, p = len(RES), sum(RES)
    print('\n%d/%d checks passed' % (p, n)); sys.exit(0 if p == n else 1)


def workdir():
    w = tempfile.mkdtemp(prefix='paperK_verify_')
    for d in ('src',) + RESULTS:
        shutil.copytree(os.path.join(ROOT, d), os.path.join(w, d), ignore=shutil.ignore_patterns('__pycache__'))
    shutil.copytree(os.path.join(ROOT, 'verification'), os.path.join(w, 'verification'), ignore=shutil.ignore_patterns('__pycache__'))
    shutil.copytree(os.path.join(ROOT, 'paper', 'sections'), os.path.join(w, 'paper', 'sections'))
    return w


def _repo(flag, env, probe):
    p = None
    if flag in sys.argv: p = sys.argv[sys.argv.index(flag) + 1]
    p = p or os.environ.get(env)
    if not p: return None
    p = os.path.abspath(p)
    return p if os.path.exists(os.path.join(p, probe)) else None


def f_repo():
    """the published repository of paper F (resolution-gap-two-sided), from --f-repo PATH or the environment variable F_REPO; None if not given.
    Its at_tools.py is the engine of the 692 recorded engine states (Theorem G1 of paper F, the window inequality, k_eff)."""
    return _repo('--f-repo', 'F_REPO', os.path.join('src', 'at_tools.py'))


def e_repo():
    """the published repository of paper E, from --e-repo PATH or the environment variable E_REPO; None if not given."""
    return _repo('--e-repo', 'E_REPO', 'README.md')


def run(w, script, args=(), frepo=None, timeout=7200):
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1')     # nothing is written into the published repository of paper F
    if frepo: env['F_REPO'] = frepo
    t0 = time.time()
    p = subprocess.run([sys.executable, '-W', 'ignore', script if os.sep in script else os.path.join('src', script)] + list(args), cwd=w, env=env, capture_output=True, text=True, timeout=timeout)
    return p.returncode, p.stdout + p.stderr, time.time() - t0


def strip(o):
    if isinstance(o, dict): return {k: strip(v) for k, v in o.items() if not any(s in k.lower() for s in IGNORE)}
    if isinstance(o, list): return [strip(v) for v in o]
    return o


def diff(a, b, path=''):
    """first difference between two records (None if equal)."""
    if isinstance(a, dict) and isinstance(b, dict):
        if set(a) != set(b): return path + ': keys ' + str(sorted(set(a) ^ set(b))[:3])
        for k in a:
            d = diff(a[k], b[k], path + '/' + k)
            if d: return d
        return None
    if isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b): return path + ': lengths %d, %d' % (len(a), len(b))
        for i, (x, y) in enumerate(zip(a, b)):
            d = diff(x, y, '%s[%d]' % (path, i))
            if d: return d
        return None
    if isinstance(a, float) and isinstance(b, float) and not isinstance(a, bool):
        if a == b or (abs(a) <= FLOAT_ATOL and abs(b) <= FLOAT_ATOL) or abs(a - b) <= FLOAT_RTOL * max(abs(a), abs(b)):
            return None
    return None if a == b else '%s: %r vs %r' % (path, a, b)


def same_record(w, name, only=None):
    """name: 'results_EH/eh_joints.json' etc."""
    a = strip(json.load(open(os.path.join(w, name), encoding='utf-8')))
    b = strip(json.load(open(os.path.join(ROOT, name), encoding='utf-8')))
    if only: a = {k: v for k, v in a.items() if k in only}; b = {k: v for k, v in b.items() if k in only}
    return diff(a, b)


def rerun_and_compare(w, script, record, label, args=(), frepo=None, only=None):
    """only: compare these top-level entries of the record (the others must still say ok, if they have an 'ok' entry)."""
    rec = record
    if os.path.exists(os.path.join(w, rec)): os.remove(os.path.join(w, rec))     # the rerun must write the record itself (some scripts skip an existing one)
    rc, out, sec = run(w, script, args, frepo)
    d = same_record(w, rec, only) if rc == 0 and os.path.exists(os.path.join(w, rec)) else 'no record written'
    if d is None:
        new = json.load(open(os.path.join(w, rec), encoding='utf-8'))
        if not new.get('ok', True): d = 'the rerun record says ok = false'
    check('%s (%s, %.0f s): exit %d, the rerun record equals %s%s' % (label, script, sec, rc, rec, '' if not d else ' -- first difference: ' + d[:200]),
          rc == 0 and d is None, '' if rc == 0 else out[-600:])
