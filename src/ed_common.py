"""Shared helpers of the chain code (written in stage ED; stages EE-EH added their result folders).  Python 3.11 grammar.
The published repository of paper F (resolution-gap-two-sided, v1.5.0) is used read-only by the engine modules (ed_engine, ee_engine and the
scripts built on them); give its path in the environment variable F_REPO (verification/run_all.py --f-repo PATH sets it).  The older lab names
EI_FREPO, EH_FREPO, ... are still read; without any of them, a sibling folder ../resolution-gap-two-sided is tried.
The comparison of records is that of the published common.py (papers I, J), copied unchanged as diff() below."""
import os, sys, json, platform
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
FREPO = (os.environ.get('F_REPO') or os.environ.get('EI_FREPO') or os.environ.get('EH_FREPO') or os.environ.get('EG_FREPO') or os.environ.get('EF_FREPO')
         or os.environ.get('EE_FREPO') or os.environ.get('ED_FREPO') or os.path.join(os.path.dirname(ROOT), 'resolution-gap-two-sided'))
RES = os.path.join(ROOT, 'results_ED')
RES_EE = os.path.join(ROOT, 'results_EE')
RES_EH = os.path.join(ROOT, 'results_EH')
RES_EI = os.path.join(ROOT, 'results_EI')
RES_EG = os.path.join(ROOT, 'results_EG')
RES_EF = os.path.join(ROOT, 'results_EF')

FLOAT_RTOL = 1e-9
FLOAT_ATOL = 1e-6
IGNORE = ('second', 'python', 'numpy', 'platform', 'elapsed', 'time')


def use_F():
    """put paper F's modules on the path (they chdir to the repository on import; we come back)."""
    cwd = os.getcwd()
    for d in ('verification', 'experiments', 'certificates', 'src'):
        p = os.path.join(FREPO, d)
        if p not in sys.path: sys.path.insert(0, p)
    import bd_case3  # noqa: F401  (its import does os.chdir(FREPO))
    os.chdir(cwd)
    np.seterr(all='ignore')


def strip(o):
    if isinstance(o, dict): return {k: strip(v) for k, v in o.items() if not any(s in k.lower() for s in IGNORE)}
    if isinstance(o, list): return [strip(v) for v in o]
    return o


def diff(a, b, path=''):
    """first difference between two records (None if equal) -- the published common.py, unchanged."""
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


def rel_all(a, b, path=''):
    """stricter companion of diff(): every float compared with the relative tolerance only (no absolute floor), since theta, Q, W2 are tiny."""
    if isinstance(a, dict) and isinstance(b, dict):
        for k in a:
            if k in b:
                d = rel_all(a[k], b[k], path + '/' + k)
                if d: return d
        return None
    if isinstance(a, list) and isinstance(b, list):
        for i, (x, y) in enumerate(zip(a, b)):
            d = rel_all(x, y, '%s[%d]' % (path, i))
            if d: return d
        return None
    if isinstance(a, float) and isinstance(b, float):
        return None if a == b or abs(a - b) <= FLOAT_RTOL * max(abs(a), abs(b)) else '%s: %r vs %r' % (path, a, b)
    return None


def env():
    return {'python': platform.python_version(), 'numpy': np.__version__}


def save(name, obj, res=None):
    res = res or RES
    os.makedirs(res, exist_ok=True)
    fn = os.path.join(res, name)
    json.dump(obj, open(fn, 'w', encoding='utf-8'), indent=1, default=float, ensure_ascii=False)
    return fn
