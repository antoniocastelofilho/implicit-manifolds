import os as _os
REPO = _os.environ.get('IMPLICIT_MANIFOLDS_ROOT') or next(_p for _p in (_os.path.abspath(_os.path.join(_os.path.dirname(__file__), *(['..'] * _i))) for _i in range(8)) if _os.path.exists(_os.path.join(_p, 'CITATION.cff')))
import sys, os, re
sys.path.insert(0, REPO + '/python/implicit_manifolds')
from polfile import read_pol, euler
R = REPO + '/book/runs/hypercubes/'
rows = {}
for line in open(R + 'timings.txt'):
    f, alg, n, k, D, t = line.split()[0], line.split()[1], line.split()[2], line.split()[3], line.split(' D=')[1].rsplit(' ', 1)[0], float(line.split()[-1])
    rows[f] = (alg, n[2:], k[2:], D, t)
cases = [
    ('sphere $S^2\\subset\\R^3$', 'sph3', ['MS', 'MHC', 'GMHC', 'GCHC']),
    ('sphere $S^3\\subset\\R^4$', 'sph4', ['MS', 'GMHC', 'GCHC']),
    ('torus $T^2\\subset\\R^4$', 'tor4', ['MS', 'GCHC']),
    ('$z=\\cos w$', 'zcosw', ['MS', 'GCHC']),
    ('$z^2=\\cos w$', 'z2cosw', ['MS', 'GCHC']),
    ('$z^2=w^5+1+i$', 'z2w5', ['MS', 'GCHC']),
    ('two spheres, $\\min$', 'two_min', ['MS', 'GMHC']),
    ('two spheres, $p=8$', 'two_p8', ['MS', 'GMHC']),
    ('three spheres, $\\cap$', 'intesf', ['MS', 'GMHC']),
    ('eight tori', 'toro8', ['MS', 'GMHC']),
    ('knot of spheres', 'knot', ['MS', 'GMHC']),
]
out = ['\\begin{tabular}{@{}llcrrrr@{}}', '\\toprule',
       'Manifold & $(n,k)$ & grid & method & time (s) & cells & $\\chi$\\\\', '\\midrule']
for name, key, algs in cases:
    first = True
    for a in algs:
        f = f'{key}_{a}.pol'
        if f not in rows:
            continue
        alg, n, k, D, t = rows[f]
        nn, kk, DD, cells = read_pol(R + f)
        counts, chi = euler(cells)
        g = 'x'.join(D.strip('[]').split())
        dims = D.strip('[]').split()
        grid = f'${dims[0]}^{len(dims)}$' if len(set(dims)) == 1 else '$' + '\\times'.join(dims) + '$'
        out.append(f"{name if first else ''} & {('$(' + n + ',' + k + ')$') if first else ''} & {grid if first else ''} & {({'MHC':'CMH','GMHC':'GCMH','GCHC':'GCCH'}).get(a,a)} & {t:.1f} & {len(cells)} & {chi}\\\\")
        first = False
    if not first:
        out.append('\\addlinespace')
out[-1:] = [] if out[-1] == '\\addlinespace' else out[-1:]
out += ['\\bottomrule', '\\end{tabular}']
os.makedirs(REPO + '/book/tables', exist_ok=True)
open(REPO + '/book/tables/octave_runs.tex', 'w').write('\n'.join(out) + '\n')
print('\n'.join(out))
