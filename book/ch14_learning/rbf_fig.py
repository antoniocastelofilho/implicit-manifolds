"""Figure and table of the RBF / Hermite RBF comparison (parses the logs of rbf_compare.py)."""
import os as _os
REPO = _os.environ.get('IMPLICIT_MANIFOLDS_ROOT') or next(_p for _p in (_os.path.abspath(_os.path.join(_os.path.dirname(__file__), *(['..'] * _i))) for _i in range(8)) if _os.path.exists(_os.path.join(_p, 'CITATION.cff')))
import ast, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 8.5, 'mathtext.fontset': 'cm'})
R = REPO + '/book/runs/ml/'
def parse(f, kinds=None):
    out = {}
    for l in open(R + f):
        if "{'kind'" not in l: continue
        t, d = l.replace('np.float64(', '').replace(')', '').split(' ', 1); r = ast.literal_eval(d)
        if kinds and r['kind'] not in kinds: continue
        out[(int(t.split('_')[1]), r['kind'])] = r
    return out
S = {'r3': {**parse('rbf_r3.log', ['poly']), **parse('rbf_r3b.log')}, 'r4': parse('rbf_r4.log'), 'r4n': parse('rbf_r4n.log')}
grids = {'r3': [(16, '$16^2{\\times}6$'), (24, '$24^2{\\times}8$'), (36, '$36^2{\\times}12$')],
         'r4': [(8, '$8^2{\\times}6^2$'), (12, '$12^2{\\times}8^2$'), (16, '$16^2{\\times}10^2$')]}
grids['r4n'] = grids['r4']
titles = {'r3': 'rippled torus, $\\mathbb{R}^3$, $\\sigma=0.002$', 'r4': 'rippled torus, $\\mathbb{R}^4$, $\\sigma=0.002$',
          'r4n': 'rippled torus, $\\mathbb{R}^4$, $\\sigma=0.01$'}
lab = {'poly': 'quadric', 'rbf': 'RBF', 'hrbf': 'Hermite RBF'}; mk = {'poly': 'o-', 'rbf': 's-', 'hrbf': '^-'}
fig, axs = plt.subplots(1, 3, figsize=(9, 2.9))
for a, key in zip(axs, ['r3', 'r4', 'r4n']):
    g = grids[key]
    for kd in ['poly', 'rbf', 'hrbf']:
        a.semilogy(range(3), [S[key][(gg, kd)]['rms_z'] for gg, _ in g], mk[kd], label=lab[kd], ms=4)
        a.semilogy(range(3), [S[key][(gg, kd)]['mx_z'] for gg, _ in g], mk[kd][0] + ':', ms=3, color=a.lines[-1].get_color(), alpha=.6)
    a.axhline(0.002 if key != 'r4n' else 0.01, color='gray', lw=.6, ls='--')
    a.set_xticks(range(3)); a.set_xticklabels([l for _, l in g]); a.set_title(titles[key], fontsize=8.5)
axs[0].set_ylabel('distance to the manifold'); axs[0].legend(fontsize=7.5, loc='upper right')
plt.tight_layout(); plt.savefig(REPO + '/book/figures/ml_rbf.pdf'); plt.close()
# table
rows = []
for key, name in [('r3', '$(3,1)$, $\\sigma{=}0.002$'), ('r4', '$(4,2)$, $\\sigma{=}0.002$'), ('r4n', '$(4,2)$, $\\sigma{=}0.01$')]:
    for i, (gg, gl) in enumerate(grids[key]):
        c = [S[key][(gg, kd)] for kd in ['poly', 'rbf', 'hrbf']]
        best = min(r['rms_z'] for r in c)
        def f(r): s = '%.2f' % (1e3 * r['rms_z']); return '\\textbf{%s}' % s if r['rms_z'] == best else s
        rows.append(('\\multirow{3}{*}{%s}' % name if i == 0 else '') + ' & ' + gl + ' & ' + ' & '.join(f(r) for r in c) + ' & '
                    + ' & '.join('%.1f' % r['t_build'] for r in c) + ' & ' + ' & '.join('%.0f' % r['t_eval'] for r in c)
                    + ' & %d\\%%' % round(100 * c[1]['rbf_cells'] / c[1]['data']) + ' \\\\')
    rows.append('\\midrule')
rows = rows[:-1]
T = r"""\begin{table}[t]
\caption{Quadric, RBF and Hermite RBF local approximations in the partition of unity. Error: root mean square
distance (times $10^{3}$) to the exact manifold of the vertices of the approximation after Gauss--Newton projection
onto $F^{-1}(0)$, which measures the error of $F$ itself; the best value of each row is in bold. Times: construction
of $F$ in seconds and one evaluation of $F$ in microseconds. Last column: percentage of the data cells that passed
the coverage test and received a radial correction. All approximations have the correct Euler characteristic $0$.}
\label{tab:ml:rbf}
\centering\footnotesize\setlength{\tabcolsep}{2.0pt}
\begin{tabular}{@{}llccccccccc r@{}}
\toprule
& & \multicolumn{3}{c}{error ($\times10^{3}$)} & \multicolumn{3}{c}{construction (s)} & \multicolumn{3}{c}{evaluation ($\mu$s)} & \\
\cmidrule(lr){3-5}\cmidrule(lr){6-8}\cmidrule(lr){9-11}
$(n,k)$, noise & grid & quad. & RBF & HRBF & quad. & RBF & HRBF & quad. & RBF & HRBF & RBF\\
\midrule
""" + '\n'.join(rows) + r"""
\bottomrule
\end{tabular}
\end{table}
"""
open(REPO + '/book/tables/ml_rbf.tex', 'w').write(T)
for k in S:
    for kk, r in sorted(S[k].items()): assert r['chi'] == 0, (k, kk)
print(T)
