import os as _os
REPO = _os.environ.get('IMPLICIT_MANIFOLDS_ROOT') or next(_p for _p in (_os.path.abspath(_os.path.join(_os.path.dirname(__file__), *(['..'] * _i))) for _i in range(8)) if _os.path.exists(_os.path.join(_p, 'CITATION.cff')))
import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from math import factorial, comb, log10
plt.rcParams.update({'font.size':9,'mathtext.fontset':'cm','axes.spines.top':False,'axes.spines.right':False})
n=np.arange(2,13)
series=[('Marching Cubes table, $2^{2^n}$ entries',[2**k*log10(2) for k in n],'#2a78d6'),
        ('Kuhn simplices per cube, $n!$',[log10(factorial(k)) for k in n],'#eb6834'),
        ('faces of a cube, $3^n$',[k*log10(3) for k in n],'#1baf7a'),
        ('edges of a cube, $n\\,2^{n-1}$',[log10(k*2**(k-1)) for k in n],'#eda100')]
fig,ax=plt.subplots(figsize=(6.2,3.4))
for lab,y,c in series:
    ax.plot(n,y,'-o',color=c,lw=2,ms=4,label=lab)
ax.set_yscale('log'); ax.set_xlabel('dimension $n$'); ax.set_ylabel('number of decimal digits')
ax.set_yticks([0.5,1,2,5,10,20,50,100,300,1000]); ax.set_yticklabels(['0.5','1','2','5','10','20','50','100','300','1000'])
ax.grid(axis='y',color='0.88',lw=0.6); ax.set_xticks(n)
ax.legend(frameon=False,fontsize=8,loc='upper left')
plt.tight_layout(); plt.savefig(REPO + '/book/figures/intro_growth.pdf')
