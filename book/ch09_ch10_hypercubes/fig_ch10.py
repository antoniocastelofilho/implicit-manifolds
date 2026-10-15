import os as _os
REPO = _os.environ.get('IMPLICIT_MANIFOLDS_ROOT') or next(_p for _p in (_os.path.abspath(_os.path.join(_os.path.dirname(__file__), *(['..'] * _i))) for _i in range(8)) if _os.path.exists(_os.path.join(_p, 'CITATION.cff')))
import sys, pickle; sys.path.insert(0,REPO + '/python/implicit_manifolds')
import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from zbuf import render, crop
OUT=REPO + '/book/figures/'
D=pickle.load(open(REPO + '/book/ch09_ch10_hypercubes/fch_tests.pkl','rb'))
def cols(v,cmap):
    v=np.asarray(v); cm=plt.get_cmap(cmap); return [cm((x-v.min())/(np.ptp(v)+1e-12))[:3] for x in v]
ims=[]; titles=[]
P=D['S4']['polys']; Q=[p[:,:3] for p in P]
ims.append(crop(render(Q,cols([q[:,2].mean() for q in Q],'viridis'),elev=25,azim=-55,size=(900,900),edge_alpha=0.45)))
titles.append('$S^4\\subset\\mathbb{R}^5$: $x_4=0.3$, $x_5=0.2$')
P=D['S5']['polys']; Q=[p[:,:3] for p in P]
ims.append(crop(render(Q,cols([q[:,2].mean() for q in Q],'viridis'),elev=25,azim=-55,size=(900,900),edge_alpha=0.45)))
titles.append('$S^5\\subset\\mathbb{R}^7$: $x_4=0.2$, $x_5=0.3$, $x_6=0.1$')
open(REPO + '/book/ch09_ch10_hypercubes/ch10_counts.txt','w').write(str({k:(v['total'],v['kept'],len(v['polys'])) for k,v in D.items()}))
# T4 in R8 (cut x6=x8=0.6): four tori, labelled by signs of x5 and x7; stereographic view of (x1..x4)
P=D['T4']['polys']
def stereo(Q,r): return r*Q[:,:3]/(r-Q[:,3:4]+1e-9)
X=np.vstack(P); r=np.sqrt((X[:,:4]**2).sum(1)).max()
groups={}
for p in P: groups.setdefault((np.sign(p[:,4].mean()),np.sign(p[:,6].mean())),[]).append(p)
print({k:len(v) for k,v in groups.items()}, r)
tor=[]
for key in sorted(groups):
    G=groups[key]; Q=[stereo(p,r*1.0001) for p in G]
    c=[plt.get_cmap('hsv')((np.arctan2(p[:,1].mean(),p[:,0].mean())+np.pi)/(2*np.pi))[:3] for p in G]
    tor.append(crop(render(Q,c,elev=35,azim=-60,size=(700,600),edge_alpha=0.5)))
fig=plt.figure(figsize=(12,7.6))
gs=fig.add_gridspec(2,4,height_ratios=[1.25,1])
for i in range(2):
    a=fig.add_subplot(gs[0,2*i:2*i+2]); a.imshow(ims[i]); a.axis('off'); a.set_title(titles[i],fontsize=11)
for j,(key,im) in enumerate(zip(sorted(groups),tor)):
    a=fig.add_subplot(gs[1,j]); a.imshow(im); a.axis('off')
    a.set_title('$T^4$: $x_5%s0$, $x_7%s0$'%('<' if key[0]<0 else '>','<' if key[1]<0 else '>'),fontsize=10)
plt.subplots_adjust(0.01,0.01,0.99,0.95,0.04,0.1); plt.savefig(OUT+'fch_tests.png',dpi=200)
from PIL import Image
im=Image.open(OUT+'fch_tests.png'); im.thumbnail((1400,900)); im.save('/tmp/claude-0/-home-claude/897a473a-5eac-5851-a062-2767f231ce50/scratchpad/v_fch.png')
