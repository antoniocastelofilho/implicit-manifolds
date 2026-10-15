import os as _os
REPO = _os.environ.get('IMPLICIT_MANIFOLDS_ROOT') or next(_p for _p in (_os.path.abspath(_os.path.join(_os.path.dirname(__file__), *(['..'] * _i))) for _i in range(8)) if _os.path.exists(_os.path.join(_p, 'CITATION.cff')))
import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from PIL import Image
F=REPO + '/book/figures/'
def load(f, box=None):
    im=Image.open(F+f).convert('RGB'); W,H=im.size
    if box: im=im.crop((int(box[0]*W),int(box[1]*H),int(box[2]*W),int(box[3]*H)))
    a=np.asarray(im).astype(float)/255; m=(np.abs(a-1).sum(2)>0.03); ys,xs=np.where(m)
    return a[ys.min():ys.max()+1, xs.min():xs.max()+1]
panels=[(load('hc_torus_r10.png'),'(a) a torus in $\\mathbb{R}^{10}$ (Chap. 9)'),
        (load('klein.png'),'(b) a Klein bottle, graph in $\\mathbb{R}^5$ (Chap. 9)'),
        (load('knot.png',(0.52,0.16,1,1)),'(c) a knot of 151 spheres (Chap. 13)'),
        (load('cahn_hilliard100.png',(0,0.1,0.5,1)),'(d) equilibria in $\\mathbb{R}^{100}$ (Chap. 10)'),
        (load('pipebox.png',(0.67,0.1,1,1)),'(e) a blended CSG model (Chap. 13)'),
        (load('mod_r4_tube.png',(0.26,0.2,0.5,1)),'(f) a section of a tube in $\\mathbb{R}^4$ (Chap. 13)')]
fig,axs=plt.subplots(2,3,figsize=(11,6.6))
for a,(im,t) in zip(axs.ravel(),panels):
    a.imshow(im); a.axis('off'); a.set_title(t,fontsize=9.5)
plt.subplots_adjust(0.01,0.01,0.99,0.95,0.05,0.15); plt.savefig(F+'intro_gallery.png',dpi=200)
