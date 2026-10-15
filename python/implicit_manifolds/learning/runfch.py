import os as _os
REPO = _os.environ.get('IMPLICIT_MANIFOLDS_ROOT') or next(_p for _p in (_os.path.abspath(_os.path.join(_os.path.dirname(__file__), *(['..'] * _i))) for _i in range(8)) if _os.path.exists(_os.path.join(_p, 'CITATION.cff')))
import sys, os, time; sys.path.insert(0,REPO + '/python/fch_lib')
import numpy as np
from fch_wrapper import FCH_CPU_Runner
OUT=REPO + '/book/runs/ml/'
def run_fch(shim, n, k, dmin, dmax, div, pre_root, tag):
    os.makedirs(OUT+tag, exist_ok=True)
    r = FCH_CPU_Runner(REPO + '/python/fch_lib/lib/')
    t0=time.time()
    r.find_first_simplex(ndim=n, kdim=k, math_func=shim, domain_min=np.asarray(dmin,float), domain_max=np.asarray(dmax,float),
                         domain_div=np.asarray(div,np.uint32), out_path=OUT+tag+'/params.bin', PRE_ROOT=np.asarray(pre_root,float), max_iter_per_exponent=3000)
    p = r.read_find_first_simplex_binary(OUT+tag+'/params.bin')
    t1=time.time()
    r.run_threadpool(math_func=shim, domain_min=p['dmin'], domain_max=p['dmax'], domain_div=p['ddiv'], first_coord=p['first_coord'],
                     first_eidx=p['first_eidx'], out_path=OUT+tag+'/results.bin', print_interval=10**7)
    t2=time.time()
    r.skeleton(in_path=OUT+tag+'/results.bin', out_path=OUT+tag+'/results.pol')
    t3=time.time()
    return OUT+tag+'/results.pol', dict(first=t1-t0, fch=t2-t1, skeleton=t3-t2)
