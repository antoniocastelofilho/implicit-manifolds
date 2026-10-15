import sys, time, importlib, shutil, os
sys.path.insert(0, '.')
name = sys.argv[1]
tag = sys.argv[2]
shutil.copy(f'examples/{name}.py', 'my_math.py') if name != 'my_math' else None
import my_math as M
from fch_wrapper import FCH_CPU_Runner
r = FCH_CPU_Runner('./lib/')
os.makedirs(f'out/{tag}', exist_ok=True)
t=time.time()
r.find_first_simplex(ndim=M.ndim,kdim=M.kdim,math_func=M.my_math,domain_min=M.dmin,domain_max=M.dmax,domain_div=M.ddiv,out_path=f'out/{tag}/params.bin',PRE_ROOT=M.pre_root,max_iter_per_exponent=3000)
t1=time.time()
p = r.read_find_first_simplex_binary(f'out/{tag}/params.bin')
r.run_threadpool(math_func=M.my_math,domain_min=p['dmin'],domain_max=p['dmax'],domain_div=p['ddiv'],first_coord=p['first_coord'],first_eidx=p['first_eidx'],out_path=f'out/{tag}/results.bin',print_interval=100000)
t2=time.time()
r.skeleton(in_path=f'out/{tag}/results.bin', out_path=f'out/{tag}/results.pol')
t3=time.time()
e = r.read_fch_output_binary(f'out/{tag}/results.bin')
print(f'RESULT {tag} n={M.ndim} k={M.kdim} div={list(M.ddiv)[:6]} first={t1-t:.2f}s fch={t2-t1:.2f}s skeleton={t3-t2:.2f}s edges={e.shape}')
