import os as _os
REPO = _os.environ.get('IMPLICIT_MANIFOLDS_ROOT') or next(_p for _p in (_os.path.abspath(_os.path.join(_os.path.dirname(__file__), *(['..'] * _i))) for _i in range(8)) if _os.path.exists(_os.path.join(_p, 'CITATION.cff')))
import numpy as np, sys
sys.path.insert(0,REPO + '/book/ch12_comparison'); from bratu import F as BF, J as BJ
def newton_lam(N, lam):
    u=np.zeros(N)
    for _ in range(50):
        x=np.r_[lam,u]; u=u-np.linalg.solve(BJ(x,N)[:,1:],BF(x,N))
    return np.r_[lam,u]
out=[]
def arr(a): return '{'+', '.join('%.17g'%v for v in a)+'}'
for N,div,umax in ((3,(24,20,20,30),(4.6,4.6,6.5)),(5,(12,10,10,12,10,10),(3.4,5.0,6.5,5.0,3.4))):
    n=N+1; x0=newton_lam(N,1.0)
    lo=[-0.131]+[-0.117-0.003*i for i in range(N)]; hi=[3.613]+[umax[i]+0.007*i for i in range(N)]
    body='\n'.join(['        const double h=1.0/%d.0;'%(N+1),'        double u[%d]; u[0]=0; u[%d]=0;'%(N+2,N+1),
        '        for(int i=1;i<=%d;++i) u[i]=x[i];'%N,
        '        for(int i=1;i<=%d;++i) f[i-1]=(u[i-1]-2*u[i]+u[i+1])/(h*h)+x[0]*std::exp(u[i]);'%N])
    out.append(f"""template <>
struct TestCaseImpl<TestCase::BOOK_BRATU{N}> {{                // discrete Bratu problem, x=(lambda,u_1..u_{N})
    static constexpr size_t NDIM = {n};
    static constexpr size_t KDIM = {N};
    static constexpr std::array<double, NDIM> DOMAIN_MIN{arr(lo)};
    static constexpr std::array<double, NDIM> DOMAIN_MAX{arr(hi)};
    static constexpr std::array<size_t, NDIM> DOMAIN_DIV{{{', '.join(map(str,div))}}};
    static constexpr std::array<double, NDIM> FIRST_POINT{arr(x0)};
    static void func(const double *x, double *f) {{
{body}
    }}
}};
""")
rng=np.random.default_rng(12)
for n,nb,div in ((3,14,(40,40,40)),(4,12,(24,24,24,24))):
    C=rng.uniform(-1.6,1.6,(nb,n)); r=rng.uniform(0.25,0.45,nb); T=0.5
    # first point: on the sphere-like level set around the first centre, along +x1
    f=lambda x: T-np.sum(np.exp(-np.sum((x-C)**2,1)/r**2))
    a=C[0].copy(); b=C[0]+np.r_[1.0,np.zeros(n-1)]
    for _ in range(80):
        mid=(a+b)/2
        if f(mid)<0: a=mid
        else: b=mid
    x0=(a+b)/2+np.r_[0,0.0013*np.arange(1,n)]*0
    terms='\n'.join(f'        {{ double d=0; for(int j=0;j<{n};++j){{ double t=x[j]-C[{i}][j]; d+=t*t; }} s+=std::exp(-d/({r[i]**2:.17g})); }}' for i in range(nb))
    Cs='{'+', '.join('{'+', '.join('%.17g'%v for v in c)+'}' for c in C)+'}'
    lo=[-2.03-0.007*i for i in range(n)]; hi=[2.01+0.005*i for i in range(n)]
    out.append(f"""template <>
struct TestCaseImpl<TestCase::BOOK_BLOBS{n}> {{                 // {nb} Gaussian blobs, level 1/2
    static constexpr size_t NDIM = {n};
    static constexpr size_t KDIM = 1;
    static constexpr std::array<double, NDIM> DOMAIN_MIN{arr(lo)};
    static constexpr std::array<double, NDIM> DOMAIN_MAX{arr(hi)};
    static constexpr std::array<size_t, NDIM> DOMAIN_DIV{{{', '.join(map(str,div))}}};
    static constexpr std::array<double, NDIM> FIRST_POINT{arr(x0)};
    static void func(const double *x, double *f) {{
        static constexpr double C[{nb}][{n}] = {Cs};
        double s=0;
{terms}
        f[0]=0.5-s;
    }}
}};
""")
    np.save(REPO + f'/book/ch12_comparison/blobs{n}.npy',np.c_[C,r])
open(REPO + '/cpp/common/test_cases/book_cases_ch12.h','w').write('#pragma once\n#include "test_case_template.h"\n#include <cmath>\n\n'+'\n'.join(out))
print('ok')
