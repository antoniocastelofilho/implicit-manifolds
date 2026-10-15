from j1a import *
def sim(G,b,s,pi,t=None):
    m=G.m;c,hs=G.box(b);Phi=refined_faces(G,b)
    steps=[(pi[j],s[pi[j]]) for j in range(m)]
    h=None
    for hh in range(1,m+1):
        if any(phi<=frozenset(steps[:hh]) for phi in Phi): h=hh;break
    v=[Fr(0)]*m;P=[tuple(v)]
    if h is None or h==m or t is None:
        for j in range(m): v[pi[j]]+=s[pi[j]];P.append(tuple(v))
    else:
        for j in range(h-1): v[pi[j]]+=s[pi[j]];P.append(tuple(v))
        v[pi[h-1]]+=s[pi[h-1]]
        for j in range(h,m): v[pi[j]]+=Fr(s[pi[j]],2)
        P.append(tuple(v))
        for j in range(h,m): v[pi[j]]+=Fr(t[j-h],2);P.append(tuple(v))
    return h,[tuple(ci+hs*vi for ci,vi in zip(c,p)) for p in P]
for m,N in [(2,4),(3,3)]:
    G=Grid(m,N); x=tuple([Fr(13,10)]*m)
    for b in G.leaves_at(x):
        if b in G.leaves: G.refine(b)
    ok=bad=0
    for F in [b for b in G.leaves if b[0]==1]:
        cF,hF=G.box(F)
        for s in itertools.product((-1,1),repeat=m):
            for pi in itertools.permutations(range(m)):
                h,P=sim(G,F,s,pi)
                if h is not None and h<m: continue
                # facet opposite v0 and the leaf across it
                a=pi[0]; y=list(cF); y[a]=cF[a]+s[a]*(hF+hF/4)
                for j in range(m):
                    if j!=a: y[j]=cF[j]+ (hF/8 if True else 0)
                Ls=[L for L in G.leaves_at(tuple(P[1][i]+ (s[a]*hF/4 if i==a else 0) for i in range(m)))]
                Cs=[L for L in Ls if L[0]==0]
                if not Cs: continue
                C=Cs[0]; cC,hC=G.box(C)
                s2=[0]*m; s2[a]=-s[a]
                for j in pi[1:]: s2[j]=1 if cF[j]>cC[j] else -1
                t=[s[j] for j in pi[1:]]
                h2,Q=sim(G,C,tuple(s2),pi,t)
                if h2==1 and set(Q[1:])==set(P[1:]): ok+=1
                else: bad+=1
    print(m,'fine->coarse pivots ok',ok,'bad',bad)
