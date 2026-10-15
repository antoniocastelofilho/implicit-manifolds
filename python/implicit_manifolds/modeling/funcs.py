import numpy as np
def knotC():
    t=4*np.pi*np.arange(151)/150
    return np.c_[(1+0.4*np.cos(1.5*t))*np.cos(t),(1+0.4*np.cos(1.5*t))*np.sin(t),0.5*np.sin(1.5*t)]
def tubeC(N=200):
    s=2*np.pi*np.arange(N)/N
    return np.c_[(1+0.35*np.cos(3*s))*np.cos(2*s),(1+0.35*np.cos(3*s))*np.sin(2*s),0.45*np.sin(3*s),0.5*np.cos(5*s)]
def uminp(fs,p): fs=np.asarray(fs); return (fs**(-p)).sum()**(-1/p)
def imaxp(fs,p): fs=np.asarray(fs); return (fs**p).sum()**(1/p)
TORC=[(6,0),(-6,0),(0,6),(0,-6),(4,4),(4,-4),(-4,4),(-4,-4)]
F={
 'two_spheres': lambda x: min(x@x,(x-[1,0,0])@(x-[1,0,0])),
 'two_spheres_p8': lambda x: uminp([x@x,(x-[1,0,0])@(x-[1,0,0])],8),
 'eight_tori': lambda x: uminp([(np.hypot(x[0]-a,x[1]-b)-2)**2+x[2]**2 for a,b in TORC],10),
 'knot': lambda x,C=knotC(): uminp(((x-C)**2).sum(1)-0.24**2+1,30),
 'intersect3': lambda x: imaxp([x@x,(x-[0,1,0])@(x-[0,1,0]),(x-[np.sqrt(3)/2,0.5,0])@(x-[np.sqrt(3)/2,0.5,0])],10),
 'r4_two_balls': lambda x: uminp([(x-[0.8,0,0,0])@(x-[0.8,0,0,0]),(x+[0.8,0,0,0])@(x+[0.8,0,0,0])],8),
 'r4_diff': lambda x: max(x@x, 2-(x[0]**2+x[1]**2)/0.2025),
 'r4_tube': lambda x,C=tubeC(): uminp(((x-C)**2).sum(1)-0.22**2+1,30),
}
