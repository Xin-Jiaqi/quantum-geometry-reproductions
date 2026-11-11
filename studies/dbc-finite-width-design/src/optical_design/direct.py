"""Full-BZ Gaussian integral with energy-aligned panels.
Adapted from round-3 adaptive_reference.py. Both optical pairs retained,
no histograms or small energy-denominator cutoffs; Gaussian 8-sigma tails omitted.
Coarse derivative scans locate monotonic intervals; density convergence is tested.
"""
from __future__ import annotations
import time
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import brentq
from .models import energy,deriv_x,kernels

def gauss_axis(breaks,order):
    a=np.unique(breaks);t,w=leggauss(order)
    x=((a[:-1,None]+a[1:,None])/2+(a[1:,None]-a[:-1,None])/2*t).ravel()
    wt=((a[1:,None]-a[:-1,None])/2*w).ravel()
    return x,wt

def regions(x,y):
    cy=np.minimum(abs(y-np.pi),abs(y+np.pi))
    c=(abs(x)<.15)&(cy<.5)
    n=(x*x+(abs(y)-2.10536)**2<.25**2)&(~c)
    return c,n,~(c|n)

def solve_brackets(xa,xb,y,band,target,p,derivative=False,iterations=38):
    fn=deriv_x if derivative else energy
    ia=np.arange(len(y));fa=fn(xa,y,p)[ia,band]-target
    a=xa.copy();b=xb.copy()
    for _ in range(iterations):
        m=(a+b)/2;fm=fn(m,y,p)[ia,band]-target
        left=fa*fm<=0;b=np.where(left,m,b);a=np.where(left,a,m);fa=np.where(left,fa,fm)
    return (a+b)/2

def integrate_point(p,photon,gamma,yorder=8,ypanels=32,xorder=12,axes=(0,0,1),shape="gaussian",xpanels=64):
    if gamma<=0:raise ValueError('positive width required')
    if shape not in ('gaussian','lorentzian'):raise ValueError('Unknown line shape')
    tic=time.perf_counter();ky=2.10536
    s=np.array([0,.00005,.0001,.0002,.0005,.001,.002,.005,.01,.02,.04,.08,.12,.2,.3,.4,.5,.7,1.])
    yb=np.r_[np.linspace(-np.pi,np.pi,ypanels+1),ky+s,ky-s,-ky+s,-ky-s,np.pi-s,-np.pi+s]
    # Additional photon-dependent branch onset points at x=0.
    coarse=np.linspace(-np.pi,np.pi,801);en=energy(np.zeros_like(coarse),coarse,p)
    extra=[]
    zs=np.array([-8,-4,-2,-1,0,1,2,4,8]) if shape=='gaussian' else np.r_[-np.array([1024,256,64,16,8,4,2,1]),0,np.array([1,2,4,8,16,64,256,1024])]
    levels=zs*gamma+photon
    for b in range(2):
        for target in levels:
            f=en[:,b]-target;ix=np.flatnonzero(f[:-1]*f[1:]<0)
            for i in ix:
                extra.append(brentq(lambda y:float(energy(0.,y,p)[b]-target),coarse[i],coarse[i+1]))
    yb=np.r_[yb,extra];yb=np.unique(np.clip(yb,-np.pi,np.pi));y,wy=gauss_axis(yb,yorder)
    xs=np.unique(np.r_[np.linspace(-np.pi,np.pi,xpanels+1),-.02,-.005,0,.005,.02])
    X,Y=np.meshgrid(xs,y,indexing='ij');v=deriv_x(X,Y,p)
    ix,iy,ib=np.nonzero(v[:-1]*v[1:]<0)
    xr=solve_brackets(xs[ix],xs[ix+1],y[iy],ib,np.zeros(len(ix)),p,derivative=True)
    per_y=[list(xs) for _ in y]
    for x,j in zip(xr,iy):per_y[j].append(float(x))
    left=[];right=[];ys=[];ids=[]
    for j,breaks in enumerate(per_y):
        xx=np.unique(breaks);left.extend(xx[:-1]);right.extend(xx[1:]);ys.extend([y[j]]*(len(xx)-1));ids.extend([j]*(len(xx)-1))
    left=np.asarray(left);right=np.asarray(right);ys=np.asarray(ys);ids=np.asarray(ids)
    el=energy(left,ys,p);er=energy(right,ys,p)
    rr=[];jj=[]
    for target in levels:
        il,ib=np.nonzero((el-target)*(er-target)<0)
        xx=solve_brackets(left[il],right[il],ys[il],ib,np.full(len(il),target),p)
        rr.extend(xx);jj.extend(ids[il])
    for x,j in zip(rr,jj):per_y[j].append(float(x))
    left=[];right=[];ys=[];weights=[]
    for j,breaks in enumerate(per_y):
        xx=np.unique(breaks);left.extend(xx[:-1]);right.extend(xx[1:]);ys.extend([y[j]]*(len(xx)-1));weights.extend([wy[j]]*(len(xx)-1))
    left=np.asarray(left);right=np.asarray(right);ys=np.asarray(ys);weights=np.asarray(weights)
    el=energy(left,ys,p);er=energy(right,ys,p)
    sel=np.any((np.maximum(el,er)>=photon-8*gamma)&(np.minimum(el,er)<=photon+8*gamma),axis=1) if shape=='gaussian' else np.ones(len(left),dtype=bool)
    left,right,ys,weights=[a[sel] for a in (left,right,ys,weights)]
    t,w=leggauss(xorder);X=(left[:,None]+right[:,None])/2+(right[:,None]-left[:,None])/2*t
    Y=np.broadcast_to(ys[:,None],X.shape);W=weights[:,None]*(right-left)[:,None]/2*w/(2*np.pi)**2
    result=np.zeros(8);pair_result=np.zeros(2);count=0
    for st in range(0,len(left),2000):
        x,z,wt=X[st:st+2000],Y[st:st+2000],W[st:st+2000]
        tr,K,D,O,E=kernels(x,z,p,axes);c,n,r=regions(x,z)
        line=(np.exp(-.5*((tr-photon)/gamma)**2)/(np.sqrt(2*np.pi)*gamma) if shape=='gaussian' else gamma/(np.pi*((tr-photon)**2+gamma**2)))
        g=line*wt[...,None]
        values=(K,np.maximum(K,0),np.minimum(K,0),K*c[...,None],K*n[...,None],K*r[...,None],np.ones_like(K),O)
        result+=np.array([np.sum(q*g) for q in values]);pair_result+=np.sum(K*g,axis=tuple(range(K.ndim-1)));count+=x.size
    return dict(lower_optical_pair=float(pair_result[0]),upper_optical_pair=float(pair_result[1]),shape=shape,xpanels=xpanels,parameters=vars(p),mass=p.mass,photon=photon,gamma=gamma,total=float(result[0]),positive=float(result[1]),negative=float(result[2]),critical_patch=float(result[3]),nodes=float(result[4]),remainder=float(result[5]),jdos=float(result[6]),optical=float(result[7]),y_nodes=len(y),kernel_nodes=count,xorder=xorder,yorder=yorder,ypanels=ypanels,seconds=time.perf_counter()-tic)
