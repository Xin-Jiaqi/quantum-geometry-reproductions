"""Stage 8b: IDEAL 1H-MX2 parallel bilayer full-atomic-motif symmetry audit.

Do not treat generated atomic positions as relaxed double-layer CIFs.
Only NumPy; no DFT. Tensor convention Gamma[axial spin, axial B, polar E].
"""
import math
import numpy as np

PARAMS={"WSe2":(3.28,1.67,3.165,"W","Se"),
        "MoS2":(3.1881,17.4639*.0894,3.0,"Mo","S")}
# WSe2 h and MoS2 interface gap are ILLUSTRATIVE, not relaxed bilayer data.


def structure(formula,shift):
    a,h,g,m,x=PARAMS[formula]
    lat=np.array(((a,0),(-a/2,a*math.sqrt(3)/2)))
    frac=[];species=[]
    for z0,offset in ((-(2*h+g)/2,(0,0)),((2*h+g)/2,shift)):
        for element,xy,z in ((m,(0,0),0),(x,(1/3,2/3),h),
                             (x,(1/3,2/3),-h)):
            frac.append([*( (np.array(offset)+xy)%1), z0+z])
            species.append(element)
    return np.array(frac),species,lat


def maps(src,elements,lat,dst,dst_elements,R,sz,tol=1e-7):
    """Test whole-atomic-motif map modulo 2D lattice and global xy origin."""
    xy=(src[:,:2]@lat@R.T)@np.linalg.inv(lat)
    z=sz*src[:,2]
    for j in range(len(dst)):
        if elements[0]!=dst_elements[j] or abs(z[0]-dst[j,2])>tol:
            continue
        translation=dst[j,:2]-xy[0]
        used=set()
        for i in range(len(src)):
            match=next((k for k in range(len(dst)) if k not in used
                and elements[i]==dst_elements[k]
                and abs(dst[k,2]-z[i])<tol
                and np.max(np.abs((xy[i]+translation-dst[k,:2]+.5)%1-.5))<tol),
                None)
            if match is None:
                break
            used.add(match)
        if len(used)==len(src):
            return True
    return False


def operations(frac,species,lat):
    refl=np.diag([-1.,1.])
    surviving=[]
    for n in range(6):
        phi=n*math.pi/3
        rot=np.array(((math.cos(phi),-math.sin(phi)),
                      (math.sin(phi),math.cos(phi))))
        for mirror in (False,True):
            R=rot@(refl if mirror else np.eye(2))
            for zflip in (1,-1):
                if maps(frac,species,lat,frac,species,R,zflip):
                    surviving.append((n,mirror,zflip,R))
    return surviving


def allows_xzx(ops):
    """Axial spin and axial Zeeman determinants cancel in rank-3 tensor."""
    basis=np.zeros((3,3,3));basis[0,2,0]=1.
    invariant=np.zeros_like(basis)
    for _,_,zflip,R2 in ops:
        R=np.eye(3);R[:2,:2]=R2;R[2,2]=zflip
        invariant+=np.einsum("ia,jb,kc,abc->ijk",
                             R,R,R,basis,optimize=True)
    return np.linalg.norm(invariant/len(ops))>1e-9


def audit():
    result={}
    for formula in PARAMS:
        stacks={}
        for name,tau in (("AA",(0,0)),("AB",(1/3,2/3)),
                         ("BA",(-1/3,-2/3)),("generic",(.13,.20))):
            f,sp,lat=structure(formula,tau)
            ops=operations(f,sp,lat)
            stacks[name]={"shift":list(tau),"order":len(ops),
                "C3":any(n==2 and not m and z==1 for n,m,z,R in ops),
                "mz":any(n==0 and not m and z==-1 for n,m,z,R in ops),
                "inversion":any(n==3 and not m and z==-1 for n,m,z,R in ops),
                "allowed_Gamma_xzx":bool(allows_xzx(ops))}
        ab,sa,L=structure(formula,(1/3,2/3))
        ba,sb,_=structure(formula,(-1/3,-2/3))
        result[formula]={"stacks":stacks,"AB_to_BA_by_mz":maps(
            ab,sa,L,ba,sb,np.eye(2),-1)}
    return result


if __name__=="__main__":
    import json
    print(json.dumps(audit(),indent=2))
