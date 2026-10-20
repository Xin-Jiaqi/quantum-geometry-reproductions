"""2D BZ integration with non-overlapping circular patches around K and K'.

All weights represent d^2k/(2pi)^2. BZ exterior uses midpoint quadrature,
nearby valley disks use equal-area polar coordinates. Valley exclusion is
performed on a torus (including G-translated copies); physical gauge
covariance under k+G is tested separately.
"""
import numpy as np

def valley_quadrature(reciprocal, exterior_n=121, radius=1.5, radial_n=50, angular_n=96, radial_power=2.):
    G=np.asarray(reciprocal,float)
    if G.shape!=(2,2) or not radius>0 or min(exterior_n,radial_n,angular_n)<1 or radial_power < 1:
        raise ValueError('bad quadrature setup')
    frac_centers=np.array([[1/3,2/3],[2/3,1/3]],float)
    centers=frac_centers@G
    replicated=centers[:,None,:]+np.array([[a,b] for a in (-1,0,1) for b in (-1,0,1)])@G
    replicated=replicated.reshape(-1,2)
    # Disks smaller than half the minimum separation of periodic valley images.
    dif=replicated[:,None,:]-replicated[None,:,:]
    dis=np.linalg.norm(dif,axis=-1)
    dis[dis<1e-9]=np.inf
    if 2*radius>=float(dis.min()):
        raise ValueError('valley disks overlap')
    uv=(np.arange(exterior_n)+.5)/exterior_n
    f=np.stack(np.meshgrid(uv,uv,indexing='ij'),-1).reshape(-1,2)
    kout=f@G
    nearest2=np.min(np.sum((kout[:,None,:]-replicated[None,:,:])**2,axis=-1),axis=1)
    kout=kout[nearest2>=radius**2]
    wout=np.full(len(kout),abs(float(np.linalg.det(G)))/(2*np.pi)**2/exterior_n**2)
    shell_edges=radius*(np.arange(radial_n+1)/radial_n)**radial_power
    r=(shell_edges[1:]+shell_edges[:-1])/2
    ring_weights=np.pi*(shell_edges[1:]**2-shell_edges[:-1]**2)/((2*np.pi)**2*angular_n)
    angles=2*np.pi*(np.arange(angular_n)+.5)/angular_n
    radial=np.stack((r[:,None]*np.cos(angles),r[:,None]*np.sin(angles)),axis=-1).reshape(-1,2)
    kin=np.concatenate([center[None,:]+radial for center in centers],axis=0)
    win=np.tile(np.repeat(ring_weights,angular_n),len(centers))
    return np.concatenate((kout,kin)),np.concatenate((wout,win)),{
        'n_exterior_kept':len(kout),'n_valley_total':len(kin),
        'exterior_area_weight':float(wout.sum()), 'valley_area_weight':float(win.sum()),
        'bz_area_weight':float(abs(np.linalg.det(G))/(2*np.pi)**2),
        'normalization_error':float(wout.sum()+win.sum()-abs(np.linalg.det(G))/(2*np.pi)**2)}
