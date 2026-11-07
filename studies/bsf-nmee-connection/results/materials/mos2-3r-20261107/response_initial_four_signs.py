"""Fixed-ion, self-consistent, layer-resolved EB spin response. PBS only."""
from pathlib import Path
import hashlib,json,os,re,shutil,socket,subprocess
import numpy as np
if not os.getenv('PBS_JOBID') or not os.getenv('PBS_NODEFILE') or socket.gethostname().startswith(('super','cms-ssh')):
    raise SystemExit('PBS compute node required')
root=Path.cwd(); label=os.environ['STRUCTURE']; gate=float(os.environ.get('GATE_VNM','0'))
mesh=int(os.environ.get('MESH','6')); electric=float(os.environ.get('PROBE','0.004')); zeeman=float(os.environ.get('ZEEMAN','0.016'))
ranks=len(Path(os.environ['PBS_NODEFILE']).read_text().splitlines())
if electric<=0 or zeeman<=0:raise ValueError('positive probe steps required')
base=root/f'{label}_F{gate:+.3f}_k{mesh}_e{electric}_b{zeeman}';base.mkdir(exist_ok=True)
lines=(root/f'POSCAR.{label}').read_text().splitlines();nat=sum(map(int,lines[6].split()));ne=24 if nat==3 else 48
cell=np.array([list(map(float,s.split())) for s in lines[2:5]])*float(lines[1]);area=np.linalg.norm(np.cross(cell[0],cell[1]))
common=f'''SYSTEM=MoS2 non-Janus NMEE {label}
ENCUT=450
PREC=Accurate
LASPH=T
LREAL=F
ISYM=-1
LSORBIT=T
SAXIS=0 0 1
MAGMOM={3*nat}*0
NBANDS={32 if nat==3 else 56}
ISMEAR=0
SIGMA=0.02
EDIFF=1e-9
NELM=240
NSW=0
IBRION=-1
LWAVE=T
LCHARG=T
KPAR={min(ranks,4)}
NCORE=1
'''
# GATE_VNM uses the conventional physical electric-field direction.
# VASP EFIELD is an electric force field with the opposite sign.
dipole=f'LDIPOL=T\nIDIPOL=3\nDIPOL=0.5 0.5 0.5\nEFIELD={-gate/10:.10f}\n'
runner='/work/app/VASP6/vasp.6.4.3_wannier90_cpu/bin/vasp_ncl'

def run(name,extra,restart=None):
    d=base/name;d.mkdir(exist_ok=True)
    expected=common+extra
    if (d/'finished.json').exists():
        if (d/'INCAR').read_text()!=expected:raise ValueError('input changed in completed directory')
        return d
    shutil.copyfile(root/f'POSCAR.{label}',d/'POSCAR')
    with (d/'POTCAR').open('wb') as f:
        for element in ('Mo_pv','S'):f.write((Path('/work/sunq/mxy/bin/potpaw_PBE.64')/element/'POTCAR').read_bytes())
    (d/'KPOINTS').write_text(f'Full Gamma mesh\n0\nGamma\n{mesh} {mesh} 1\n0 0 0\n')
    (d/'INCAR').write_text(expected)
    if restart:
        for f in ('WAVECAR','CHGCAR'):shutil.copyfile(restart/f,d/f)
    with (d/'vasp.log').open('w') as f:
        subprocess.run(['mpirun','-np',str(ranks),'-ppn',str(ranks),'-hostfile',os.environ['PBS_NODEFILE'],runner],cwd=d,stdout=f,stderr=subprocess.STDOUT,check=True)
    out=(d/'OUTCAR').read_text(errors='replace');log=(d/'vasp.log').read_text(errors='replace')
    if 'General timing' not in out or 'aborting loop because EDIFF is reached' not in out or re.search(r'^\s*(DAV|RMM|SDA|DMP|CGA|CG)\s*:\s*240\s',log,re.M):
        raise RuntimeError(f'SCF not converged: {name}')
    if name.startswith('E') and 'p_tot' not in out+log:raise RuntimeError('missing electric enthalpy')
    e=(d/'EIGENVAL').read_text().splitlines();nel,nk,nb=map(int,e[5].split());j=6;v=[];c=[]
    if nel!=ne:raise ValueError('unexpected electron number')
    for k in range(nk):
        while not e[j].strip():j+=1
        j+=1;v.append(float(e[j+ne-1].split()[1]));c.append(float(e[j+ne].split()[1]));j+=nb
    gap=min(c)-max(v)
    if gap<=0:raise RuntimeError('insulating response protocol invalid: gap closed')
    record={'gap_eV':gap,'node':socket.gethostname(),'pbs':os.environ['PBS_JOBID']}
    (d/'finished.json').write_text(json.dumps(record))
    print(name,record,flush=True)
    return d

def spin(d):
    with (d/'CHGCAR').open() as f:
        header=[next(f) for _ in range(8)]
        for _ in range(sum(map(int,header[6].split()))):next(f)
        dims=next(s for s in f if s.strip()).split();nx,ny,nz=map(int,dims)
        if nz%2:raise ValueError('even z grid required for symmetric layer integration')
        weights=np.zeros(nz);weights[1:nz//2]=1;weights[[0,nz//2]]=.5
        answer=[]
        for axis in range(4):
            if axis:
                for line in f:
                    if line.split()==dims:break
                else:raise RuntimeError('missing spin density')
            vals=[]
            while len(vals)<nx*ny*nz:vals.extend(map(float,next(f).split()))
            plane=np.asarray(vals).reshape(nz,ny,nx).mean(axis=(1,2))/nz
            total=float(plane.sum());bottom=float(weights@plane)
            answer.append([total,bottom,total-bottom])
        if abs(answer[0][0]-ne)>.01:raise ValueError('CHGCAR charge normalization')
        return answer

start=run('start','ALGO=Normal\nISTART=0\nICHARG=2\nBEXT=0 0 0\n')
rows={};gaps={}
for bs in (1,-1):
    seed=run(f'seedB{bs:+d}',dipole+f'ALGO=Normal\nISTART=1\nICHARG=1\nBEXT=0 0 {bs*zeeman}\n',start)
    for es in (1,-1):
        name=f'E{es:+d}B{bs:+d}'
        d=run(name,dipole+f'ALGO=All\nISTART=1\nICHARG=1\nBEXT=0 0 {bs*zeeman}\nEFIELD_PEAD={es*electric} 0 0\nLPEAD_RELAX=T\nEDIFF_PEAD=1e-9\nSKIP_EDOTP=F\nISEARCH=1\nTIME=.30\n',seed)
        rows[name]=spin(d);gaps[name]=json.loads((d/'finished.json').read_text())['gap_eV']
mu=2.00231930436256*5.7883818060e-5/2
factor=1e10*mu/area
cp=(np.array(rows['E+1B+1'])[1:]-np.array(rows['E-1B+1'])[1:])*factor/(2*electric*zeeman)
cm=-(np.array(rows['E+1B-1'])[1:]-np.array(rows['E-1B-1'])[1:])*factor/(2*electric*zeeman)
c=(cp+cm)/2
if np.max(abs(c[:,0]-c[:,1]-c[:,2]))>1e-7:raise ValueError('layer sum')
result={'structure':label,'gate_physical_Vnm':gate,'vasp_EFIELD':-gate/10,'probe_VA':electric,'zeeman_eV':zeeman,'mesh':mesh,'area_A2':float(area),'C_xyz_total_bottom_top_1e-9':c.tolist(),'B_even_residual_xyz':((cp-cm)/2).tolist(),'one_B_sign_C_xyz':[cp.tolist(),cm.tolist()],'spin_raw':rows,'gaps_eV':gaps,'pbs':os.environ['PBS_JOBID'],'node':socket.gethostname(),'structure_sha256':hashlib.sha256((root/f'POSCAR.{label}').read_bytes()).hexdigest(),'mono_halfspace_warning':label=='mono'}
(base/'answer.json').write_text(json.dumps(result,indent=2));print('ANSWER',json.dumps(result),flush=True)
