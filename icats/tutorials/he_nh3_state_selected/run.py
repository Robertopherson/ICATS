"""Public-CLI preparation followed by inexpensive rigid-top dynamics.

Quick example: python run.py --samples 128 --workers 4
Paper batches: python run.py --paper-task 0 --workers 8
"""
import argparse
from concurrent.futures import ProcessPoolExecutor
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

import numpy as np
from scipy.integrate import solve_ivp
from icats.constants import amu2au,au2ang,ev2au
import icats
from model import MU,CM,CONSTANTS,INERTIA,rhs,invariants

HERE=Path(__file__).resolve().parent
ENERGY=93.8265
SEED=2026111017
CASES=((3,2,0.),(3,2,.8),(5,4,0.),(5,4,.8))


def inputs(folder,j,k,width,lo,hi,n,seed):
    # Masses and geometry reproduce exactly the benchmark's rigid Hamiltonian.
    m=np.array([14.003074,1.007825,1.007825,1.007825])*1822.888486209
    mhe=MU*m.sum()/(m.sum()-MU)
    radius=np.sqrt(INERTIA[2]/(3*m[1]))
    height=np.sqrt((INERTIA[0]-.5*INERTIA[2])/(m[0]*(3*m[1])/m.sum()))
    angle=np.arange(3)*2*np.pi/3
    ref=np.vstack(([0,0,height],np.column_stack((radius*np.cos(angle),radius*np.sin(angle),np.zeros(3)))))
    ref-=np.average(ref,axis=0,weights=m)
    xyz='4\nModel geometry matched to the MOLSCAT rotational constants\n'
    xyz+=''.join(f'{e} '+ ' '.join(f'{a:.16g}' for a in r*au2ang)+'\n' for e,r in zip(['N','H','H','H'],ref))
    (folder/'nh3.xyz').write_text(xyz)
    (folder/'he.xyz').write_text('1\nHe\nHe 0 0 0\n')
    (folder/'nh3.txt').write_text('xyz = nh3.xyz\nmasses = '+' '.join(f'{x/amu2au:.16g}' for x in m)+
        f'\nrotor-state = {j} {k}\nrotor-k-rule = nh3-para\nrotor-k-width = {width}\n')
    (folder/'he.txt').write_text(f'xyz = he.xyz\nmasses = {mhe/amu2au:.16g}\n')
    (folder/'input.txt').write_text(f'''mol = 0 he.txt
mol = 1 nh3.txt
Nsamp = {n}
workers = 1
seed = {seed}
seed-mode = fixed
run-mode = fresh
continue = False
Trot = 0
Tvib = 0
vib-mode = rigid
maxv = 0
collision-energy = {ENERGY/CM/ev2au:.16g}
Rz = {60*au2ang:.16g}
orbital-sampling = flat-b
minb = {lo*au2ang:.16g}
maxb = {hi*au2ang:.16g}
maxl = 0
wang = False
impact-phi = 0
output-frame = internal
orientation-mode = isotropic
audit-initial-sample = True
fileout = out
dirout = outputs
printout = 0 0 1 0
output-format = npz
units-out = au
plothist = False
hist_initial = False
hist_sampled = False
keepinfo = False
save-frequency = 0
progress = quiet
''')
    return ref,m,mhe


def prepare(folder,j,k,width,lo,hi,n,seed):
    folder.mkdir(parents=True,exist_ok=True)
    sources=[Path(__file__),HERE/'model.py']+sorted(Path(icats.__file__).parent.glob('*.py'))
    settings=dict(j=j,k=k,width=width,bmin=lo,bmax=hi,samples=n,seed=seed,
                  source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sources})
    manifest=folder/'manifest.json'
    if manifest.exists() and json.loads(manifest.read_text())!=settings:
        raise ValueError('Existing batch has different inputs or source; use a new output directory')
    manifest.write_text(json.dumps(settings,indent=2)+'\n')
    target=folder/'initial.npz'
    if target.exists():
        with np.load(target) as z:
            assert z['states'].shape==(n,18)
            return z['states'],z['b']
    ref,m,mhe=inputs(folder,j,k,width,lo,hi,n,seed)
    env=os.environ.copy()
    env['PYTHONPATH']=str(Path(icats.__file__).resolve().parent.parent)+os.pathsep+env.get('PYTHONPATH','')
    with (folder/'generation.log').open('w') as log:
        subprocess.run([sys.executable,'-m','icats.cli_init','input.txt'],cwd=folder,
                       env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
    exports=sorted((folder/'rd_input'/'outputs').glob('out_*.npz'),key=lambda p:int(p.stem.split('_')[-1]))
    assert len(exports)==n, (len(exports),n)
    states=[];bs=[];geometry_error=0.
    for path in exports:
        with np.load(path) as z:x,v=z['x'],z['v']
        com=np.average(x[1:],axis=0,weights=m)
        vel=np.average(v[1:],axis=0,weights=m)
        xx=x[1:]-com;vv=v[1:]-vel
        u,_,vt=np.linalg.svd(ref.T@(m[:,None]*xx))
        q=(u@vt).T
        assert np.linalg.det(q)>0
        geometry_error=max(geometry_error,float(abs(xx-ref@q.T).max()))
        ell=q.T@np.sum(m[:,None]*np.cross(xx,vv),axis=0)
        r=x[0]-com;p=MU*(v[0]-vel)
        b=np.linalg.norm(np.cross(r,p))/np.linalg.norm(p)
        assert lo-1e-7<=b<=hi+1e-7
        assert r@p<0
        assert abs(ell@ell-j*(j+1))<1e-7
        assert abs(p@p/(2*MU)*CM-ENERGY)<1e-6
        if width==0:assert abs(abs(ell[2])-k)<1e-7
        states.append(np.r_[r,p,ell,q.ravel()]);bs.append(b)
    assert geometry_error<1e-8
    states=np.asarray(states);bs=np.asarray(bs)
    np.savez_compressed(target,states=states,b=bs)
    with np.load(target) as z:assert np.array_equal(states,z['states'])
    (folder/'initial_audit.json').write_text(json.dumps(dict(samples=n,
        geometry_error=geometry_error,mean_K2=float(np.mean(states[:,8]**2)),target_K2=k*k,
        bmin=float(bs.min()),bmax=float(bs.max()),model_reduced_mass_au=MU),indent=2)+'\n')
    # Full-precision exports are now retained compactly as initial.npz.
    for path in exports:path.unlink()
    for pattern in ('work_sys_*.pkl','dat_*.pkl','info_*.pkl'):
        for path in folder.rglob(pattern):path.unlink()
    return states,bs


def propagate(args):
    index,state,b,width,lo,hi=args
    e0,j0,_=invariants(state)
    rstart=np.linalg.norm(state[:3]);speed=np.linalg.norm(state[3:6])/MU
    def outgoing(t,y):return np.linalg.norm(y[:3])-rstart
    outgoing.terminal=True;outgoing.direction=1
    start=time.perf_counter()
    for rtol in (2e-9,2e-10,2e-11):
        for flight in (5.,20.,80.):
            sol=solve_ivp(rhs,(0,flight*rstart/speed),state,method='DOP853',
                          rtol=rtol,atol=rtol*.01,events=outgoing)
            if not sol.success:raise RuntimeError(sol.message)
            if len(sol.t_events[0]):break
        else:raise RuntimeError('Outgoing boundary not reached')
        val=[invariants(y) for y in sol.y.T]
        de=max(abs(e-e0) for e,_,_ in val)
        dj=max(np.linalg.norm(j-j0) for _,j,_ in val)
        orth=max(o for _,_,o in val)
        if de<1e-8 and dj<1e-5 and orth<1e-6:break
    else:raise RuntimeError('Conservation tolerances not met')
    y=sol.y[:,-1];ell=y[6:9]
    return dict(sample=index,width=width,b_bohr=b,initial_K=float(state[8]),
        initial_Erot_cm=float(CONSTANTS@state[6:9]**2*CM),final_K=float(ell[2]),
        final_abs_K=float(abs(ell[2])),final_j_cont=float((np.sqrt(1+4*(ell@ell))-1)/2),
        final_Erot_cm=float(CONSTANTS@ell**2*CM),final_R_bohr=float(np.linalg.norm(y[:3])),
        max_energy_error_Eh=de,max_J_error_hbar=dj,max_orthogonality_error=orth,
        seconds=time.perf_counter()-start,nfev=sol.nfev,rtol=rtol,
        area_weight_A2=2*np.pi*b*(hi-lo)*(.529177210903**2))


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output',type=Path,default=Path('results'))
    ap.add_argument('--samples',type=int,default=128)
    ap.add_argument('--workers',type=int,default=1)
    ap.add_argument('--seed',type=int,default=SEED)
    ap.add_argument('--j',type=int,default=3);ap.add_argument('--k',type=int,default=2)
    ap.add_argument('--width',type=float,default=.8)
    ap.add_argument('--bmin',type=float,default=0.);ap.add_argument('--bmax',type=float,default=16.)
    ap.add_argument('--paper-task',type=int,choices=range(256))
    ap.add_argument('--prepare-only',action='store_true')
    args=ap.parse_args()
    if args.samples<1 or args.workers<1:ap.error('Positive samples/workers required')
    if args.paper_task is not None:
        case,rest=divmod(args.paper_task,64);shell,batch=divmod(rest,4)
        args.j,args.k,args.width=CASES[case]
        args.bmin,args.bmax=shell,shell+1
        args.samples=4096
        args.seed=args.seed+(case//2)*10000+shell*100+batch
        folder=args.output/f'j{args.j}_k{args.k}_s{args.width:g}'/f'b{shell:02d}_part{batch:02d}'
    else:folder=args.output
    states,b=prepare(folder,args.j,args.k,args.width,args.bmin,args.bmax,args.samples,args.seed)
    if args.prepare_only:return
    target=folder/'final.csv.gz'
    if target.exists():
        with gzip.open(target,'rt') as stream:assert sum(1 for _ in csv.DictReader(stream))==args.samples
        print('Verified existing result',target);return
    jobs=[(i,y,bb,args.width,args.bmin,args.bmax) for i,(y,bb) in enumerate(zip(states,b))]
    with ProcessPoolExecutor(max_workers=args.workers) as pool:rows=list(pool.map(propagate,jobs,chunksize=8))
    temporary=target.with_suffix('.tmp')
    with gzip.open(temporary,'wt') as stream:
        writer=csv.DictWriter(stream,fieldnames=rows[0]);writer.writeheader();writer.writerows(rows)
    temporary.replace(target)
    print(f'Completed {target}: N={len(rows)}, trajectory seconds={sum(r["seconds"] for r in rows):.1f}',flush=True)


if __name__=='__main__':main()
