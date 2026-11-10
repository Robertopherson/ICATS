"""Plot completed tutorial trajectories against bundled quantum tables."""
import argparse
import csv
import gzip
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

HERE=Path(__file__).resolve().parent
A2=.529177210903**2


def estimate(values):
    return float(np.mean(values)),float(np.std(values,ddof=1)/np.sqrt(len(values))) if len(values)>1 else np.nan


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('inputs',nargs='+',type=Path)
    ap.add_argument('--output',type=Path,default=Path('plots'))
    ap.add_argument('--no-energy-filter',action='store_true')
    ap.add_argument('--no-symmetry-filter',action='store_true')
    args=ap.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    cases={};coverage=[]
    paths=sorted({p.resolve() for root in args.inputs for p in root.rglob('final.csv.gz')})
    if not paths:raise ValueError('No completed final.csv.gz files found')
    for path in paths:
        manifest=json.loads((path.parent/'manifest.json').read_text())
        key=manifest['j'],manifest['k'],manifest['width']
        with gzip.open(path,'rt') as stream:data=np.atleast_1d(np.genfromtxt(stream,delimiter=',',names=True))
        if len(data)!=manifest['samples']:raise ValueError(f'Incomplete file {path}')
        if len(np.unique(data['sample']))!=len(data):raise ValueError(f'Duplicate sample IDs: {path}')
        if not all(np.isfinite(data[name]).all() for name in data.dtype.names):raise ValueError('Nonfinite result')
        assert data['max_energy_error_Eh'].max()<1e-8
        assert data['max_J_error_hbar'].max()<1e-5
        assert data['max_orthogonality_error'].max()<1e-6
        lo,hi=manifest['bmin'],manifest['bmax']
        cases.setdefault(key,{}).setdefault((lo,hi),[]).append(data)
        coverage.append(dict(path=str(path),initial_j=key[0],initial_K=key[1],width=key[2],
                             bmin=lo,bmax=hi,n=len(data),seed=manifest['seed']))
    seen=set()
    for c in coverage:
        key=tuple(c[k] for k in ('initial_j','initial_K','width','bmin','bmax','seed'))
        if key in seen:raise ValueError('Repeated seeds in the same preparation interval')
        seen.add(key)
    op=[];sections=[]
    for (j,k,s),strata in cases.items():
        intervals=sorted(strata)
        if any(b[0]<a[1] for a,b in zip(intervals,intervals[1:])):
            raise ValueError('Overlapping preparation intervals cannot be added as independent strata')
        partial=intervals[0][0]!=0 or intervals[-1][1]<16 or any(b[0]!=a[1] for a,b in zip(intervals,intervals[1:]))
        if partial:print(f'WARNING: incomplete b=0..16 coverage for {(j,k,s)}: {intervals}')
        sums={jf:[] for jf in range(1,7) if jf!=j}
        for (lo,hi),batches in strata.items():
            d=np.concatenate(batches);b=d['b_bohr']
            jf=np.floor(d['final_j_cont']+.5).astype(int)
            kf=np.floor(d['final_K']+.5).astype(int)
            allowed=np.ones(len(d),dtype=bool)
            if not args.no_symmetry_filter:allowed &= (abs(kf)<=jf)&(kf%3!=0)
            if not args.no_energy_filter:
                e=9.9441*jf*(jf+1)+(6.2294-9.9441)*kf*kf
                allowed &= e <= d['initial_Erot_cm']+93.8265
            for target in sums:
                event=(jf==target)&allowed
                sigma,se=estimate(d['area_weight_A2']*event)
                sums[target].append((sigma,se))
                edges=np.unique(np.r_[lo,np.arange(np.ceil(lo/.5)*.5,hi,.5),hi])
                for a,c in zip(edges[:-1],edges[1:]):
                    mask=(b>=a)&(b<c)
                    if not mask.any():continue
                    p,err=estimate(2*b[mask]/(a+c)*event[mask])
                    op.append(dict(initial_j=j,initial_K=k,width=s,final_j=target,b=(a+c)/2,
                                   opacity=p,se=err,n=int(mask.sum())))
        for target,est in sums.items():
            sections.append(dict(initial_j=j,initial_K=k,width=s,final_j=target,
                sigma_A2=sum(x[0] for x in est),se_A2=np.sqrt(sum(x[1]**2 for x in est)),partial=partial))
    def write(name,rows):
        with (args.output/name).open('w') as stream:
            writer=csv.DictWriter(stream,fieldnames=rows[0]);writer.writeheader();writer.writerows(rows)
    write('coverage.csv',coverage);write('opacities.csv',op);write('cross_sections.csv',sections)
    q=np.genfromtxt(HERE/'quantum_opacities.csv',delimiter=',',names=True)
    qi=np.genfromtxt(HERE/'quantum_cross_sections.csv',delimiter=',',names=True)
    states=sorted({(j,k) for j,k,s in cases})
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axs=plt.subplots(len(states),4,figsize=(13,3.2*len(states)+1),squeeze=False,
                         gridspec_kw={'width_ratios':[1,1,1,1.25]})
    ymax=0.;cmax=0.
    for row,(j,k) in enumerate(states):
        if (j,k) not in ((3,2),(5,4)):raise ValueError('Bundled reference covers (3,2) and (5,4) only')
        targets=(1,2,4) if j==3 else (3,4,6)
        widths=sorted(s for jj,kk,s in cases if (jj,kk)==(j,k))
        colors={0.:'#177da3',.8:'#bd4c39',.4:'#95629c',.2:'#568348'}
        for col,target in enumerate(targets):
            ax=axs[row,col]
            qq=q[(q['initial_j']==j)&(q['j_final']==target)]
            ax.plot(qq['b'],qq['p'],'k-',lw=1.5,label='Quantum')
            ymax=max(ymax,float(qq['p'].max()))
            for s in widths:
                g=sorted([r for r in op if r['initial_j']==j and r['width']==s and r['final_j']==target],key=lambda r:r['b'])
                xx,yy,ee=([r[key] for r in g] for key in ('b','opacity','se'))
                ax.errorbar(xx,yy,yerr=ee,color=colors.get(s),marker='o',ms=2.5,capsize=2,lw=1,label=f'Classical s={s:g}')
                ymax=max(ymax,float(np.nanmax(np.array(yy)+np.array(ee))))
            ax.text(.96,.94,f"j'={target}",transform=ax.transAxes,ha='right',va='top')
            ax.set_xlim(0,8);ax.set_xticks([0,2,4,6] if col<2 else [0,2,4,6,8]);ax.grid(alpha=.2)
            if col:ax.tick_params(labelleft=False)
            else:ax.set_ylabel('Opacity')
            if row==len(states)-1:ax.set_xlabel('b (bohr)')
            else:ax.tick_params(labelbottom=False)
        axs[row,0].set_title(f'Incoming (j, |K|) = ({j}, {k})',loc='left',fontsize=11,pad=12)
        ax=axs[row,3];js=[n for n in range(1,7) if n!=j]
        qq=[float(qi[(qi['initial_j']==j)&(qi['final_j']==jf)]['sigma_A2'].sum()) for jf in js]
        ax.plot(js,qq,'k_',ms=13,mew=2);cmax=max(cmax,max(qq))
        for i,s in enumerate(widths):
            g=[next(r for r in sections if r['initial_j']==j and r['width']==s and r['final_j']==jf) for jf in js]
            yy=[r['sigma_A2'] for r in g];ee=[r['se_A2'] for r in g]
            ax.errorbar(np.array(js)+(i-(len(widths)-1)/2)*.14,yy,yerr=ee,fmt='o',ms=3,capsize=3,color=colors.get(s))
            cmax=max(cmax,float(np.nanmax(np.array(yy)+np.array(ee))))
        ax.set_xticks(js);ax.set_xlabel("Final j'");ax.set_ylabel(r'Cross section ($\AA^2$)');ax.grid(axis='y',alpha=.2)
    for ax in axs[:,:3].flat:ax.set_ylim(0,.05*np.ceil(1.05*ymax/.05))
    for ax in axs[:,3]:ax.set_ylim(0,1.12*cmax)
    handles,labels=axs[0,0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='upper center',bbox_to_anchor=(.5,.94),ncol=3,frameon=False)
    fig.suptitle('He + NH3: ICATS initial conditions and rigid-model dynamics')
    filters=('para allowed; ' if not args.no_symmetry_filter else 'all nearest bins; ')
    filters+=('energy-open filter' if not args.no_energy_filter else 'no energy filter')
    fig.text(.5,.025,filters+'. All incoming trajectories in denominator; error bars: 1 SE.\n'
             'Cross sections use the full sampled b range. See coverage.csv for completed batches.',ha='center',fontsize=9)
    fig.tight_layout(rect=(0,.09,1,.86))
    for row in range(len(states)):
        left=axs[row,0].get_position().x0;right=axs[row,2].get_position().x1
        width=(right-left-.012)/3
        for col in range(3):
            ax=axs[row,col];p=ax.get_position();ax.set_position([left+col*(width+.006),p.y0,width,p.height])
    fig.savefig(args.output/'opacities.pdf');fig.savefig(args.output/'opacities.png',dpi=160)
    print('Wrote',args.output/'opacities.pdf')


if __name__=='__main__':main()
