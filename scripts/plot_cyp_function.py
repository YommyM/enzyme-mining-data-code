#!/usr/bin/env python3
"""Plot individual Figure 5A/B/C panels from the deposited final analysis tables."""
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/cyp_mplconfig')
import argparse,json
from pathlib import Path
from collections import Counter
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from cyp_common import DATA,OUTPUT,read_csv
BLUE='#0072B2';GRAY='#808080'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':12,'axes.labelsize':12,'xtick.labelsize':12,'ytick.labelsize':12,'legend.fontsize':12,'pdf.fonttype':42,'svg.fonttype':'none','axes.linewidth':.8,'xtick.major.width':.8,'ytick.major.width':.8})

def finish(fig,ax,letter,layout=True):
    ax.tick_params(direction='in');ax.spines[['top','right']].set_visible(False)
    if layout: fig.tight_layout()
    OUTPUT.mkdir(parents=True,exist_ok=True)
    for ext in ['png','pdf','svg']:fig.savefig(OUTPUT/f'Figure_5{letter}.{ext}',dpi=400)
    plt.close(fig)

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--results',type=Path,default=DATA/'results');a=p.parse_args();data=a.results
    fig,ax=plt.subplots(figsize=(4.4,3.2));rows=read_csv(data/'A_ARI.csv')
    for rep,col,marker in [('global',GRAY,'o'),('local',BLUE,'s')]:
        rr=[r for r in rows if r['representation']==rep]
        ax.plot([int(r['clusters']) for r in rr],[float(r['ARI']) for r in rr],c=col,marker=marker,lw=1.4,ms=3.5,label='Global' if rep=='global' else 'Local ±3')
    ax.axvline(7,color='.72',ls=(0,(3,3)),lw=.8);ax.set(xlabel='Clusters',ylabel='Adjusted Rand index',xticks=range(2,11),ylim=(0,.55),yticks=np.arange(0,.6,.1));ax.legend(frameon=False);finish(fig,ax,'A')
    rows=[r for r in read_csv(data/'B_class_AUC.csv') if r['condition']=='all_references']
    short=['Hydroxylation','Coupling/cycl.','Dealkylation','Decarboxylation','Methyl → COOH','Nitration','Rearrangement']
    fig=plt.figure(figsize=(4.0,2.9));ax=fig.add_axes([.52,.30,.435,.58]);lab=fig.add_axes([.055,.30,.45,.58])
    lab.set(xlim=(0,1),ylim=(6.5,-.5));lab.axis('off')
    for j,(r,name) in enumerate(zip(rows,short)):
        if j%2==0:
            ax.axhspan(j-.46,j+.46,color='#F3F6F8',lw=0);lab.axhspan(j-.46,j+.46,color='#F3F6F8',lw=0)
        lab.text(0,j,name,va='center');lab.text(.98,j,str(r['queries']),ha='right',va='center')
        g,l=float(r['global_AUC']),float(r['local_AUC']);ax.plot([g,l],[j,j],color='.7',lw=1)
        ax.scatter(g,j,s=20,c=GRAY,marker='o',label='Global' if j==0 else None);ax.scatter(l,j,s=23,c=BLUE,marker='s',label='Local ±3' if j==0 else None)
    lab.text(.98,1.03,'n',transform=lab.transAxes,ha='right',va='bottom')
    ax.set(yticks=[],ylim=(6.5,-.5),xlim=(.45,1.04),xticks=[.5,.75,1],xlabel='Ranking AUC');ax.spines['left'].set_visible(False);ax.axvline(.5,color='.65',ls='--',lw=.8)
    handles,labels=ax.get_legend_handles_labels();fig.legend(handles,labels,loc='lower center',bbox_to_anchor=(.5,.01),frameon=False,ncol=2,handlelength=.7,handletextpad=.3)
    finish(fig,ax,'B',layout=False)
    rows=read_csv(data/'BC_query_AUC.csv');conditions=['all_references','exclude_cluster','exclude_cluster_features'];eligible={r['accession'] for r in rows if r['condition']==conditions[-1]}
    tests=json.loads((data/'C_significance.json').read_text())['results'];rng=np.random.default_rng(20261117);jitter={a:rng.uniform(.055,.145) for a in sorted(eligible)}
    fig,ax=plt.subplots(figsize=(4.7,3.6))
    for h,cond in enumerate(conditions):
        rr=sorted([r for r in rows if r['condition']==cond and r['accession'] in eligible],key=lambda r:r['accession']);counts=Counter(r['category'] for r in rr);weights=np.array([1/len(counts)/counts[r['category']] for r in rr])
        for rep,col,offset in [('global',GRAY,-.20),('local',BLUE,.20)]:
            vals=np.array([float(r[rep+'_AUC']) for r in rr]);x=h+offset;order=np.argsort(vals);v=vals[order];w=weights[order];q1,med,q3=np.interp([.25,.5,.75],np.cumsum(w)-w/2,v)
            lo=min(q1,np.min(v[v>=q1-.5*(q3-q1)]));hi=max(q3,np.max(v[v<=q3+.5*(q3-q1)]))
            ax.plot([x+.013]*2,[lo,hi],color=col,lw=.8)
            ax.hlines([lo,hi],x-.014,x+.040,color=col,lw=.8)
            ax.add_patch(Rectangle((x-.017,q1),.060,q3-q1,facecolor='white',edgecolor=col,lw=.8,zorder=4));ax.plot([x-.017,x+.043],[med]*2,color=col,lw=1,zorder=5)
            ax.scatter([x+jitter[r['accession']] for r in rr],vals,s=7,c=col,alpha=.48,linewidths=0)
            ax.scatter(x+.013,weights@vals,s=22,marker='D',c=col,edgecolor='white',linewidth=.4,zorder=6)
        pv=tests[h]['p_Holm'];symbol='***' if pv<.001 else ('**' if pv<.01 else ('*' if pv<.05 else 'ns'))
        ax.plot([h-.2,h-.2,h+.2,h+.2],[1.035,1.055,1.055,1.035],c='.3',lw=.8);ax.text(h,1.065,symbol,ha='center',va='bottom')
    ax.axhline(.5,color='#B8BEC4',ls='--',lw=.8);ax.set(xlim=(-.42,2.43),ylim=(-.025,1.18),xticks=[0,1,2],xticklabels=['All','External\ncandidates','External\nprofiles'],yticks=[0,.5,1],ylabel='Query AUC');finish(fig,ax,'C')
    print('Wrote individual Figure 5A, 5B and 5C panels; no manuscript assembly performed.')

if __name__=='__main__':main()
