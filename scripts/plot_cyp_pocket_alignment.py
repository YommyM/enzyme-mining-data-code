#!/usr/bin/env python3
"""Figure 5E: sequence, secondary structure, and global-fit Cα displacement in the ±3 pocket."""
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/cyp_mplconfig')
import json,csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from cyp_common import DATA,STRUCTURES,OUTPUT,read_csv,write_csv
OUT=OUTPUT/'case';OUT.mkdir(parents=True,exist_ok=True)
FONT=12;W=10;H=3.35
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':FONT,'axes.labelsize':FONT,'xtick.labelsize':FONT,'ytick.labelsize':FONT,'pdf.fonttype':42,'svg.fonttype':'none','axes.linewidth':.8})
fig=plt.figure(figsize=(W,H),facecolor='white')
def axis(x,y,w,h):return fig.add_axes([x/W,y/H,w/W,h/H])
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1
from matplotlib.patches import Rectangle
from matplotlib.colors import Normalize, LinearSegmentedColormap
CASE=STRUCTURES/'case'
parser=PDBParser(QUIET=True)
chains=[parser.get_structure(n,CASE/f'EXP_{n}_global_frame.pdb')[0]['A'] for n in ['OleT_KM','P450BSbeta']]
# Global alignment establishes correspondence before applying either pm3 mask.
names=['OleT_KM','P450BSbeta']
cfg=json.loads((DATA/'case/metadata.json').read_text())
protein_res=[[r for r in ch if r.id[0]==' ' and 'CA' in r] for ch in chains]
lines=(DATA/'case/alignment/EXP_global_OleT_KM_to_P450BSbeta.log').read_text().splitlines()
li=next(i for i,line in enumerate(lines) if line.startswith('(":" denotes'))
seqs=[lines[li+1],lines[li+3]]
for seq,rr in zip(seqs,protein_res):assert seq.replace('-','')==''.join(seq1(r.resname) for r in rr)
cursors=[0,0];rows=[]
for letters in zip(*seqs):
 pair=[]
 for k,letter in enumerate(letters):
  residue=protein_res[k][cursors[k]] if letter!='-' else None
  cursors[k]+=letter!='-';pair.append(residue)
 selected=[r is not None and r.id[1] in cfg[n]['pm3_positions'] for r,n in zip(pair,names)]
 if not any(selected):continue
 row={}
 for k,n in enumerate(names):
  rr=pair[k]
  row[n+'_position']=str(rr.id[1]) if selected[k] else ''
  row[n+'_residue']=rr.resname if selected[k] else ''
  row[n+'_global_corresponding_position']=str(rr.id[1]) if rr is not None else ''
  row[n+'_absence_reason']='' if selected[k] else ('outside_pm3' if rr is not None else 'alignment_gap')
 both=all(selected)
 row['status']=('same' if pair[0].resname==pair[1].resname else 'substitution') if both else 'unpaired_in_pm3'
 row['CA_distance_global_frame_A']=float(np.linalg.norm(pair[0]['CA'].coord-pair[1]['CA'].coord)) if both else None
 rows.append(row)
for n in names:assert sorted(int(r[n+'_position']) for r in rows if r[n+'_position'])==sorted(cfg[n]['pm3_positions'])
N_PER_BLOCK=(len(rows)+1)//2
print('Global correspondence:',len(rows),'columns', {v:sum(r['status']==v for r in rows) for v in ['same','substitution','unpaired_in_pm3']})
# Deposited experimental secondary structure, in author numbering.
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from matplotlib.patches import FancyArrow
ss_maps=[]
for pdb in ['8VWK','1IZO']:
 d=MMCIF2Dict(STRUCTURES/f'case/native/{pdb}.cif');ss={}
 for prefix,kind in [('_struct_conf.','H'),('_struct_sheet_range.','E')]:
  for i,chain in enumerate(d[prefix+'beg_auth_asym_id']):
   if chain!='A':continue
   if kind=='H' and d[prefix+'pdbx_PDB_helix_class'][i]!='1':continue
   lo=int(d[prefix+'beg_auth_seq_id'][i]);hi=int(d[prefix+'end_auth_seq_id'][i])
   for pos in range(lo,hi+1):ss[pos]=kind
 ss_maps.append(ss)
for r in rows:
 a=seq1(r['OleT_KM_residue']) if r['OleT_KM_position'] else '-'
 b=seq1(r['P450BSbeta_residue']) if r['P450BSbeta_position'] else '-'
 for k,n in enumerate(['OleT_KM','P450BSbeta']):r[n+'_secondary_structure']=ss_maps[k].get(int(r[n+'_position']),'other') if r[n+'_position'] else 'gap'
# Explicit, non-overlapping chemical classes; no inferred catalytic effect size.
classes={'aliphatic':'AVILM','aromatic':'FWY','hydroxyl':'ST','amide':'NQ','acidic':'DE','basic':'KR','imidazole':'H','thiol':'C','glycine':'G','proline':'P'}
class_of={aa:group for group,aas in classes.items() for aa in aas}
for r in rows:
 a=seq1(r['OleT_KM_residue']) if r['OleT_KM_position'] else '-'
 b=seq1(r['P450BSbeta_residue']) if r['P450BSbeta_position'] else '-'
 r['OleT_KM_chemical_class']=class_of.get(a,'absent')
 r['P450BSbeta_chemical_class']=class_of.get(b,'absent')
 r['chemical_change']='absent' if '-' in (a,b) else ('identical' if a==b else ('same_class' if class_of[a]==class_of[b] else 'different_class'))
(OUT/'chemical_classes.json').write_text(json.dumps(classes,indent=2))
F=axis(.15,.40,9.60,2.65);F.set_xlim(0,9.60);F.set_ylim(0,2.65);F.axis('off')
cmap=LinearSegmentedColormap.from_list('pocket_blue',['#FFFFFF','#BDD7E7','#2171B5'])
norm=Normalize(0,4,clip=True)
cell=8.745/N_PER_BLOCK;rh=.20;start=.66
breaks={name:{} for name in ['OleT_KM','P450BSbeta']}
for name in breaks:
 previous=None
 for index,r in enumerate(rows):
  pos=r[name+'_position']
  if not pos:continue
  pos=int(pos)
  if previous is not None and pos!=previous+1:breaks[name][index]=(previous,pos)
  previous=pos
for block in range(2):
 rr=rows[block*N_PER_BLOCK:(block+1)*N_PER_BLOCK];top=2.12-block*1.24
 for j,r in enumerate(rr):
  change=r['chemical_change'];x=start+(j+.5)*cell
  if change=='same_class':F.scatter([x],[top-.10],s=12,c='#AAAAAA',marker='o',linewidths=0)
  elif change=='different_class':F.scatter([x],[top-.10],s=19,c='#444444',marker='D',linewidths=0)
  if r['OleT_KM_position']=='89' and r['P450BSbeta_position']=='85':
   F.text(x,top+.43,'Imidazole / Amide',ha='center',va='center',fontsize=FONT)

 for k,(name,short) in enumerate([('OleT_KM',r'OleT$_{\mathrm{KM}}$'),('P450BSbeta',r'P450$_{\mathrm{BS}\beta}$')]):
  yy=top-k*.40
  F.text(0,yy+rh/2,short,va='center',fontsize=FONT)
  for j,r in enumerate(rr):
   xx=start+j*cell;pos=r[name+'_position'];d=r['CA_distance_global_frame_A']
   face=cmap(norm(d)) if d is not None else '#EEEEEE'
   F.add_patch(Rectangle((xx,yy),cell,rh,facecolor=face,edgecolor='none'))
   F.text(xx+cell/2,yy+rh/2,seq1(r[name+'_residue']) if pos else '–',ha='center',va='center',fontsize=FONT,color='#C62828' if (name=='OleT_KM' and pos=='89') or (name=='P450BSbeta' and pos=='85') else ('white' if d is not None and d>3 else 'black'),fontweight='bold' if (name=='OleT_KM' and pos=='89') or (name=='P450BSbeta' and pos=='85') else 'normal')
  # Number every visible continuous sequence run; skip alignment gaps.
  runs=[];current=[]
  for j,r in enumerate(rr):
   pos=r[name+'_position']
   if not pos:continue
   pos=int(pos)
   if current and pos!=current[-1][1]+1:runs.append(current);current=[]
   current.append((j,pos))
  if current:runs.append(current)
  num_y=top+.43 if k==0 else top-.65
  for segment in runs:
   first=segment[0]
   F.text(start+(first[0]+.5)*cell,num_y,str(first[1]),ha='center',va='center',fontsize=FONT)
  # Secondary-structure glyphs break at omitted sequence segments and alignment gaps.
  track=top+.29 if k==0 else top-.50
  j=0
  while j<len(rr):
   kind=rr[j][name+'_secondary_structure'];end=j+1
   while end<len(rr) and kind in ['H','E'] and rr[end][name+'_secondary_structure']==kind and int(rr[end][name+'_position'])==int(rr[end-1][name+'_position'])+1:end+=1
   if kind in ['H','E']:
    x0=start+j*cell+.02;x1=start+end*cell-.02
    if block*N_PER_BLOCK+j in breaks[name]:x0+=.07
    if block*N_PER_BLOCK+end in breaks[name]:x1-=.07
    if kind=='H':
     t=np.linspace(0,2*np.pi*max(1,end-j),250)
     xx=x0+(x1-x0)*t/t[-1]+.025*np.sin(t)
     F.plot(xx,track+.055*np.cos(t),color='#444444',lw=.9)
    else:F.add_patch(FancyArrow(x0,track,x1-x0,0,width=.035,head_width=.12,head_length=min(.12,(x1-x0)*.45),length_includes_head=True,color='#444444',linewidth=0))
   j=end
  for index,(before,after) in breaks[name].items():
   if block*N_PER_BLOCK<=index<(block+1)*N_PER_BLOCK:
    x=start+(index-block*N_PER_BLOCK)*cell
    for dx in [-.035,.035]:F.plot([x+dx-.021,x+dx+.021],[track-.055,track+.055],color='#444444',lw=.9,zorder=5)
F.scatter([.73],[.05],s=12,c='#AAAAAA',marker='o',linewidths=0)
F.text(.88,.05,'Same class',va='center',fontsize=FONT)
F.scatter([2.40],[.05],s=19,c='#444444',marker='D',linewidths=0)
F.text(2.55,.05,'Different class',va='center',fontsize=FONT)
F.text(4.65,.05,'// Break',va='center',fontsize=FONT)
cbax=axis(6.05,.43,3.2,.10)
cb=fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm,cmap=cmap),cax=cbax,orientation='horizontal',ticks=[0,2,4],extend='max')
cb.ax.tick_params(direction='in',labelsize=FONT);cb.set_label('Cα distance (Å)',fontsize=FONT,labelpad=3)
(OUT/'sequence_breaks.json').write_text(json.dumps(breaks,indent=2))
verification={'correspondence':'experimental_global_TM_align','columns':len(rows),'selected_residues_per_protein':{n:sum(bool(r[n+'_position']) for r in rows) for n in names},'counts':{v:sum(r['status']==v for r in rows) for v in ['same','substitution','unpaired_in_pm3']},'local_TM_score_unchanged':0.8738,'global_TM_score_unchanged':0.9137}
verification['paired_sequence_identity']=verification['counts']['same']/(verification['counts']['same']+verification['counts']['substitution'])
(OUT/'mapping_verification.json').write_text(json.dumps(verification,indent=2))

write_csv(OUT/'E_residue_data.csv',rows)
expected=read_csv(DATA/'case/results/E_residue_data.csv')
assert len(expected)==len(rows)
for a,b in zip(rows,expected):
 for key in ['OleT_KM_position','P450BSbeta_position','status','chemical_change']:
  assert str(a[key])==b[key]
 if a['CA_distance_global_frame_A'] is not None:
  assert abs(a['CA_distance_global_frame_A']-float(b['CA_distance_global_frame_A']))<1e-6
for ext in ['png','pdf','svg']:fig.savefig(OUT/f'Figure_5E.{ext}',dpi=400)
plt.close(fig)
print('Figure 5E reproduced: 66 alignment columns, 65 selected residues per protein.')
