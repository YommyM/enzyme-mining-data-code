#!/usr/bin/env python3
"""Figure 5D–F case: ligand/heme positioning, pocket displacements and deposited water geometry."""
from pathlib import Path
import json,csv,hashlib
import numpy as np
from Bio.PDB import PDBParser,MMCIFParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB.vectors import Vector,calc_dihedral
from cyp_common import DATA,STRUCTURES,OUTPUT
OUT=OUTPUT/'case/geometry';OUT.mkdir(parents=True,exist_ok=True)
CASE=STRUCTURES/'case'
CFG=json.loads((DATA/'case/metadata.json').read_text());names=['OleT_KM','P450BSbeta'];pdbs=['8VWK','1IZO']
def dist(a,b):return float(np.linalg.norm(a.coord-b.coord))
def angle(u,v):return float(np.degrees(np.arccos(np.clip(np.dot(u,v)/np.linalg.norm(u)/np.linalg.norm(v),-1,1))))
def heavy(res):return [a for a in res if a.element not in ['H','D']]
def side(res):return [a for a in heavy(res) if a.name not in ['N','CA','C','O','OXT']]
def near(aa,bb):
 x=min((dist(a,b),a.name,b.name) for a in aa for b in bb)
 return {'distance_A':x[0],'protein_atom':x[1],'target_atom':x[2]}
def save(name,rows):
 with (OUT/name).open('w',encoding='utf-8-sig',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
models=[];geometry=[];keyrows=[];contacts=[];argrows=[];waters=[];keypolar=[]
for n,pdb,positions in zip(names,pdbs,[[89,244,173,174,175,176,303,83],[85,242,170,171,172,173,292,79]]):
 c=CFG[n];ch=PDBParser(QUIET=True).get_structure(n,CASE/f'EXP_{n}_global_frame.pdb')[0]['A']
 lig=next(r for r in ch if r.id[1]==c['ligand_residue']);heme=next(r for r in ch if r.id[1]==c['heme_residue']);fe=heme['FE']
 d=MMCIF2Dict(STRUCTURES/f'case/native/{pdb}.cif')
 bonds=[(d['_chem_comp_bond.atom_id_1'][i],d['_chem_comp_bond.atom_id_2'][i]) for i,v in enumerate(d['_chem_comp_bond.comp_id']) if v==c['ligand']]
 assert all((a,b) in bonds or (b,a) in bonds for a,b in [('C1','O1'),('C1','O2'),('C1','C2'),('C2','C3')])
 ns=np.array([heme[a].coord for a in ['NA','NB','NC','ND']],float);center=ns.mean(axis=0);_,_,vh=np.linalg.svd(ns-center);normal=vh[-1]
 if np.dot(lig['C3'].coord-fe.coord,normal)<0:normal=-normal
 xyz=np.array([lig[a].coord for a in ['C1','O1','O2']],float);_,_,cvh=np.linalg.svd(xyz-xyz.mean(axis=0))
 g={'protein':n,'pdb':pdb,'ligand':lig.resname,'plane_angle_carboxyl_vs_heme_deg':min(angle(normal,cvh[-1]),180-angle(normal,cvh[-1])),'C1_C2_C3_C4_dihedral_deg':float(np.degrees(calc_dihedral(*(Vector(lig[a].coord) for a in ['C1','C2','C3','C4']))))}
 for a in ['C1','C2','C3','O1','O2']:
  v=lig[a].coord-fe.coord;z=float(v@normal)
  g[a+'_Fe_distance_A']=dist(lig[a],fe);g[a+'_height_from_Fe_along_heme_normal_A']=z;g[a+'_radial_offset_A']=float(np.sqrt(max(0,np.dot(v,v)-z*z)));g[a+'_axis_angle_deg']=angle(v,normal)
 geometry.append(g)
 for pos in positions:
  r=ch[pos]
  for target,aa in [('whole_ligand',heavy(lig)),('carboxyl',[lig[a] for a in ['C1','O1','O2']]),('alpha_C2',[lig['C2']]),('beta_C3',[lig['C3']]),('heme',heavy(heme)),('iron',[fe])]:
   for mode,bb in [('sidechain',side(r)),('all_heavy',heavy(r))]:
    keyrows.append(dict(protein=n,position=pos,residue=r.resname,selection=mode,target=target,**near(bb,aa)))
  if pos==positions[0]:
   for a in side(r):
    if a.element not in ['N','O','S']:continue
    for b in heavy(heme):
     if b.element in ['N','O','S'] and dist(a,b)<=4.5:keypolar.append(dict(protein=n,position=pos,atom=a.name,heme_atom=b.name,distance_A=dist(a,b)))
 for r in ch:
  if r.id[0]!=' ' or 'CA' not in r:continue
  for target,aa in [('carboxyl',[lig[a] for a in ['C1','O1','O2']]),('alpha_C2',[lig['C2']]),('beta_C3',[lig['C3']])]:
   m=near(heavy(r),aa)
   if m['distance_A']<=4.5:contacts.append(dict(protein=n,position=r.id[1],residue=r.resname,target=target,in_pm3=r.id[1] in c['pm3_positions'],**m))
 arg=ch[positions[1]]
 for a in ['NE','NH1','NH2']:
  for b in ['O1','O2']:argrows.append(dict(protein=n,arg_position=positions[1],arginine_atom=a,ligand_atom=b,distance_A=dist(arg[a],lig[b])))
 # Use native coordinates for deposited waters (not present in display models).
 raw=MMCIFParser(QUIET=True).get_structure(pdb,STRUCTURES/f'case/native/{pdb}.cif')[0]
 rc=raw['A'];rl=next(r for r in rc if r.resname==c['ligand'] and r.id[1]==c['ligand_residue']);rh=next(r for r in rc if r.resname=='HEM' and r.id[1]==c['heme_residue']);kr=rc[positions[0]]
 for chain in raw:
  for water in chain:
   if water.resname not in ['HOH','WAT']:continue
   wa=heavy(water)[0];pol=[a for a in side(kr) if a.element in ['N','O']]
   m=near(pol,[wa])
   if m['distance_A']<=4.0:
    waters.append(dict(protein=n,water_chain=chain.id,water_number=water.id[1],key_atom=m['protein_atom'],key_water_distance_A=m['distance_A'],water_Fe_A=dist(wa,rh['FE']),water_carboxyl_O_A=min(dist(wa,rl[a]) for a in ['O1','O2']),occupancy=wa.occupancy,B_factor=wa.bfactor))
 models.append((ch,lig,heme))
save('ligand_geometry.csv',geometry);save('key_region_distances.csv',keyrows);save('reaction_head_contacts_4p5A.csv',contacts);save('arginine_oxygen_distances.csv',argrows)
if waters:save('His_Gln_nearby_waters_4A.csv',waters)
if keypolar:save('His_Gln_heme_polar_contacts_4p5A.csv',keypolar)
# Heme-aligned frame distinguishes heme placement from global fold alignment.
x=np.array([models[0][2][a].coord for a in ['FE','NA','NB','NC','ND']],float);y=np.array([models[1][2][a].coord for a in ['FE','NA','NB','NC','ND']],float)
u,_,vt=np.linalg.svd((x-x.mean(0)).T@(y-y.mean(0)));rot=u@np.diag([1,1,np.linalg.det(u@vt)])@vt;t=y.mean(0)-x.mean(0)@rot
comparisons=[]
for a,b in [(89,85),(244,242),(173,170),(174,171),(175,172),(176,173),(303,292),(83,79)]:
 r1=models[0][0][a];r2=models[1][0][b]
 l1=PDBParser(QUIET=True).get_structure('a',CASE/'EXP_OleT_KM_local_frame.pdb')[0]['A'];l2=PDBParser(QUIET=True).get_structure('b',CASE/'EXP_P450BSbeta_local_frame.pdb')[0]['A']
 comparisons.append(dict(OleT_position=a,BSbeta_position=b,OleT_residue=r1.resname,BSbeta_residue=r2.resname,CA_global_A=dist(r1['CA'],r2['CA']),CA_heme_frame_A=float(np.linalg.norm(r1['CA'].coord@rot+t-r2['CA'].coord)),CA_local_fit_A=dist(l1[a]['CA'],l2[b]['CA'])))
save('corresponding_residue_displacements.csv',comparisons)
ligdis=[]
for a in ['C1','C2','C3']:
 ligdis.append(dict(atom=a,global_displacement_A=dist(models[0][1][a],models[1][1][a]),heme_frame_displacement_A=float(np.linalg.norm(models[0][1][a].coord@rot+t-models[1][1][a].coord))))
save('ligand_head_displacements.csv',ligdis)
summary={'geometry':geometry,'nearby_water_candidates':waters,'His_Gln_heme_polar_candidates':keypolar,'heme_N4_Fe_fit_RMSD_A':float(np.sqrt(np.mean(np.sum((x@rot+t-y)**2,axis=1)))),'residue_displacements':comparisons,'ligand_head_displacements':ligdis,'ligand_atom_identity_checked_from_CCD_bonds':True,'hydrogens_or_reactive_oxo_added':False}
(OUT/'summary.json').write_text(json.dumps(summary,indent=2))


# Deposited water-network candidates: heavy-atom geometry only.
water_geometry=[]
for name,water_id,position,oxygen in [('OleT_KM',614,89,'O1'),('P450BSbeta',663,85,'O2')]:
 c=CFG[name];chain=MMCIFParser(QUIET=True).get_structure(name,STRUCTURES/f"case/native/{c['PDB']}.cif")[0]['A']
 wat=next(r for r in chain if r.id[1]==water_id and r.resname=='HOH')['O']
 lig=next(r for r in chain if r.id[1]==c['ligand_residue'] and r.resname==c['ligand'])
 atom=chain[position]['NE2']
 water_geometry.append(dict(pdb=c['PDB'],residue=position,water=water_id,ligand_O=oxygen,NE2_water_O_angle_deg=angle(atom.coord-wat.coord,lig[oxygen].coord-wat.coord),water_B=wat.bfactor,NE2_B=atom.bfactor,ligand_O_B=lig[oxygen].bfactor,water_occ=wat.occupancy))
(OUT/'water_geometry.json').write_text(json.dumps(water_geometry,indent=2))
expected=json.loads((DATA/'case/results/summary.json').read_text())
for a,b in zip(summary['geometry'],expected['geometry']):
 for key in a:
  if isinstance(a[key],float):assert abs(a[key]-b[key])<1e-5,(key,a[key],b[key])
print('Experimental case geometry reproduced; no hydrogens or reactive oxo were added.')
