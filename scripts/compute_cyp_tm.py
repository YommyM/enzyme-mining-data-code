#!/usr/bin/env python3
"""Recompute symmetric TM-scores from the Figure 5 coordinates using external TM-align."""
import argparse,csv,itertools,json,re,subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np
from Bio.PDB import MMCIFParser
from Bio.PDB.Polypeptide import protein_letters_3to1_extended
from Bio.Data.IUPACData import protein_letters_1to3
from cyp_common import DATA,STRUCTURES,OUTPUT,load_inputs,write_csv

def ca_trace(source,destination):
    residues=list(MMCIFParser(QUIET=True).get_structure('model',source).get_residues());lines=[]
    for i,r in enumerate(residues,1):
        name=protein_letters_1to3.get(protein_letters_3to1_extended.get(r.resname,'X'),'UNK').upper();x,y,z=r['CA'].coord
        lines.append(f'ATOM  {i:5d}  CA  {name:>3s} A{i:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           C\n')
    destination.parent.mkdir(parents=True,exist_ok=True);destination.write_text(''.join(lines)+'TER\nEND\n')

def pair_scores(a,b,paths,folder,executable):
    scores=[]
    for u,v in [(a,b),(b,a)]:
        stem=folder/f'{u}__{v}';output=subprocess.run([executable,str(paths[u]),str(paths[v]),'-m',str(stem.with_suffix('.matrix'))],capture_output=True,text=True,check=True).stdout
        stem.with_suffix('.log').write_text(output)
        values=[float(x) for x in re.findall(r'^TM-score=\s*([0-9.]+)',output,re.M)]
        if len(values)!=2:raise ValueError('Expected both length-normalized TM-scores')
        scores.extend(values)
    return dict(accession_a=a,accession_b=b,TM_mean=float(np.mean(scores)),forward_length_a=scores[0],forward_length_b=scores[1],reverse_length_b=scores[2],reverse_length_a=scores[3])

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--tmalign',default='TMalign');p.add_argument('--scope',choices=['case','global','local','all'],default='case');p.add_argument('--workers',type=int,default=1);p.add_argument('--accessions',nargs='+',help='Optional reference subset for verification')
    args=p.parse_args();out=OUTPUT/'tm_align';out.mkdir(parents=True,exist_ok=True)
    if args.scope=='case':
        expected=json.loads((DATA/'case/TM_scores.json').read_text());results=[]
        for region in ['global','local']:
            folder=out/f'case_{region}';folder.mkdir(exist_ok=True)
            names=['OleT_KM','P450BSbeta'];paths={n:STRUCTURES/f'case/tm_inputs/EXP_{n}_{region}.pdb' for n in names}
            r=pair_scores(*names,paths,folder,args.tmalign);r['representation']=region
            r['difference_from_packaged']=r['TM_mean']-expected['EXP_'+region]['TM_mean'];results.append(r)
        write_csv(out/'case_scores.csv',results);print(results);return
    _,ids,matrices=load_inputs();selected=args.accessions or ids
    if not set(selected)<=set(ids):raise ValueError('Unknown accession')
    for rep in (['global','local'] if args.scope=='all' else [args.scope]):
        folder=out/rep;folder.mkdir(exist_ok=True);paths={}
        for a in selected:
            dest=folder/f'{a}.pdb';ca_trace(STRUCTURES/('global' if rep=='global' else 'local_pm3')/f'{a}.cif',dest);paths[a]=dest
        pairs=list(itertools.combinations(selected,2))
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            rows=list(pool.map(lambda ab:pair_scores(*ab,paths,folder,args.tmalign),pairs))
        matrix=np.eye(len(selected));reference=matrices[0 if rep=='global' else 1]
        for r in rows:
            a,b=r['accession_a'],r['accession_b'];i,j=selected.index(a),selected.index(b);matrix[i,j]=matrix[j,i]=r['TM_mean'];r['difference_from_packaged']=r['TM_mean']-reference[ids.index(a),ids.index(b)]
        write_csv(folder/'pair_scores.csv',rows)
        with (folder/'matrix.csv').open('w',newline='') as f:
            w=csv.writer(f);w.writerow(['accession']+selected);w.writerows([[a]+list(matrix[i]) for i,a in enumerate(selected)])
        print(rep,len(selected),'proteins; maximum absolute score difference:',max((abs(r['difference_from_packaged']) for r in rows),default=0))

if __name__=='__main__':main()
