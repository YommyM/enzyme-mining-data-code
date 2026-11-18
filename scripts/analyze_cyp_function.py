#!/usr/bin/env python3
"""Reproduce Figure 5A–C: profile clustering, class-wise AUC and cluster-paired inference."""
import json
from collections import Counter, defaultdict
import numpy as np
from scipy.cluster.hierarchy import linkage, cut_tree
from scipy.spatial.distance import pdist, squareform
from scipy.sparse.csgraph import connected_components
from Bio import SeqIO
from sklearn.metrics import adjusted_rand_score, roc_auc_score
from cyp_common import DATA, OUTPUT, read_csv, write_csv, load_inputs

CONDITIONS = ['all_references', 'exclude_cluster', 'exclude_cluster_features']
CLASSES = ['hydroxylation', 'coupling_cyclization', 'oxidative_dealkylation',
           'oxidative_decarboxylation', 'overall:methyl_to_carboxylate', 'nitration', 'rearrangement']

def bootstrap(g, l, clusters, taxa):
    groups = defaultdict(list)
    for c in sorted(set(clusters)):
        ix = np.flatnonzero(clusters == c)
        assert len(set(taxa[ix])) == 1
        groups[taxa[ix][0]].append(ix)
    rng = np.random.default_rng(2026111702); samples = []; attempts = 0
    while len(samples) < 5000:
        attempts += 1
        if attempts > 50000: raise RuntimeError('Too few evaluable bootstrap draws')
        draw = np.concatenate([np.concatenate([u[i] for i in rng.integers(0, len(u), len(u))]) for u in groups.values()])
        gg, ll = g[draw], l[draw]
        if not (np.isfinite(gg).any(0).all() and np.isfinite(ll).any(0).all()): continue
        gm, lm = np.nanmean(gg, axis=0).mean(), np.nanmean(ll, axis=0).mean()
        samples.append([gm, lm, lm-gm])
    a = np.array(samples)
    return {**{name+'_ci95': np.quantile(a[:,j], [.025,.975]).tolist() for j,name in enumerate(['global','local','delta'])},
            'valid_bootstraps': len(samples), 'attempts': attempts}

def significance(query_rows, eligible, cluster_by_id):
    units = sorted({cluster_by_id[a] for a in eligible}); contrib = []; results = []
    for condition in CONDITIONS:
        rr = sorted([r for r in query_rows if r['condition']==condition and r['accession'] in eligible], key=lambda r:r['accession'])
        counts = Counter(r['category'] for r in rr)
        w = np.array([1/len(counts)/counts[r['category']] for r in rr])
        g,l = [np.array([r[rep+'_AUC'] for r in rr]) for rep in ['global','local']]
        v = np.array([sum(w[j]*(l[j]-g[j]) for j,r in enumerate(rr) if cluster_by_id[r['accession']]==c) for c in units])
        contrib.append(v)
        results.append(dict(condition=condition,queries=len(rr),classes=len(counts),sequence_clusters=len(units),global_macro_AUC=float(w@g),local_macro_AUC=float(w@l),delta=float(w@(l-g))))
    contrib = np.array(contrib).T; rng=np.random.default_rng(2026111802); N=200000; hits=np.zeros(3,int)
    observed=np.abs(contrib.sum(0))
    for _ in range(N//10000):
        null=rng.choice([-1,1],size=(10000,len(units)))@contrib
        hits+=(np.abs(null)>=observed-1e-12).sum(0)
    p=(hits+1)/(N+1);order=np.argsort(p);holm=np.empty(3)
    holm[order]=np.minimum(1,np.maximum.accumulate(p[order]*np.arange(3,0,-1)))
    for j,r in enumerate(results):r.update(p_two_sided=float(p[j]),p_Holm=float(holm[j]),extreme_draws=int(hits[j]))
    return dict(test='Two-sided paired sequence-cluster sign-flip of class-balanced AUC; fixed reference set',draws=N,seed=2026111802,multiple_testing='Holm across three conditions',results=results)

def main():
    out=OUTPUT/'analysis';out.mkdir(parents=True,exist_ok=True)
    meta,ids,mats=load_inputs();clusters=np.array([int(r['sequence_cluster']) for r in meta]);taxa=np.array([r['taxonomic_group'] for r in meta])
    # Reconstruct the exclusion partition from the supplied domain alignment.
    aligned={r.id:str(r.seq).upper() for r in SeqIO.parse(DATA/'aligned_domains_86.fasta','fasta')}
    adjacency=np.eye(len(ids),dtype=bool)
    for i in range(len(ids)):
        for j in range(i):
            a,b=aligned[ids[i]],aligned[ids[j]];den=min(len(a.replace('-','')),len(b.replace('-','')))
            identity=sum(x==y and x!='-' for x,y in zip(a,b))/den
            coverage=sum(x!='-' and y!='-' for x,y in zip(a,b))/den
            adjacency[i,j]=adjacency[j,i]=identity>=.7 and coverage>=.8
    n,derived=connected_components(adjacency,directed=False)
    assert n==73 and np.array_equal(derived[:,None]==derived,clusters[:,None]==clusters)
    # Recreate memberships and query selection directly from deposited reaction annotations.
    memberships=defaultdict(set)
    for r in read_csv(DATA/'reaction_annotations.csv'):
        memberships[r['accession']].add(r['reaction_type'] if r['reaction_type']!='multistep_composite' else 'overall:'+r['reaction_subtype'])
    for r in meta:assert memberships[r['accession']]==set(r['reaction_labels'].split(';'))
    single={a:next(iter(v)) for a,v in memberships.items() if len(v)==1}
    counts=Counter(single.values());selected={a:c for a,c in single.items() if counts[c]>=2}
    cohort=read_csv(DATA/'single_label_queries_58.csv');assert selected=={r['accession']:r['label'] for r in cohort}
    sx=[ids.index(r['accession']) for r in cohort];truth=[r['label'] for r in cohort]
    full=[squareform(pdist(m)) for m in mats];curves=[]
    for rep,d in zip(['global','local'],full):
        tree=linkage(squareform(d[np.ix_(sx,sx)],checks=False),method='average')
        for k in range(2,11):curves.append(dict(distance='profile',representation=rep,clusters=k,proteins=58,ARI=float(adjusted_rand_score(truth,cut_tree(tree,n_clusters=k).ravel()))))
    query_rows=[];methods={c:[np.full((86,7),np.nan) for _ in range(2)] for c in CONDITIONS}
    for row in cohort:
        i=ids.index(row['accession']);j=CLASSES.index(row['label']);targets=np.array([row['label'] in memberships[a] for a in ids]);external=clusters!=clusters[i]
        for condition,allowed in zip(CONDITIONS,[np.arange(86)!=i,external,external]):
            y=targets[allowed]
            if len(np.unique(y))<2:continue
            for k,m in enumerate(mats):
                d=np.linalg.norm(m[allowed][:,external]-m[i,external],axis=1) if condition=='exclude_cluster_features' else full[k][i,allowed]
                methods[condition][k][i,j]=roc_auc_score(y,-d)
            query_rows.append(dict(condition=condition,accession=ids[i],category=row['label'],same_class_candidates=int(y.sum()),other_candidates=int(len(y)-y.sum()),global_AUC=methods[condition][0][i,j],local_AUC=methods[condition][1][i,j]))
    perclass=[]
    for j,c in enumerate(CLASSES):
        for condition in [CONDITIONS[0]]:
            values=[np.nanmean(a[:,j]) for a in methods[condition]]
            perclass.append(dict(category=c,queries=counts[c],condition=condition,global_AUC=values[0],local_AUC=values[1]))
    valid=np.isfinite(methods[CONDITIONS[-1]][0]);active=valid.any(0);summaries=[]
    for condition in CONDITIONS:
        g,l=[np.where(valid,a,np.nan)[:,active] for a in methods[condition]]
        summaries.append(dict(condition=condition,matched_queries=int(valid.any(1).sum()),matched_classes=int(active.sum()),global_AUC=float(np.nanmean(g,axis=0).mean()),local_AUC=float(np.nanmean(l,axis=0).mean()),**bootstrap(g,l,clusters,taxa)))
    eligible={ids[i] for i in np.flatnonzero(valid.any(1))};tests=significance(query_rows,eligible,dict(zip(ids,clusters)))
    # Compare every plotted statistic with the frozen Figure 5 source tables.
    expected=read_csv(DATA/'results/A_ARI.csv')
    for a,b in zip(curves,expected):assert abs(a['ARI']-float(b['ARI']))<1e-12
    expected={(r['condition'],r['accession']):r for r in read_csv(DATA/'results/BC_query_AUC.csv')}
    for r in query_rows:
        for rep in ['global','local']:assert abs(r[rep+'_AUC']-float(expected[r['condition'],r['accession']][rep+'_AUC']))<1e-12
    previous=json.loads((DATA/'results/C_significance.json').read_text())
    for a,b in zip(tests['results'],previous['results']):assert abs(a['p_Holm']-b['p_Holm'])<1e-12
    previous=json.loads((DATA/'results/C_summary.json').read_text())
    for a,b in zip(summaries,previous):
        for key in ['global_ci95','local_ci95','delta_ci95']:assert np.allclose(a[key],b[key],atol=1e-12,rtol=0)
    write_csv(out/'A_ARI.csv',curves);write_csv(out/'B_class_AUC.csv',perclass);write_csv(out/'BC_query_AUC.csv',query_rows)
    (out/'C_summary.json').write_text(json.dumps(summaries,indent=2));(out/'C_significance.json').write_text(json.dumps(tests,indent=2))
    print('Figure 5A–C reproduced: 86 references; 58 A/B queries; 56 C queries. All source-result checks passed.')

if __name__=='__main__':main()
