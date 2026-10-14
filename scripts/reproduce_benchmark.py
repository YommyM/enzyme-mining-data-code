#!/usr/bin/env python3
"""Verify manuscript and window statistics using the packaged matrices only."""

from pathlib import Path
import csv
import numpy as np
import plot_omega_local_advantage as methods

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'data/omega_reference'
RANGE = ROOT / 'data/omega_pocket_range_current'


def read(path):
    return list(csv.DictReader(path.open()))


def matrix(path, ids):
    rows = list(csv.reader(path.open()))
    assert rows[0][1:] == ids and [r[0] for r in rows[1:]] == ids
    values = np.array([[float(x) for x in r[1:]] for r in rows[1:]])
    assert np.allclose(values, values.T) and np.allclose(np.diag(values), 1)
    return values


def main():
    rows = read(DATA / 'reference_function_reclassification.csv')
    ids = [r['structure_id'] for r in rows]
    levels = {r['structure_id'].lower(): r for r in read(DATA / 'reference_function_levels.csv')}
    memberships = methods.functional_memberships(rows, levels)
    primary = {i.lower(): levels[i.lower()]['substrate_class'] for i in ids}
    matrices = {'global': DATA / 'global_tm_original_workflow_58.csv'}
    matrices.update({f'window_{w}': RANGE / f'local_tm_original_mapping_window_{w}_58.csv'
                     for w in [3,5,7,10,15,20]})
    records = {name: methods.profile_benchmark(matrix(path, ids), ids, rows, memberships, primary)
               for name, path in matrices.items()}
    metrics = [f'recovery_{k}' for k in range(1,6)] + ['average_precision']
    expected = {(r['representation'], r['structure_id']): r
                for r in read(DATA / 'original_workflow_multilabel_profile_benchmark.csv')}
    for name in ['global', 'window_3']:
        original = 'global' if name == 'global' else 'local_window_3'
        for row in records[name]:
            saved = expected[(original, row['structure_id'])]
            for metric in metrics:
                assert abs(row[metric]-float(saved[metric])) < 1e-12, (name, row['structure_id'], metric)
    expected_windows = {r['representation']: r for r in read(RANGE / 'pocket_range_summary.csv')}
    summary=[]
    for name, queries in records.items():
        for metric in metrics:
            scores = methods.family_summary(queries, metric)
            baseline = methods.family_summary(records['global'], metric)
            families = list(scores)
            values = [scores[f] for f in families]
            mean=float(np.mean(values))
            p=methods.exact_sign_flip([scores[f]-baseline[f] for f in families])
            low, high=methods.bootstrap_ci(values)
            if name != 'global':
                saved = expected_windows[name]
                assert abs(mean-float(saved[metric])) < 1e-12, (name,metric)
                assert abs(p-float(saved['p_'+metric])) < 1e-12, (name,metric,'p')
            summary.append(dict(representation=name,metric=metric,mean=mean,
                                ci_low=low,ci_high=high,p_vs_global=p))
    out=ROOT/'figures/benchmark_check';out.mkdir(parents=True,exist_ok=True)
    with (out/'reproduced_summary.csv').open('w',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(summary[0]));writer.writeheader();writer.writerows(summary)
    print('PASS: all primary per-query metrics and six-window mean scores/P values match saved results.')
    for row in summary:
        if row['metric']=='average_precision':
            print(row['representation'], 'AP=',round(row['mean'],6),'P=',round(row['p_vs_global'],6))


if __name__=='__main__':
    main()
