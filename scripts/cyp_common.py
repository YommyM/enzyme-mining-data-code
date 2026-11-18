"""Shared I/O for the Figure 5 CYP analyses."""
from pathlib import Path
import csv
import numpy as np
ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'data/cyp_reference'
STRUCTURES = ROOT / 'structures/cyp'
OUTPUT = ROOT / 'figures/cyp'

def read_csv(path):
    with Path(path).open(encoding='utf-8-sig') as handle:
        return list(csv.DictReader(handle))

def write_csv(path, rows):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        return
    with path.open('w', encoding='utf-8-sig', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

def load_inputs():
    metadata = read_csv(DATA / 'protein_metadata_86.csv')
    ids = [r['accession'] for r in metadata]
    matrices = []
    for name in ['global', 'local']:
        with (DATA / f'{name}_tm_86.csv').open() as handle:
            rows = list(csv.reader(handle))
        assert rows[0][1:] == ids == [r[0] for r in rows[1:]]
        matrix = np.array([r[1:] for r in rows[1:]], dtype=float)
        assert np.allclose(matrix, matrix.T) and np.allclose(np.diag(matrix), 1)
        matrices.append(matrix)
    return metadata, ids, matrices
