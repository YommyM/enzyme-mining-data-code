#!/usr/bin/env python3
"""Reproduce Figure 5 local fragments from template-first alignments and ±3 windows."""
import argparse
import numpy as np
from Bio.PDB import MMCIFParser, MMCIFIO, Select
from Bio.PDB.Polypeptide import protein_letters_3to1_extended as AA
from cyp_common import DATA, STRUCTURES, OUTPUT, read_csv, write_csv

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=str, default=str(OUTPUT / 'extracted_local'))
    args = parser.parse_args()
    from pathlib import Path
    output = Path(args.output); output.mkdir(parents=True, exist_ok=True)
    template = read_csv(DATA / 'template_residues.csv')
    seeds = {int(r['observed_index']) for r in template if r['seed_6A'] == 'True'}
    raw = MMCIFParser(QUIET=True).get_structure('template', STRUCTURES / 'template/2ZWU.cif')[0]['A']
    ligand = next(r for r in raw if r.resname == 'CAM')
    ligand_xyz = np.array([a.coord for a in ligand if a.element not in ['H', 'D']])
    # Heavy-atom contact anchors are rechecked against the deposited CAM complex.
    for row in template:
        residue = raw[int(row['author_position'])]
        xyz = np.array([a.coord for a in residue if a.element not in ['H', 'D']])
        nearest = np.linalg.norm(xyz[:, None] - ligand_xyz, axis=2).min()
        assert (nearest <= 6) == (int(row['observed_index']) in seeds)
    assert len(seeds) == 14
    inventory = []; mappings = []
    for meta in read_csv(DATA / 'protein_metadata_86.csv'):
        accession = meta['accession']
        model = MMCIFParser(QUIET=True).get_structure(accession, STRUCTURES / 'global' / f'{accession}.cif')
        residues = list(model[0].get_residues())
        lines = (DATA / 'template_alignments' / f'{accession}.txt').read_text().splitlines()
        start = next(i for i, line in enumerate(lines) if 'denotes residue pairs' in line)
        template_seq, target_seq = lines[start + 1].strip(), lines[start + 3].strip()
        assert template_seq.replace('-', '') == ''.join(r['amino_acid'] for r in template)
        assert target_seq.replace('-', '') == ''.join(AA[r.resname] for r in residues)
        i = j = -1; selected = set()
        for a, b in zip(template_seq, target_seq):
            i += a != '-'; j += b != '-'
            if a == '-' or i not in seeds: continue
            window = list(range(max(0, j - 3), min(len(residues), j + 4))) if j >= 0 else []
            selected.update(window)
            mappings.append(dict(accession=accession, template_observed_index=i,
                target_resolved_index=j, mapping_status='direct' if b != '-' else ('preceding_target' if j >= 0 else 'omitted'),
                window_resolved_indices=';'.join(map(str, window))))
        chosen = [residues[i] for i in sorted(selected)]
        packaged = list(MMCIFParser(QUIET=True).get_structure('local', STRUCTURES / 'local_pm3' / f'{accession}.cif').get_residues())
        assert [r.id for r in chosen] == [r.id for r in packaged]
        for a, b in zip(chosen, packaged):
            assert np.array_equal(np.array([x.coord for x in a]), np.array([x.coord for x in b]))
        class Selection(Select):
            def accept_residue(self, residue): return id(residue) in selected_objects
        selected_objects = {id(r) for r in chosen}
        writer = MMCIFIO(); writer.set_structure(model); writer.save(str(output / f'{accession}.cif'), Selection())
        inventory.append(dict(accession=accession, residues=len(chosen), packaged_coordinates_identical=True))
    write_csv(output / 'inventory.csv', inventory); write_csv(output / 'anchor_mapping.csv', mappings)
    print(f'Validated and extracted {len(inventory)} local structures; 14 contact anchors per target.')

if __name__ == '__main__': main()
