# Protein Structures

Analysis-ready coordinates in PDB and mmCIF formats.

| Directory | Number of coordinate files | Contents |
|---|---:|---|
| `lal/global/reference/` | 29 | Lal reference structures |
| `lal/global/candidate/` | 74 | Lal candidate structures |
| `omega_ta/global/reference/` | 58 | Reference structures: 49 experimental records and 9 predicted models |
| `omega_ta/global/candidate/` | 52 | Candidate models |
| `omega_ta/local_pm3/reference/` | 58 | Reference active-site fragments |
| `omega_ta/local_pm3/candidate/` | 52 | Candidate active-site fragments |
| `cyp/global/` | 86 mmCIF | CYP P450-domain reference structures |
| `cyp/local_pm3/` | 86 mmCIF | CYP template-anchor ±3 fragments |
| `cyp/template/` | 1 mmCIF + 1 PDB | P450cam 2ZWU complex and resolved template CA trace |
| `cyp/case/` | 2 mmCIF + 8 PDB | Native experimental complexes, global/local superpositions, and TM-align CA inputs |

`structure_index.csv` lists the dataset, protein identifier, reference/candidate/template/case role, structural representation, residue count, and file path.

Local fragments comprise mapped template anchors and three adjacent resolved residues on each side (±3), with overlapping intervals merged. Coordinates and residue numbering follow the manuscript analysis.

CYP sequence numbering, selected domains and author-numbered experimental case positions are documented in the [CYP dataset README](../data/cyp_reference/README.md). The case complexes retain their deposited ligands and water identities.

## CYP coordinate organization

The CYP addition contains 184 coordinate files: 172 reference files (86 global domains and 86 local fragments), two template files, and ten experimental case files. These are coordinate representations of 86 reference proteins plus the template and two experimental case proteins; the file count is not a count of distinct enzymes.

| Path within `cyp/` | Interpretation |
|---|---|
| `global/<accession>.cif` | AFDB-derived P450-domain coordinates, with original residue numbering |
| `local_pm3/<accession>.cif` | Selected mapped-anchor ±3 residues, retaining all atoms |
| `template/2ZWU.cif` | Deposited P450cam–CAM complex used to identify the 14 contact anchors |
| `template/P00183_experimental_CA.pdb` | Resolved template CA trace used for alignment |
| `case/native/8VWK.cif`, `case/native/1IZO.cif` | Deposited experimental complexes with ligands, heme and waters |
| `case/EXP_*_global_frame.pdb` | Case coordinates in the global superposition, used for D/E and pocket comparisons |
| `case/EXP_*_local_frame.pdb` | Case coordinates in the independently fitted local superposition |
| `case/tm_inputs/EXP_*_global.pdb`, `case/tm_inputs/EXP_*_local.pdb` | CA-only inputs defining the exact global/local TM-score calculations |
| `structure_index.csv` | CYP-only coordinate inventory; also included in the main structure index |

Local-frame case files retain the complete selected global domain and its experimental ligands; the local fit uses only the 65-residue selection specified in `data/cyp_reference/case/metadata.json`. Panel E establishes residue correspondence using the global alignment, then displays the union of the two local selections (66 columns). Do not interpret the local-frame file's total residue count as the local-window length.

The case labels distinguish the documented common substrate from the crystallographic ligands: 8VWK contains MYR, whereas 1IZO contains PAM. Water-network calculations use the deposited water atoms and do not introduce hydrogens or a modeled reactive oxo. See [analysis commands](../scripts/README.md#running-the-cyp-analysis).
