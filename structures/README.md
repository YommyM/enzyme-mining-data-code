# Protein Structures

Analysis-ready coordinates in PDB format.

| Directory | Number of PDB files | Contents |
|---|---:|---|
| `lal/global/reference/` | 29 | Lal reference structures |
| `lal/global/candidate/` | 74 | Lal candidate structures |
| `omega_ta/global/reference/` | 58 | Reference structures: 49 experimental records and 9 predicted models |
| `omega_ta/global/candidate/` | 52 | Candidate models |
| `omega_ta/local_pm3/reference/` | 58 | Reference active-site fragments |
| `omega_ta/local_pm3/candidate/` | 52 | Candidate active-site fragments |

`structure_index.csv` lists the dataset, protein identifier, reference/candidate role, structural representation, residue count, and file path.

Local fragments comprise mapped template anchors and three adjacent resolved residues on each side (±3), with overlapping intervals merged. Coordinates and residue numbering follow the manuscript analysis.
