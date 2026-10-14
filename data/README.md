# Data Files

## Lal references and candidates

Directory: `lals_reference/`

| File or directory | Contents |
|---|---|
| `reference_sequences_29.fasta` | Reference sequences |
| `candidate_sequences_74.fasta` | Candidate sequences |
| `protein_metadata_103.csv` | Protein identifiers, reference/candidate roles, experimental selection, sequence lengths, and FASTA headers |
| `nr_structure_set_103.csv` | Protein identities and status labels for Figure 2 |
| `reference_sequence_identity_29.csv` | Reference sequence identity matrix for Figure S1A |
| `global_tm_similarity_nr_103.csv` | Global structural similarity matrix for Figure 2B |
| `profile_hmms/` | Eight cluster-specific profiles and one pooled profile, with sequences and alignments |

## ω-Transaminase references and candidates

Directory: `omega_reference/`

| File or directory | Contents |
|---|---|
| `reference_sequences_58.fasta` | Reference sequences |
| `omega_ta_candidates_52.fasta` | Candidate sequences |
| `reference_metadata_58.csv` | Reference identifiers, deposited functional labels, and sequence-cluster memberships |
| `candidate_metadata_52.csv` | Candidate identifiers, sequence lengths, database descriptions, and organisms |
| `reference_function_reclassification.csv` | Protein entities, biochemical descriptions, substrate scope, evidence, and benchmark roles |
| `reference_function_levels.csv` | Primary substrate groups and reaction classes |
| `sequence_similarity_58.csv` | Sequence similarity matrix for Figure S1B |
| `global_tm_original_workflow_58.csv` | Global structural similarity matrix for Figure 3A |
| `local_tm_original_workflow_window_3_58.csv` | Local ±3 structural similarity matrix for Figure 3B |
| `original_workflow_multilabel_profile_benchmark.csv` | Per-query recovery and average precision for Figure 3C,D |
| `foldtype_4.fasta` | Fold IV reference identities |
| `profile_hmms/` | Six retrieval profiles with sequences and alignments |
| `profile_source_mapping.csv` | Profile identities and source-file correspondence |

Benchmark fields: `representation`, global or local comparison; `structure_id`, query identifier; `family`, primary substrate group; `recovery_1`–`recovery_5`, top-k recovery; `average_precision`, average precision over the full ranking. Overlapping biochemical memberships are defined by `functional_memberships()` in the Figure 3 script.

## Active-site window analysis

Directory: `omega_pocket_range_current/`

- `local_tm_original_mapping_window_*_58.csv`: local structural similarity matrices.
- `pocket_range_summary.csv`: group-balanced scores, confidence intervals, and P values.
- `pocket_range_per_query.csv`: per-query results.
- `fragment_sizes.csv`: residue counts for each protein and window.
- `window_3_reproduction_validation.csv`: consistency checks for the ±3 fragments.

## Experimental data

Directory: `experimental_validation/`

| File or directory | Contents |
|---|---|
| `lal_products.csv` | Lal candidate identifiers, dipeptide products, and HRMS responses |
| `omega_ta_summary.csv` | TA1–TA10 identifiers, UniProt descriptions, predicted reactions, expression, and qualitative activity |
| `omega_ta_reactions.csv` | Amino donors, keto acceptors, expected products, and detection outcomes |
| `ta_peak_annotations.csv` | Peak identities, retention times, and standard assignments |
| `extracted_plots/` | Analytical assets for Figures S3–S7 |

HRMS responses are in instrument units. Activity symbols denote qualitative detection. Displayed HPLC traces are independently normalized.

## Matrix format

The first row and column contain protein identifiers in matching order. Structural matrices contain TM-align scores, with unit diagonals and mirrored pairwise entries. Global and local matrices retain the manuscript scoring convention.
