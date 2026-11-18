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

## CYP references and experimental case

Directory: [`cyp_reference/`](cyp_reference/README.md). These are the inputs and final results for Figure 5.

| File or directory | Contents |
|---|---|
| `protein_metadata_86.csv` | All 86 reference proteins, organisms, domain boundaries, reaction memberships, sequence clusters and query roles |
| `reaction_annotations.csv` | 402 evidence-linked reaction records, covering 281 enzyme–substrate associations and 127 standardized substrates |
| `single_label_queries_58.csv` | The seven-class, 58-query evaluation set for A/B |
| `aligned_domains_86.fasta` | Domain alignment for reconstructing the 73 reference sequence clusters |
| `global_tm_86.csv`, `local_tm_86.csv` | Global-domain and ±3 pocket TM-score matrices; the same 86 identifiers and order |
| `template_residues.csv`, `template_alignments/` | P450cam 2ZWU contact anchors and target mappings |
| `anchor_mapping_86.csv`, `local_region_inventory_86.csv` | Selected ±3 positions, gap handling and region sizes |
| `results/` | A clustering scores, B class-level AUCs, B/C query-level AUCs, C bootstrap intervals and corrected significance tests |
| `case/metadata.json` | Experimental protein, ligand/heme and selected-residue identities |
| `case/experimental_residue_number_mapping.csv`, `case/alignment/` | UniProt/author numbering correspondence and experimental TM-align records |
| `case/results/`, `case/TM_scores.json` | Panel E sequence/structure data and D/F geometric analyses |

All 86 proteins contribute reference dimensions and retrieval candidates. A/B evaluate 58 single-label queries in seven classes; C uses the same 56 eligible queries in six classes across all three conditions. Multi-label proteins remain in the reference pool. The 14 coarse annotation classes and seven evaluation classes are connected by explicit net-reaction refinement and the single-label/minimum-two-member selection, as detailed in the [dataset README](cyp_reference/README.md).

The case consists of OleT_KM (A0A0B0D9P4, 8VWK chain A) and P450BSβ (O31440, 1IZO chain A). Each has 65 selected local residues. Their experimental global/local TM-scores are 0.913725/0.87375. Protein annotations describe their myristate activities; the crystal ligands are MYR (myristic acid) and PAM (palmitoleic acid). Original and aligned coordinates are in `../structures/cyp/`.

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
