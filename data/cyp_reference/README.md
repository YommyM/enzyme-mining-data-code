# CYP reference data — Figure 5

The reference set contains 86 CYP proteins, 281 enzyme–substrate associations, 127 standardized substrates, and 402 reaction-annotation records. All 86 proteins contribute reference dimensions and retrieval candidates. Figure 5A/B evaluate 58 single-label queries in seven classes; Figure 5C evaluates the same 56 eligible queries in six classes in each condition.

## Files

| File or directory | Contents |
|---|---|
| `protein_metadata_86.csv` | Protein identifiers, organisms, source groups, sequence clusters, complete functional memberships, query roles, and modeled P450-domain boundaries |
| `aligned_domains_86.fasta` | Aligned P450-domain sequences used to reconstruct sequence-cluster memberships |
| `reaction_annotations.csv` | Substrate structures, reaction types/subtypes, products, Rhea identifiers, evidence tiers, publications and source URLs |
| `single_label_queries_58.csv` | Figure 5A/B query membership and the seven labels |
| `global_tm_86.csv`, `local_tm_86.csv` | Symmetric structural similarity matrices, with matching accession-labelled rows/columns |
| `template_residues.csv` | Resolved P450cam residues, author/UniProt numbering, minimum CAM distance, and 14 contact-anchor flags |
| `template_alignments/` | The 86 template-first TM-align sequence alignments used for anchor mapping |
| `anchor_mapping_86.csv` | Mapped anchors, gap handling and selected ±3 resolved-residue indices |
| `local_region_inventory_86.csv` | Selected positions, fragment sizes and confidence summaries |
| `results/A_ARI.csv` | Global/local ARI curves for 2–10 clusters |
| `results/B_class_AUC.csv` | Seven class-specific mean AUCs and query counts |
| `results/BC_query_AUC.csv` | Per-query AUCs and candidate counts for each condition |
| `results/C_summary.json` | Matched-query macro AUCs and paired bootstrap intervals |
| `results/C_significance.json` | Sequence-cluster paired sign-flip tests and Holm-adjusted p values |
| `case/` | Experimental structures' metadata, residue numbering, alignments and geometry results |

Matrix entries are the mean of the two length-normalized TM-scores in each of the two input orders (four reported scores). TM-align version 20220412 was used. The global coordinates are the selected P450 domains, rather than complete multi-domain proteins. Sequence clusters are the connected components of pairs passing identity ≥0.70 and coverage ≥0.80 in the supplied domain alignment. Identity is the number of identical nongap positions divided by the shorter ungapped sequence length; coverage is the number of jointly nongap positions divided by the same denominator. The analysis script reconstructs and checks all 73 clusters; stored numeric cluster IDs are retained for reproducible random resampling.

## Functional labels and cohorts

The annotations contain 14 coarse reaction types. The generic `multistep_composite` class is represented by its explicitly annotated net-reaction subtype (`overall:<subtype>`); this produces 22 labels in the complete protein-membership matrix. No component activities are inferred from a composite reaction.

Select proteins with exactly one resulting label, then retain labels represented by at least two such proteins. The seven Figure 5A/B classes are hydroxylation (33), coupling/cyclization (9), oxidative dealkylation (6), oxidative decarboxylation (3), methyl-to-carboxylate (3), nitration (2), and rearrangement (2). The other reference proteins retain their complete, potentially overlapping memberships and remain in the reference pool.

Each protein is represented by its vector of TM-scores against all 86 references. Euclidean distances between these vectors are used for average-linkage clustering (A) and ranking candidates (B/C). AUC is the probability that an annotated same-class candidate precedes an unannotated-class candidate, with half credit for a distance tie. The query itself is excluded. Averages are taken within each class and then equally across classes.

## Local structure definition

P450cam (UniProt P00183, PDB 2ZWU chain A) is the single template. Its 14 anchors have at least one protein heavy atom within 6 Å of a CAM heavy atom. Map those anchors through the archived template-first global TM-align alignment, expand by three resolved sequence neighbors on each side in the target, and merge overlapping windows. For an anchor aligned to a gap, use the preceding resolved target residue; omit it if none exists. Keep all atoms of selected residues. `extract_cyp_local.py` checks the 6 Å contacts and reproduces all 86 supplied fragments.

## Figure 5C conditions and statistics

- **All:** all reference candidates except the query itself; all 86 profile dimensions.
- **External candidates:** remove candidates in the query's sequence cluster; retain all 86 profile dimensions.
- **External profiles:** also remove profile dimensions corresponding to references in the query's sequence cluster, then recalculate distances.

The two rearrangement queries have no same-class candidate after their shared cluster is excluded. All three C conditions therefore use the same 56 queries, six classes and 54 query sequence clusters. Box summaries use equal total weight per class; the short whiskers reach the furthest observation within 0.5 IQR of the box. All query points are retained. Diamonds indicate class-balanced macro AUC.

Confidence intervals use 5,000 paired, source-stratified query-sequence-cluster bootstrap samples with seed 2026111702. Candidate rankings are held fixed. Significance uses 200,000 two-sided paired cluster sign flips with seed 2026111802: global/local labels are exchanged jointly for all queries in a sequence cluster. The statistic is the class-balanced mean AUC difference. Monte Carlo p = (extreme draws + 1)/(draws + 1); Holm correction covers the three conditions. The plotted `**` means adjusted p < 0.01. B contains no p-value annotation.

## Experimental case — Figure 5D/E/F

OleT_KM (A0A0B0D9P4, PDB 8VWK chain A) and P450BSβ (O31440, PDB 1IZO chain A) each have 65 selected local residues. Their experimental global/local TM-scores are 0.913725/0.87375. The aligned coordinates, independently fitted local coordinates, original deposited mmCIF files and CA-only TM-align inputs are under `structures/cyp/case/`.

The two enzymes have documented myristate activities (decarboxylation and hydroxylation, respectively). The crystallographic complexes contain different ligands: MYR is myristic acid, C14:0, in 8VWK; PAM is palmitoleic acid, C16:1 cis-Δ9, in 1IZO. Their atom identities are checked against the deposited chemical-component bond definitions.

Figure E uses the **experimental global** alignment to establish sequence correspondence, then takes the union of the two ±3 selections: 66 columns, 64 paired positions, and 65 residues per protein. A dash includes positions outside the selected region, not just sequence gaps. Deposited α-helices and β-strands use author residue numbering; `//` marks omitted sequence intervals. Background color is the Cα distance in the global superposition. Chemical-class symbols describe substitutions, and His89/Gln85 is highlighted in red. The CSV records each residue and the reason for each missing selected position.

Case geometry includes carboxyl/α/β-carbon placement relative to heme, Arg contacts, His89/Gln85 and neighboring residues, global/local/heme-frame displacements, and candidate water bridges. HOH A614 in 8VWK and HOH A663 in 1IZO are deposited waters. Reported heavy-atom distances and angles describe geometric contacts; no hydrogens or reactive oxo were added. The case metadata and reaction annotations retain traceable evidence for interpretation.

## Reproduction

Run from the package root:

```bash
python scripts/extract_cyp_local.py
python scripts/analyze_cyp_function.py
python scripts/plot_cyp_function.py
python scripts/analyze_cyp_case.py
python scripts/plot_cyp_pocket_alignment.py
```

The analysis and extraction scripts validate their outputs against the supplied final tables/coordinates. To plot newly recalculated functional results, use `python scripts/plot_cyp_function.py --results figures/cyp/analysis`.

TM-align is optional for re-running structural alignment; it is not needed for the supplied-matrix analysis:

```bash
python scripts/compute_cyp_tm.py --scope case --tmalign /path/to/TMalign
python scripts/compute_cyp_tm.py --scope all --workers 4 --tmalign /path/to/TMalign
```

The full calculation has 3,655 unordered reference pairs per representation, with both input orders evaluated. Outputs go to `figures/cyp/`. Scripts generate only the individual A/B/C/E panels and analysis outputs; the author-composed D/F views and final page assembly are not generated.

## Package navigation and generated panels

- [Package overview and figure–script mapping](../../README.md)
- [All data collections](../README.md)
- [Analysis scripts and execution order](../../scripts/README.md#running-the-cyp-analysis)
- [Coordinate organization](../../structures/README.md#cyp-coordinate-organization)

`plot_cyp_function.py` exports Figure 5A/B/C separately; `plot_cyp_pocket_alignment.py` exports Figure 5E. Generated files go to `figures/cyp/`, while the supplied inputs and frozen results remain in this directory. The D/F case views use the experimental coordinates and analysis tables described above. The package contains no manuscript assembly or narrative-report generation script.

`results/reproduction_checks.json` records the performed checks: all 86 local fragments, A/B/C statistical results, case geometry, and panel E correspondence were reproduced. Experimental-case TM-align scores and the two corresponding predicted-reference pairs were recalculated; the complete 3,655-pair structural matrices were supplied from the final analysis rather than recomputed during packaging. File provenance and SHA-256 checksums are listed in the package-level manifests.
