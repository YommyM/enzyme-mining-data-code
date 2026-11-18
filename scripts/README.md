# Analysis Scripts

| Script | Output |
|---|---|
| `plot_lal_uniprot_branch.py` | Figure 2 |
| `plot_omega_local_advantage.py` | Figure 3 and per-query functional benchmark |
| `plot_experimental_validation.py` | Figure 4 |
| `plot_reference_spaces.py` | Figure S1 |
| `plot_omega_window_sensitivity.py` | Figure S2 |
| `plot_experimental_evidence_si.py` | Figures S3–S4 |
| `replot_ta_vector_chromatograms.py` | Figures S5–S7 |
| `reproduce_benchmark.py` | Recalculated benchmark scores and statistical tests |
| `extract_cyp_local.py` | Validate 14 contact anchors and reproduce the 86 ±3 local fragments |
| `compute_cyp_tm.py` | Optional global/local TM-align calculations and experimental case scores |
| `analyze_cyp_function.py` | Reproduce Figure 5A–C clustering, AUCs, bootstrap intervals and paired tests |
| `plot_cyp_function.py` | Individual Figure 5A/B/C panels |
| `analyze_cyp_case.py` | Experimental case ligand/heme geometry, pocket displacements and water contacts |
| `plot_cyp_pocket_alignment.py` | Individual Figure 5E sequence/structure panel |
| `cyp_common.py` | Package-relative I/O shared by CYP scripts |

Scripts use package-relative paths and write outputs to `figures/`. The Figure 3 analysis uses Euclidean distances between complete structural similarity profiles, average-linkage clustering, overlapping substrate-group memberships, group-balanced scores, bootstrap confidence intervals, and exact paired sign-flip tests.

The benchmark script verifies per-query results and the six active-site-window comparisons against the supplied result tables. Installation and execution commands are provided in the root README.

CYP scripts use only the Figure 5 cohort and final methods. They do not generate narrative reports, manuscript page assemblies, or ChimeraX views. Case analysis outputs support the author-created D/F views. See the [CYP dataset README](../data/cyp_reference/README.md) for data definitions and statistical methods.

## Running the CYP analysis

From the package root, install NumPy, SciPy, Matplotlib, Biopython and scikit-learn, then run:

```bash
python scripts/extract_cyp_local.py
python scripts/analyze_cyp_function.py
python scripts/plot_cyp_function.py --results figures/cyp/analysis
python scripts/analyze_cyp_case.py
python scripts/plot_cyp_pocket_alignment.py
```

| Script | Main inputs | Output under `figures/cyp/` |
|---|---|---|
| `extract_cyp_local.py` | 86 global domains, 2ZWU CAM contacts, archived template-first alignments | `extracted_local/`: fragments and extraction checks |
| `analyze_cyp_function.py` | 86 × 86 TM matrices, reaction records, domain alignment, query metadata | `analysis/`: A/B/C statistics, bootstrap intervals and Holm-adjusted p values |
| `plot_cyp_function.py` | Final result tables, or the directory passed with `--results` | Individual `Figure_5A`, `Figure_5B`, `Figure_5C` PNG/PDF/SVG files |
| `analyze_cyp_case.py` | Native and aligned experimental coordinates, case metadata | `case/geometry/`: ligand/heme, residue-displacement and water-contact tables |
| `plot_cyp_pocket_alignment.py` | Experimental global alignment, selected residues, deposited secondary structure | `case/Figure_5E` PNG/PDF/SVG and residue-level analysis tables |

The functional script reconstructs the 73 sequence clusters and the 58-query/seven-class selection. It checks all plotted ARIs and query AUCs, C bootstrap intervals and corrected p values against the supplied results. C uses the matched 56-query/six-class set. The extraction script checks all 86 fragments against the supplied coordinates, and the case scripts check the geometry and panel E correspondence.

Optional coordinate-level TM-score reproduction requires an external TM-align executable (study version 20220412):

```bash
python scripts/compute_cyp_tm.py --scope case --tmalign /path/to/TMalign
python scripts/compute_cyp_tm.py --scope all --workers 4 --tmalign /path/to/TMalign
```

Both input orders and both length normalizations are averaged. The full run evaluates 3,655 unordered pairs per representation. Its logs, rotation matrices and scores are written to `figures/cyp/tm_align/`. Existing Figure 5 panels can be regenerated from the supplied matrices without TM-align.
