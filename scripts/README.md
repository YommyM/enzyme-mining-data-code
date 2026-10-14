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

Scripts use package-relative paths and write outputs to `figures/`. The Figure 3 analysis uses Euclidean distances between complete structural similarity profiles, average-linkage clustering, overlapping substrate-group memberships, group-balanced scores, bootstrap confidence intervals, and exact paired sign-flip tests.

The benchmark script verifies per-query results and the six active-site-window comparisons against the supplied result tables. Installation and execution commands are provided in the root README.
