# BEiGUIDE

[![Tests](https://github.com/agmcfarland/BEiGUIDE/actions/workflows/test-build.yml/badge.svg)](https://github.com/agmcfarland/BEiGUIDE/actions/workflows/test-build.yml)
[![codecov](https://codecov.io/gh/agmcfarland/BEiGUIDE/graph/badge.svg?token=NPALNGNUFJ)](https://codecov.io/gh/agmcfarland/BEiGUIDE)

# Background

**Base Editor iGUIDE (BEiGUIDE)** quantifies base editing at both on- and off-target cut sites detected by iGUIDE when using CRISPR base editor-Cas9 fusion enzymes.

Following an iGUIDE analysis, BEiGUIDE identifies candidate edit sites, quantifies nucleotide composition surrounding each cut site, and performs statistical testing to determine whether observed base conversions significantly exceed the expected sequencing error rate.

For additional details, please see our preprint:

https://www.biorxiv.org/content/10.1101/2025.08.26.667396v1.abstract

---

# Installation

```r
devtools::install_github("agmcfarland/BEiGUIDE")
```

---

# Workflow

A typical BEiGUIDE analysis consists of two steps.

1. Quantify edits from iGUIDE output.

```r
BEiGUIDE::quantify_edits(
  base_directory = "path/to/base_directory",
  output_directory = "path/to/output_directory",
  analysis_name = "my_analysis",
  abundance_cutoff = 5,
  expected_PAM_proximal_cut_distance = 3,
  expected_PAM_distal_cut_distance = 17,
  allowed_aln_start_cut_site_PAM_distal = 4,
  allowed_aln_start_cut_site_PAM_proximal = 4,
  reference_genome_path = "path/to/genome.fasta",
  editable_base = "A",
  expected_edit = "G",
  null_proportion = 0.001,
  binomial_p_value_threshold = 0.05,
  binomial_direction = "greater",
  p_value_adj_method = "BH",
  n_processors = 6,
  n_edit_site_processors = 3,
  overwrite = FALSE
)
```

2. Generate diagnostic plots and run-level summaries.

```r
BEiGUIDE::quantify_edits_analysis(
  quantify_edits_output_path = "path/to/quantify_edits"
)
```

---

# `quantify_edits()` Parameters

| Parameter | Type | Description | Default |
|------------|------|-------------|---------|
| `base_directory` | Character | Base directory containing iGUIDE output | **Required** |
| `output_directory` | Character | Directory where results will be written | **Required** |
| `analysis_name` | Character | Name of analysis output directory | `"quantify_edits"` |
| `abundance_cutoff` | Numeric | Minimum read count required for an edit site | `3` |
| `expected_PAM_proximal_cut_distance` | Numeric | Expected PAM-proximal cut distance | `3` |
| `expected_PAM_distal_cut_distance` | Numeric | Maximum PAM-distal distance examined | `17` |
| `allowed_aln_start_cut_site_PAM_distal` | Numeric | Allowed alignment start within the protospacer | `3` |
| `allowed_aln_start_cut_site_PAM_proximal` | Numeric | Allowed alignment start outside the protospacer | `3` |
| `reference_genome_path` | Character | FASTA reference genome used to retrieve reference sequence | `""` |
| `editable_base` | Character | Base expected to undergo editing | `"A"` |
| `expected_edit` | Character | Expected edited base | `"G"` |
| `null_proportion` | Numeric | Expected background sequencing error rate | `0.001` |
| `binomial_p_value_threshold` | Numeric | Significance threshold for binomial testing | `0.05` |
| `binomial_direction` | Character | Alternative hypothesis (`greater`, `less`, `two.sided`) | `"greater"` |
| `p_value_adj_method` | Character | Multiple-testing correction method passed to `p.adjust()` | `"BH"` |
| `n_processors` | Numeric | Number of processors for BAM parsing | `4` |
| `n_edit_site_processors` | Numeric | Number of processors for per-edit-site analysis | `3` |
| `overwrite` | Logical | Overwrite existing output directory | `TRUE` |

---

# `quantify_edits_analysis()`

After `quantify_edits()` has completed, this function generates publication-ready summaries.

For each edit site passing the abundance filter it:

- generates nucleotide heatmaps
- generates on-target-only plots
- generates sequencing depth plots
- generates stacked percentage plots

It also produces:

**Significant base edits found and what they are**

- `significant_base_edits.csv`
- `significant_base_edits.rds`

**Summary of total cuts and base edits observed**

- `run_base_edit_summary.csv`
- `run_base_edit_summary.rds`

**Cut sites that passed abundance thresholds for base edit analysis but didn't have alignments starting/ending within the allowed cut window**

- `edit_site_fail.csv`
- `edit_site_fail.rds`

All plots are written to the `plots/` directory inside the analysis output folder.

---

# Output

Running `quantify_edits()` produces:

- `run_parameters.csv` / `.rds`
- `edit_sites_overview.csv` / `.rds`
- one directory for each detected edit site containing:
  - processed alignments
  - base composition tables
  - statistical results
  - intermediate data
- analysis log file

Running `quantify_edits_analysis()` additionally produces:

- PDF diagnostic plots
- combined significant base edit tables
- run-level summary tables

---

# Notes

`n_processors` controls how many BAM files are processed simultaneously. Increasing this value increases RAM usage substantially. If memory becomes limiting, reduce this parameter.

`n_edit_site_processors` controls parallel processing of edit sites after they have been identified. This stage is typically much less memory intensive than BAM processing.

---

# Future additions

- Command-line wrapper
- Docker image