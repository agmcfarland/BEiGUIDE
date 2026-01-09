# BEiGUIDE

[![Tests](https://github.com/agmcfarland/BEiGUIDE/actions/workflows/test-build.yml/badge.svg)](https://github.com/agmcfarland/BEiGUIDE/actions/workflows/test-build.yml)
[![codecov](https://codecov.io/gh/agmcfarland/BEiGUIDE/graph/badge.svg?token=NPALNGNUFJ)](https://codecov.io/gh/agmcfarland/BEiGUIDE)

# Background

`Base Editor iGUIDE (BEiGUIDE)` quantifies base edits in on and off-target cut sites detected by iGUIDE when a base editor-Cas9 fusion enzyme is used.

After running [iGUIDE](https://github.com/cnobles/iGUIDE), `BEiGUIDE` analyzes each position of the protospacer to identify statistically significant changes in the observed nucleotide frequency relative to the expected nucleotide. [Please see our pre-print for more details](https://www.biorxiv.org/content/10.1101/2025.08.26.667396v1.abstract)

# Installation

```R
devtools::install_github('agmcfarland/BEiGUIDE')
```

# Running within R

```R
BEiGUIDE::quantify_edits(
  base_directory = "path/to/base_directory",
  output_directory = "path/to/output_directory",
  analysis_name = "my_analysis",
  abundance_cutoff = 5,
  end_distance_from_cut_site = 25,
  cut_site_start_distance_within_gRNA = 4,
  cut_site_start_distance_outside_gRNA = 4,
  reference_genome_path = 'path/to/genome.fasta',
  editable_base = 'A',
  expected_edit = 'G',
  binomial_p_value_threshold = 0.05,
  binomial_direction = 'greater',
  n_processors = 6,
  overwrite = FALSE
)

BEiGUIDE::quantify_edits_analysis(
	quantify_edits_output_path = '/path/to/quantify_edits/output'
	)
```

# quantify_edits() Parameters

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| base_directory | Character | Base directory where input data is located | REQUIRED |
| output_directory | Character | Directory where output data will be saved | REQUIRED |
| analysis_name | Character | Name of the analysis | "quantify_edits" |
| abundance_cutoff | Numeric | Minimum read count for including an edit site | 3 |
| expected_cut_distance_from_pam | Numeric | Distance from PAM where the cut site is expected | 3 |
| end_distance_from_cut_site | Numeric | Maximum distance from the cut site to consider bases | 20 |
| cut_site_start_distance_within_gRNA | Numeric | Distance within the gRNA to start considering cut sites | 3 |
| cut_site_start_distance_outside_gRNA | Numeric | Distance outside the gRNA to start considering cut sites | 3 |
| reference_genome_path | Character | Path to a FASTA file used to retrieve reference genome sequence | REQUIRED |
| editable_base | Character | Original base targeted for editing | 'A' |
| expected_edit | Character | Base that the editable base is expected to be converted into | 'G' |
| binomial_p_value_threshold | Numeric | P-value threshold for binomial test of significance | 0.05 |
| binomial_direction | Character | Direction of binomial test ('greater', 'less', 'two.sided') | 'greater' |
| n_processors | Numeric | Number of CPU cores to use for parallel processing | 4 |
| overwrite | Logical | Whether to overwrite existing analysis output | TRUE |


# Notes

`n_processesors` controls how many BAM files to load into memory when combining them prior to removing unneeded information, this can cause RAM usage to explode. If a run fails due to memory, lower this. 

Edit sites are processed automatically with n-1 cores. Should not use more than a Gb of memory.

# Future additions

* Wrapper to run from command line
* Docker image


