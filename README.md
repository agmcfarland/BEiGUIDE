# BEiGUIDE

<!-- Badges start -->
[![Tests](https://github.com/agmcfarland/BEiGUIDE/actions/workflows/test-build.yml/badge.svg)](https://github.com/agmcfarland/BEiGUIDE/actions/workflows/test-build.yml)
[![codecov](https://codecov.io/gh/agmcfarland/BEiGUIDE/graph/badge.svg?token=NPALNGNUFJ)](https://codecov.io/gh/agmcfarland/BEiGUIDE)
<!-- Badges end -->

# Background

`Base editor iGUIDE (BEiGUIDE)` allows for quantification of base edits in on and off-target cut sites identified by iGUIDE.

# Installation

```R
devtools::install_github('agmcfarland/BEiGUIDE')
```

# Running as a script

```sh
Rscript -e "BEiGUIDE::quantify_edits('path/to/base_directory', 'path/to/output_directory')"
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
```

# Parameters

```md
base_directory:
	Character. The base directory where input data is located.

output_directory:
	Character. The directory where output data will be saved.

analysis_name:
	Character. The name of the analysis (default: 'quantify_edits').

abundance_cutoff:
	Numeric. The minimum abundance cutoff for including an edit site (default: 3).

end_distance_from_cut_site:
	Numeric. The maximum distance from the cut site to consider bases (default: 20).

cut_site_start_distance_within_gRNA:
	Numeric. The distance within gRNA to start considering cut sites (default: 3).

cut_site_start_distance_outside_gRNA:
	Numeric. The distance outside gRNA to start considering cut sites (default: 3).

reference_genome_path
	Character. Path to a fasta file used to generate the base iGUIDE result (default: '').

editable_base
	Character. Base that can be edited by base editor (default: A).

expected_edit
	Character. The base that edtiable base will be turned into by base editor (default: G).

binomial_p_value_threshold
	Numeric. Significance threshold for binomial test of proportions to determine if percentage of edited base is significant (default: 0.05). 

binomial_direction
	Character. One of c('greater', 'less', 'two.sided'). The direction of the binomial test of proportions (default: 'greater').

n_processors:
	Numeric. The number of processors to use for parallel processing (default: 4).

overwrite:
	Logical. Whether to overwrite existing analysis output (default: TRUE).
```