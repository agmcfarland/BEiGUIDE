#' Quantify Edits in Sequencing Data
#'
#' This function quantifies edits in sequencing data by:
#' \itemize{
#'   \item Parsing and filtering aligned reads from BAM files.
#'   \item Identifying edit sites based on read abundance thresholds.
#'   \item Calculating base composition percentages at each edit site.
#'   \item Performing binomial tests to determine statistical significance of edits.
#' }
#' Results are written to disk as RDS and CSV files along with logs and per-site plots.
#'
#' @param base_directory Character. Base directory where input data is located.
#' @param output_directory Character. Directory where output data will be saved.
#' @param analysis_name Character. Name for this analysis (default: `'quantify_edits'`).
#' @param abundance_cutoff Numeric. Minimum read count for including an edit site (default: `3`).
#' @param expected_cut_distance_from_pam Numeric. Distance from PAM where the cut site is expected (default: `3`).
#' @param end_distance_from_cut_site Numeric. Maximum distance from the cut site to consider bases (default: `20`).
#' @param cut_site_start_distance_within_gRNA Numeric. Distance within the gRNA to start considering cut sites (default: `3`).
#' @param cut_site_start_distance_outside_gRNA Numeric. Distance outside the gRNA to start considering cut sites (default: `3`).
#' @param reference_genome_path Character. Path to a FASTA file used to retrieve reference genome sequence (default: `''`).
#' @param editable_base Character. The original base that is targeted for editing (default: `'A'`).
#' @param expected_edit Character. The expected base after editing (default: `'G'`).
#' @param binomial_p_value_threshold Numeric. P-value threshold for determining significance using a binomial test (default: `0.05`).
#' @param binomial_direction Character. Direction of the binomial test: `'greater'`, `'less'`, or `'two-sided'` (default: `'greater'`).
#' @param n_processors Numeric. Number of CPU cores to use for parallel processing of BAM files (default: `4`).
#' @param n_edit_site_processors Numeric. Number of CPU cores to use for parallel processing of edit sites with significant levels of base editing.
#' @param overwrite Logical. Whether to overwrite existing analysis output (default: `TRUE`).
#'
#' @return This function does not return an object. It writes the following to `output_directory/analysis_name`:
#' \itemize{
#'   \item `run_parameters.csv` / `run_parameters.rds` — analysis configuration
#'   \item `edit_sites_overview.csv` / `edit_sites_overview.rds` — table of edit site metadata
#'   \item One directory per edit site containing processed alignments, base composition tables, statistical results, and plots
#'   \item A log file of all steps and parameters used
#' }
#'
#' @import logr dplyr stringr BSgenome GenomicRanges Biostrings parallel
#'
#' @examples
#' \dontrun{
#' quantify_edits(
#'   base_directory = "path/to/base_directory",
#'   output_directory = "path/to/output_directory",
#'   analysis_name = "my_analysis",
#'   abundance_cutoff = 5,
#'   expected_cut_distance_from_pam = 3,
#'   end_distance_from_cut_site = 25,
#'   cut_site_start_distance_within_gRNA = 4,
#'   cut_site_start_distance_outside_gRNA = 4,
#'   reference_genome_path = "path/to/genome.fasta",
#'   editable_base = "C",
#'   expected_edit = "T",
#'   binomial_p_value_threshold = 0.01,
#'   binomial_direction = "two-sided",
#'   n_processors = 6,
#'   n_edit_site_processors = 3,
#'   overwrite = FALSE
#' )
#' }
#' @export

quantify_edits <- function(
		base_directory,
		output_directory,
		analysis_name = 'quantify_edits',
		abundance_cutoff = 3,
		expected_cut_distance_from_pam = 3,
		end_distance_from_cut_site = 20,
		cut_site_start_distance_within_gRNA = 3,
		cut_site_start_distance_outside_gRNA = 3,
		reference_genome_path = '',
		editable_base = 'A',
		expected_edit = 'G',
		binomial_p_value_threshold = 0.05,
		binomial_direction = 'greater',
		n_processors = 4,
		n_edit_site_processors = 3,
		overwrite = TRUE
) {

  library(BSgenome)
  library(logr)

  run_params <- data.frame(
  	'base_directory' = base_directory,
  	'output_directory' = output_directory,
  	'analysis_name' = analysis_name,
  	'abundance_cutoff' = abundance_cutoff,
  	'expected_cut_distance_from_pam' = expected_cut_distance_from_pam,
  	'end_distance_from_cut_site' = end_distance_from_cut_site,
  	'cut_site_start_distance_within_gRNA' = cut_site_start_distance_within_gRNA,
  	'cut_site_start_distance_outside_gRNA' = cut_site_start_distance_outside_gRNA,
  	'reference_genome_path' = reference_genome_path,
  	'editable_base' = editable_base,
  	'expected_edit' = expected_edit,
  	'binomial_p_value_threshold' = binomial_p_value_threshold,
  	'binomial_direction' = binomial_direction,
  	'n_processors' = n_processors,
  	'n_edit_site_processors' = n_edit_site_processors,
  	'overwrite' = overwrite
  )

  ## Troubleshooting start ##
  # run_params <- data.frame(
  #   'base_directory' = '/data/BEiGUIDE/data-raw/230705_MN01490_0144_A000H5KVFN',
  #   'output_directory' = '/data/BEiGUIDE/tests',
  #   'analysis_name' = 'quantify_edits',
  #   'abundance_cutoff' = 5,

  #   'expected_cut_distance_from_pam' = 3,
  #   'end_distance_from_cut_site' = 20,
  #   'cut_site_start_distance_within_gRNA' = 3,
  #   'cut_site_start_distance_outside_gRNA' = 3,
  #   'reference_genome_path' = '/data/iGUIDE/genomes/hg38.fasta',

  #   'editable_base' = 'A',
  #   'expected_edit' = 'G',

  #   'binomial_p_value_threshold' = 0.15,
  #   'binomial_direction' = 'greater',

  #   'n_processors' = 2,
  #   'overwrite' = TRUE
  # )
  ## Troubleshooting End ##

  run_params$analysis_output <- file.path(run_params$output_directory, run_params$analysis_name)
  manage_run_directory(run_params = run_params)
  saveRDS(run_params, file.path(run_params$analysis_output, 'run_parameters.rds'))
  write.csv(run_params, file.path(run_params$analysis_output, 'run_parameters.csv'), row.names = F)

  dir.create(file.path(run_params$analysis_output, 'edit_sites'))
  dir.create(file.path(run_params$analysis_output, 'plots'))

  logr::log_open(file_name = file.path(run_params$analysis_output, 'BEiGUIDE_quantify_edits'), logdir = FALSE)

  logr::log_print(run_params)

  df_ft_data <- pull_ft_data_tables(base_directory = run_params$base_directory)

  df_annotations <- pull_combo_overview_table(base_directory = run_params$base_directory)

  df_edit_sites <- make_edit_site_table(
  	ft_data_table = df_ft_data,
  	spec_info_combo_overview_table = df_annotations,
  	abundance_cutoff = run_params$abundance_cutoff)

  saveRDS(df_edit_sites, file.path(run_params$analysis_output, 'edit_sites_overview.rds'))
  write.csv(df_edit_sites, file.path(run_params$analysis_output, 'edit_sites_overview.csv'), row.names = F)

  bam_files <- list_bam_files(base_directory = run_params$base_directory)

  logr::log_print(bam_files)

  df_bam <- parallel_bam_file_list_to_table(
  	bam_file_list = bam_files,
  	scan_param_what_list = c('qname', 'rname', 'strand', 'pos', 'qwidth', 'seq', 'cigar', 'flag'),
  	number_of_cpu = run_params$n_processors
  )

  # Filter once for speed.
  df_bam <- df_bam %>%
    dplyr::filter(
      !stringr::str_detect(cigar, 'I'),
      !stringr::str_detect(cigar, 'S'),
      !stringr::str_detect(cigar, 'D')
    )

  genome_sequence <- load_reference_genome(run_params$reference_genome_path)

  df_edit_sites_pass_abundance <- df_edit_sites %>%
  	dplyr::filter(pass_abundance_filter == T)

  # edit_sites to process
  edit_site_list <- split(
  	df_edit_sites_pass_abundance,
  	1:nrow(df_edit_sites_pass_abundance)
  )

  base::assign("df_bam", df_bam, envir = .GlobalEnv)
  base::assign("run_params", run_params, envir = .GlobalEnv)
  base::assign("genome_sequence", genome_sequence, envir = .GlobalEnv)

  logr::log_print(paste("Starting parallel characterization of", length(edit_site_list), "sites using", run_params$n_edit_site_processors, "cores."))

	results <- parallel::mclapply(
		edit_site_list,
		characterize_edit_site,
		mc.cores = run_params$n_edit_site_processors
	)

	logr::log_print('Finished')

	logr::log_close()

	}



