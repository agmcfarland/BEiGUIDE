#' Characterize Editing Outcomes at a Specific Genomic Site
#'
#' This function processes a single edit site by filtering BAM alignments,
#' calculating base compositions, identifying significant editing events via
#' binomial testing, and saving the resulting data structure as an RDS file.
#'
#' @param edit_site_row A single-row data frame containing site metadata.
#' Must include columns: \code{specimen}, \code{chromosome}, \code{position},
#' \code{strand}, and \code{unique_edit_site_id}.
#'
#' @details
#' This function relies on several objects existing in the global environment.
#' \itemize{
#'   \item \code{df_bam}: The master data frame of BAM alignments.
#'   \item \code{run_params}: A list of configuration parameters (e.g., thresholds, distances).
#'   \item \code{genome_sequence}: The reference genome object used for base extraction.
#' }
#'
#' The function performs the following steps:
#' \enumerate{
#'   \item Spatial filtering of BAM alignments around the cut site.
#'   \item Conversion of alignments to base-position frequency tables.
#'   \item Coordinate annotation for protospacers and PAM sequences.
#'   \item Binomial significance testing for editing levels against a threshold.
#'   \item Strand-aware base complementation.
#'   \item Saving a serialized \code{.rds} object containing the full site characterization.
#' }
#'
#' @return Returns the \code{unique_edit_site_id} (character/numeric) of the processed
#' site if successful. Returns \code{NULL} if no alignments are found after filtering.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Example usage within a loop or mclapply:
#' result <- characterize_edit_site(df_edit_sites[1, ])
#' }
characterize_edit_site <- function(edit_site_row) {

  edit_site <- list(
    'description' = edit_site_row
  )

  # Filter BAM
  df_bam_filtered <- filter_bam_alignments(
    df_bam = df_bam,
    specimen = edit_site$description$specimen,
    chromosome = edit_site$description$chromosome,
    edit_site_position = edit_site$description$position,
    edit_site_strand = edit_site$description$strand,
    allowed_aln_start_cut_site_PAM_distal = run_params$allowed_aln_start_cut_site_PAM_distal,
    allowed_aln_start_cut_site_PAM_proximal = run_params$allowed_aln_start_cut_site_PAM_proximal
  )

  if (nrow(df_bam_filtered) == 0) {
    edit_site$fail <- 'no rows after bam filtering'
    saveRDS(edit_site, file.path(run_params$analysis_output, 'edit_sites', paste0('edit_site_id_', edit_site$description$unique_edit_site_id, '.rds')))
    return(NULL)
  }

  # Process Alignments
  df_aln_pos <- aln_to_base_and_position_table(df_bam_filtered = df_bam_filtered)


  # Standardize coordinates, add reference
  edit_site$protospacer <- cut_site_to_protospacer_coordinates(
    cut_site = edit_site$description$position,
    strand = edit_site$description$strand,
    protospacer_distance_from_cut_PAM_distal = run_params$expected_PAM_distal_cut_distance,
    protospacer_distance_from_cut_PAM_proximal = run_params$expected_PAM_proximal_cut_distance,
    pam_length = 3
  )

  absolute_start_stop_coords <- protospacer_coords_to_absolute_start_stop_coords(
    protospacer_coords = edit_site$protospacer,
    strand = edit_site$description$strand
    )

  edit_site$genomic_reference <- extract_reference_bases(
    genomic_sequence = genome_sequence,
    chromosome = edit_site$description$chromosome,
    start = absolute_start_stop_coords$start,
    stop = absolute_start_stop_coords$stop
  )

  # Calculate percentage of each base
  df_aln_pos_filtered <- calculate_base_percentages(
    df_aln_pos = df_aln_pos,
    start = absolute_start_stop_coords$start,
    stop = absolute_start_stop_coords$stop
  )

  # Significance Testing
  df_base_percentages_with_reference <- binomial_prop_edit_test(
    df_base_counts = df_aln_pos_filtered,
    null_proportion = run_params$null_proportion,
    alternative_ = run_params$binomial_direction
  ) %>%
    dplyr::left_join(edit_site$genomic_reference, by = 'position')

  # Strand Correction
  if (edit_site$description$strand == '-') {
    df_base_percentages_with_reference$reference_base <- as.character(unlist(sapply(df_base_percentages_with_reference$reference_base, base_complement)))
    df_base_percentages_with_reference$base <- as.character(unlist(sapply(df_base_percentages_with_reference$base, base_complement)))
  }

  edit_site$significant_edits <- test_for_significant_edits(
    df_edit_significance = df_base_percentages_with_reference,
    editable_base = run_params$editable_base,
    expected_edit = run_params$expected_edit,
    pvalue_threshold = run_params$binomial_p_value_threshold,
    adj_method = run_params$p_value_adj_method
  )

  # Final Formatting
  edit_site$base_composition <- df_base_percentages_with_reference %>%
    dplyr::select(position, position_depth, base, base_count, percentage, reference_base)

  edit_site$standardized_locus <- standardize_cut_site_direction(
    strand = edit_site$description$strand,
    relative_zero = edit_site$protospacer$protospacer_start,
    genomic_positions = edit_site$genomic_reference$position,
    bases = edit_site$genomic_reference$reference_base
  )

  # Save and Clean up
  saveRDS(edit_site, file.path(run_params$analysis_output, 'edit_sites', paste0('edit_site_id_', edit_site$description$unique_edit_site_id, '.rds')))

  return(edit_site$description$unique_edit_site_id) # Return ID to track completion

  }
