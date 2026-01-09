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

  # 1. Filter BAM
  df_bam_filtered <- filter_bam_alignments(
    df_bam = df_bam,
    specimen = edit_site$description$specimen,
    chromosome = edit_site$description$chromosome,
    edit_site_position = edit_site$description$position,
    edit_site_strand = edit_site$description$strand,
    allowed_bases_within_gRNA = run_params$cut_site_start_distance_within_gRNA,
    allowed_bases_outside_gRNA = run_params$cut_site_start_distance_outside_gRNA
  )

  if (nrow(df_bam_filtered) == 0) {
    return(NULL) # Equivalent of 'next'
  }

  # 2. Process Alignments
  df_aln_pos <- aln_to_base_and_position_table(df_bam_filtered = df_bam_filtered)

  df_aln_pos_filtered <- calculate_base_percentages(
    df_aln_pos = df_aln_pos,
    edit_site_strand = edit_site$description$strand,
    edit_site_position = edit_site$description$position,
    bases_from_cut_site = run_params$end_distance_from_cut_site
  )

  # 3. Handle Coordinates & Reference
  edit_site$protospacer <- cut_site_to_protospacer_coordinates(
    cut_site = edit_site$description$position,
    strand = edit_site$description$strand,
    end_distance_from_cut_site = run_params$end_distance_from_cut_site,
    expected_cut_distance_from_pam = run_params$expected_cut_distance_from_pam,
    pam_length = 3
  )

  if (edit_site$description$strand == '+') {
    start_pos <- edit_site$protospacer$protospacer_end
    stop_pos <- edit_site$protospacer$PAM_start
  } else {
    start_pos <- edit_site$protospacer$PAM_start
    stop_pos <- edit_site$protospacer$protospacer_end
  }

  edit_site$genomic_reference <- extract_reference_bases(
    genomic_sequence = genome_sequence,
    chromosome = edit_site$description$chromosome,
    start = start_pos,
    stop = stop_pos
  )

  # 4. Significance Testing
  df_base_percentages_with_reference <- binomial_prop_edit_test(
    df_base_counts = df_aln_pos_filtered,
    p_value_threshold = run_params$binomial_p_value_threshold,
    alternative_ = run_params$binomial_direction
  ) %>%
    dplyr::left_join(edit_site$genomic_reference, by = 'position')

  # 5. Strand Correction
  if (edit_site$description$strand == '-') {
    df_base_percentages_with_reference$reference_base <- as.character(unlist(sapply(df_base_percentages_with_reference$reference_base, base_complement)))
    df_base_percentages_with_reference$base <- as.character(unlist(sapply(df_base_percentages_with_reference$base, base_complement)))
  }

  edit_site$significant_edits <- test_for_significant_edits(
    df_edit_significance = df_base_percentages_with_reference,
    editable_base = run_params$editable_base,
    expected_edit = run_params$expected_edit,
    pvalue_threshold = run_params$binomial_p_value_threshold
  )

  # 6. Final Formatting
  edit_site$base_composition <- df_base_percentages_with_reference %>%
    dplyr::select(position, position_depth, base, base_count, percentage, reference_base)

  edit_site$standardized_locus <- standardize_cut_site_direction(
    strand = edit_site$description$strand,
    relative_zero = edit_site$protospacer$protospacer_start,
    genomic_positions = edit_site$genomic_reference$position,
    bases = edit_site$genomic_reference$reference_base
  )

  # 7. Save and Clean up
  saveRDS(edit_site, file.path(run_params$analysis_output, 'edit_sites', paste0('edit_site_id_', edit_site$description$unique_edit_site_id, '.rds')))

  return(edit_site$description$unique_edit_site_id) # Return ID to track completion

  }
