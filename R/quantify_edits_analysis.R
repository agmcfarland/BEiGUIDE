#' Quantify and Summarize Base Edits from a BEiGUIDE Run
#'
#' This function reads the output files from a BEiGUIDE quantify_edits run,
#' generates diagnostic plots for each edit site passing abundance filters,
#' collects significant base edits across sites, and creates run-level
#' summary tables.
#'
#' Specifically, it:
#' \itemize{
#'   \item Loads run parameters and edit site overview metadata.
#'   \item Filters edit sites by abundance filter.
#'   \item For each passing site, loads detailed results, generates plots
#'         (heatmap, on-target-only, depth, and stacked percentage), and saves them as PDFs.
#'   \item Collects all significant edits into a combined table, saves as CSV and RDS.
#'   \item Generates a run-level summary table showing counts of sites with on/off-target edits
#'         by annotation, and saves it as CSV and RDS.
#' }
#'
#' @param quantify_edits_output_path Character string. Path to the `quantify_edits` output directory
#'   from a BEiGUIDE run. This directory should contain the files
#'   `run_parameters.rds` and `edit_sites_overview.rds`.
#'
#' @return This function is used for its side effects:
#'   \itemize{
#'     \item Saved PDF diagnostic plots in `<analysis_output>/plots/`
#'     \item `significant_base_edits.csv` and `.rds` in `<analysis_output>/`
#'     \item `run_base_edit_summary.csv` and `.rds` in `<analysis_output>/`
#'   }
#' It does not return a value.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' quantify_edits_analysis('/path/to/quantify_edits')
#' }
quantify_edits_analysis <- function(quantify_edits_output_path) {

  # quantify_edits_output_path = '/data/BEiGUIDE/tests/quantify_edits'

  run_params <- readRDS(file.path(quantify_edits_output_path, 'run_parameters.rds'))

  df_edit_sites <- readRDS(file.path(quantify_edits_output_path, 'edit_sites_overview.rds'))

  df_edit_sites_pass_abundance <- df_edit_sites %>%
    dplyr::filter(pass_abundance_filter == T)

  df_base_edits <- data.frame()
  for (edit_site_row in split(df_edit_sites_pass_abundance, 1:nrow(df_edit_sites_pass_abundance))) {#break}

    edit_site_data <- file.path(run_params$analysis_output, 'edit_sites', paste0('edit_site_id_', edit_site_row$unique_edit_site_id, '.rds'))

    if (file.exists(edit_site_data)) {
      edit_site <- readRDS(edit_site_data)
    }

    plot_results <- generate_edit_site_diagnostic_plots(edit_site)

    ggsave(
      file.path(run_params$analysis_output, 'plots', paste0('heatmap_edit_site_', edit_site$description$unique_edit_site_id, '.pdf')),
      plot_results$heatmap,
      width = 8,
      height = 4
      )

    ggsave(
      file.path(run_params$analysis_output, 'plots', paste0('on_target_only_edit_site_', edit_site$description$unique_edit_site_id, '.pdf')),
      plot_results$on_target_only,
      width = 8,
      height = 4
    )

    ggsave(
      file.path(run_params$analysis_output, 'plots', paste0('depth_edit_site_', edit_site$description$unique_edit_site_id, '.pdf')),
      plot_results$depth,
      width = 8,
      height = 4
    )

    ggsave(
      file.path(run_params$analysis_output, 'plots', paste0('stacked_percentage_edit_site_', edit_site$description$unique_edit_site_id, '.pdf')),
      plot_results$stacked_percentage,
      width = 8,
      height = 4
    )

    df_base_edits <- rbind(
      df_base_edits,
      edit_site$significant_edits$significant %>%
        dplyr::mutate(unique_edit_site_id = edit_site$description$unique_edit_site_id)
    )
  }

  df_base_edits_extended <- df_base_edits %>%
    dplyr::rename(
      base_edit_position = position
    ) %>%
    dplyr::left_join(
      df_edit_sites,
      by = 'unique_edit_site_id'
    )

  write.csv(df_base_edits_extended, file.path(run_params$analysis_output, 'significant_base_edits.csv'), row.names = F)
  saveRDS(df_base_edits_extended, file.path(run_params$analysis_output, 'significant_base_edits.rds'))

  df_base_edits_mod <- df_base_edits %>%
    dplyr::group_by(unique_edit_site_id) %>%
    dplyr::mutate(
      on_target_edits = sum(on_target_edit) > 0,
      off_target_edits = sum(!on_target_edit) > 0
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(unique_edit_site_id, on_target_edits, off_target_edits) %>%
    dplyr::distinct()

  df_overview <- expand.grid(
    'unique_edit_site_id' = df_edit_sites$unique_edit_site_id,
    'on_target_edits' = F,
    'off_target_edits' = F
  ) %>%
    dplyr::filter(!unique_edit_site_id %in% df_base_edits_mod$unique_edit_site_id) %>%
    dplyr::bind_rows(df_base_edits_mod) %>%
    dplyr::left_join(
      df_edit_sites,
      by = 'unique_edit_site_id'
    )

  df_overview_summarized <- df_overview %>%
    dplyr::group_by(annotation) %>%
    dplyr::mutate(
      edit_sites = dplyr::n(),
      edit_sites_passing_abundance_filter = sum(pass_abundance_filter),
      edit_sites_with_on_target_edits = sum((on_target_edits == T & off_target_edits == F)),
      edit_sites_with_off_target_edits = sum((on_target_edits == F & off_target_edits == T)),
      edit_sites_with_on_and_off_target_edits = sum((on_target_edits == T & off_target_edits == T))
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(annotation, target.seq, edit_sites, edit_sites_passing_abundance_filter, edit_sites_with_on_target_edits, edit_sites_with_off_target_edits, edit_sites_with_on_and_off_target_edits) %>%
    dplyr::distinct() %>%
    dplyr::arrange(annotation)

  write.csv(df_overview_summarized, file.path(run_params$analysis_output, 'run_base_edit_summary.csv'), row.names = F)
  saveRDS(df_overview_summarized, file.path(run_params$analysis_output, 'run_base_edit_summary.rds'))


  }
