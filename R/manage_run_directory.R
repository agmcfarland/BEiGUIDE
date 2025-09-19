#' Manage Run Directory for quantify_edits
#'
#' This function manages the run directory for an analysis by ensuring the output directory exists, handling existing analysis output directories based on the `overwrite` parameter, and creating the analysis output directory if it does not exist.
#'
#' @param run_params A data frame containing run parameters, including `output_directory`, `analysis_output`, and `overwrite`.
#'
#' @return None. The function manages the directory structure and throws errors if conditions are not met.
#'
#' @examples
#' \dontrun{
#' run_params <- data.frame(
#'   output_directory = "path/to/output_directory",
#'   analysis_output = "path/to/output_directory/analysis_name",
#'   overwrite = TRUE
#' )
#' manage_run_directory(run_params)
#' }
manage_run_directory <- function(run_params) {
  # browser()
  if (!dir.exists(run_params$output_directory)) {
    dir.create(run_params$output_directory, recursive = TRUE)
    print(paste0('making base output directory: ', run_params$output_directory))
  }

  if (dir.exists(run_params$analysis_output)) {
    if (run_params$overwrite) {
      unlink(run_params$analysis_output, recursive = TRUE)
    } else {
      stop(paste('analysis directory exists and overwrite is set to FALSE ', run_params$analysis_output))
    }
  }

  if (!dir.exists(run_params$analysis_output)) {
    print(paste0('making analysis output directory: ', run_params$output_directory))
    dir.create(run_params$analysis_output)
  }
  }

