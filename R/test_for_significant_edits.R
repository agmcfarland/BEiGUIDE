#' Test for Significant Edits
#'
#' Identifies edited bases that differ from the reference base. Filters out 
#' unedited bases, missing calls (N), and zero-count bases, calculates multiple-testing
#' corrected p-values on the remaining edits, and flags on-target vs. off-target events.
#'
#' @param df_edit_significance A data frame containing at least the columns:
#'   `position`, `base`, `reference_base`, `position_depth`, `base_count`,
#'   `p_value`, `conf_low`, `conf_high`, and `percentage`.
#' @param editable_base A character indicating the reference base considered editable
#'   (default is `'A'`).
#' @param expected_edit A character indicating the expected edited base
#'   (default is `'G'`).
#' @param pvalue_threshold Numeric. Significance threshold for adjusted p-values
#'   (default is `0.05`).
#' @param adj_method A character specifying the p-value adjustment method passed to
#'   p.adjust() (default is `'BH'`).
#'
#' @return A named list with two data frames:
#' \describe{
#'   \item{`all`}{All non-reference base calls with on-target flag and adjusted p-values.}
#'   \item{`significant`}{Subset of `all` where `adj_p_value` <= `pvalue_threshold`.}
#' }
#'
#' @export
test_for_significant_edits <- function(df_edit_significance, 
                                       editable_base = 'A', 
                                       expected_edit = 'G', 
                                       pvalue_threshold = 0.05, 
                                       adj_method = 'BH') {

  # Subset to required columns
  df_edit_significance <- df_edit_significance[, c(
    'position', 'base', 'reference_base', 'position_depth', 
    'base_count', 'p_value', 'conf_low', 'conf_high', 'percentage'
  )]

  # Filter to actual edit events and calculate multiple testing correction
  df_edit_significance <- df_edit_significance %>%
    dplyr::filter(
      base != reference_base,
      base != 'N',
      base_count > 0
    ) %>%
    dplyr::mutate(
      on_target_edit = ifelse(reference_base == editable_base & base == expected_edit, TRUE, FALSE),
      adj_p_value = stats::p.adjust(p_value, method = adj_method)
    )

  return(
    list(
      'all' = df_edit_significance,
      'significant' = df_edit_significance %>% dplyr::filter(adj_p_value <= pvalue_threshold)
    )
  )
}