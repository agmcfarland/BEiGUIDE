#' Test for Significant Edits
#'
#' Given a data frame of per-base edit significance statistics, this function
#' identifies edited bases that differ from the reference. If the strand is negative
#' then the sequence of the strand is complemented. Edits (e.g. A→G) are detected relative to the 5'->3' orientation.
#' Function will flags on-target and off-target edits where an on-target matches the editable base to expected edit
#' while off target is any other combination where an observed base and reference base pair differ.
#' Function returns both the full filtered data frame where each row is a base-reference base pair that are different
#' and a filtered dataframe where edits have a p_value lower than the provided threshold.
#'
#' @param df_edit_significance A data frame containing at least the columns:
#'   `position`, `base`, `reference_base`, `position_depth`, `base_count`,
#'   `p_value`, `conf_low`, `conf_high`, and `percentage`.
#' @param strand A character (`'+'` or `'-'`) indicating the strand orientation.
#'   If `'-'`, the reference and observed bases will be reverse-complemented.
#' @param editable_base A character indicating the reference base considered editable
#'   (default is `'A'`).
#' @param expected_edit A character indicating the expected edited base
#'   (default is `'G'`).
#' @param pvalue_threshold Numeric. Significance threshold for p-values
#'   (default is `0.05`).
#'
#' @return A named list with two data frames:
#' \describe{
#'   \item{`all`}{All non-reference base calls with on-target flag.}
#'   \item{`significant`}{Subset of `all` where `p_value` ≤ `pvalue_threshold`.}
#' }
#'
#' @details
#' This function is useful for downstream filtering of base editing results.
#' For negative-strand data, bases are converted to their complements before
#' comparison to the reference.
#'
#' @examples
#' \dontrun{
#' df <- data.frame(
#'   position = 1:3,
#'   base = c('G', 'A', 'T'),
#'   reference_base = c('A', 'A', 'T'),
#'   position_depth = c(20, 20, 20),
#'   base_count = c(10, 2, 1),
#'   p_value = c(0.01, 0.2, 0.5),
#'   conf_low = c(0.2, 0.01, 0.0),
#'   conf_high = c(0.8, 0.3, 0.1),
#'   percentage = c(50, 10, 5)
#' )
#'
#' test_for_significant_edits(df, strand = '+')
#' }
#'
#' @export
test_for_significant_edits <- function(df_edit_significance, editable_base = 'A', expected_edit = 'G', pvalue_threshold = 0.05) {

  df_edit_significance <- df_edit_significance[c('position', 'base', 'reference_base', 'position_depth', 'base_count', 'p_value', 'conf_low', 'conf_high', 'percentage')]

  df_edit_significance <- df_edit_significance %>%
    dplyr::filter(base != reference_base) %>%
    dplyr::mutate(
      on_target_edit = ifelse(reference_base == editable_base & base == expected_edit, TRUE, FALSE)
    )

  return(
    list(
      'all' = df_edit_significance,
      'significant' = df_edit_significance %>% dplyr::filter(p_value <= pvalue_threshold)
    )
  )

  }
