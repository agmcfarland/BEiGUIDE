#' Binomial Proportion Edit Test
#'
#' Performs a binomial test for each base at each genomic position to assess whether
#' the observed base count is significantly different from an expected baseline
#' proportion (default 5%). For each base, this function reports the p-value and
#' confidence interval from the binomial test, along with the percentage of reads
#' supporting that base.
#'
#' @param df_base_counts A data frame containing base counts and coverage depth for each position.
#' It must include the columns:
#' \itemize{
#'   \item `position` — Genomic position (numeric or integer).
#'   \item `base` — Base identity as a character (A, C, G, T, or N).
#'   \item `base_count` — Number of reads supporting that base.
#'   \item `position_depth` — Total read depth at that position.
#' }
#' @param p_value_threshold A numeric value specifying the expected proportion under the null
#' hypothesis (default is 0.05).
#' @param alternative_ A character string specifying the alternative hypothesis to test.
#' One of `"two.sided"`, `"less"`, or `"greater"`. Default is `"two.sided"`.
#'
#' @return A data frame with the following columns:
#' \itemize{
#'   \item `position` — The genomic position tested.
#'   \item `base` — The base tested.
#'   \item `base_count` — Number of reads supporting that base.
#'   \item `position_depth` — Total depth at that position.
#'   \item `p_value` — The p-value from the binomial test.
#'   \item `conf_low` — Lower bound of the binomial confidence interval.
#'   \item `conf_high` — Upper bound of the binomial confidence interval.
#'   \item `percentage` — Percentage of reads supporting the base (0–100).
#' }
#'
#' @details
#' If any of the bases A, C, G, T, or N are missing from the input at a given position,
#' they are added with a count of 0 before testing. The binomial test is run using
#' \code{stats::binom.test} on each base at each position.
#'
#' @import dplyr tidyr
#'
#' @export
binomial_prop_edit_test <- function(df_base_counts, p_value_threshold = 0.05, alternative_ = c('two.sided', 'less', 'greater')) {
  if (!alternative_ %in% c('two.sided', 'less', 'greater')) {
    stop('Specify one of two.sided, less, greater')
  }

  df_base_counts <- df_base_counts[c('position', 'base', 'base_count', 'position_depth')]

  df_edit_evidence <- df_base_counts %>%
    # Fill missing base counts
    dplyr::group_by(position) %>%
    tidyr::pivot_wider(names_from = 'base', values_from = 'base_count', values_fill = 0) %>%
    dplyr::ungroup()

  # Add columns in case that base is not in the data and not available for pivot_wider to act on
  if (!'N' %in% colnames(df_edit_evidence)) {
    df_edit_evidence$`N` <- 0
  }
  if (!'A' %in% colnames(df_edit_evidence)) {
    df_edit_evidence$`A` <- 0
  }
  if (!'T' %in% colnames(df_edit_evidence)) {
    df_edit_evidence$`T` <- 0
  }
  if (!'C' %in% colnames(df_edit_evidence)) {
    df_edit_evidence$`C` <- 0
  }
  if (!'G' %in% colnames(df_edit_evidence)) {
    df_edit_evidence$`G` <- 0
  }

  # Convert back to long form
  df_edit_evidence <- df_edit_evidence %>%
    tidyr::pivot_longer(cols = c(`A`, `C`, `G`, `T`, `N`), names_to = 'base', values_to = 'base_count')

  # Test whether edit is significantly greater than 5%
  df_edit_evidence <- df_edit_evidence %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      binom_test = list(stats::binom.test(
        x = base_count,
        n = position_depth,
        p = p_value_threshold,
        alternative = alternative_
      )),
      p_value = binom_test$p.value,
      conf_low = binom_test$conf.int[1],
      conf_high = binom_test$conf.int[2]
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(-binom_test)

   # Add percentage column
   df_edit_evidence <- df_edit_evidence %>%
    dplyr::mutate(percentage = 100 * (base_count/position_depth))

  return(df_edit_evidence)
  }
