#' Binomial Proportion Edit Test
#'
#' Performs a binomial test for each base at each genomic position to assess whether
#' the observed base count is significantly different from an expected baseline
#' proportion (default 0.1%).
#'
#' @param df_base_counts A data frame containing columns: `position`, `base`,
#'   `base_count`, and `position_depth`.
#' @param null_proportion Numeric. The expected background/error rate under the
#'   null hypothesis (default: `0.001`).
#' @param alternative_ One of `"greater"`, `"two.sided"`, or `"less"`. Default is `"greater"`.
#'
#' @return A data frame containing the tested positions, bases, counts, depths,
#'   binomial p-values, confidence intervals, and percentages.
#'
#' @export
binomial_prop_edit_test <- function(df_base_counts,
                                    null_proportion = 0.001,
                                    alternative_ = c('greater', 'two.sided', 'less')) {

  alternative_ <- match.arg(alternative_)

  df_base_counts <- df_base_counts[, c('position', 'base', 'base_count', 'position_depth')]

  # Complete the grid so all bases (A, C, G, T, N) are represented per position
  # Position depth of the base is also recorded.
  # Results in a long-form dataframe
  df_edit_evidence <- df_base_counts %>%
    dplyr::group_by(position) %>%
    tidyr::complete(base = c('A', 'C', 'G', 'T', 'N'), fill = list(base_count = 0)) %>%
    dplyr::mutate(position_depth = max(position_depth, na.rm = T)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(percentage = 100 * (base_count / position_depth))

  # Run binomial test on all base combinations
  df_edit_evidence <- df_edit_evidence %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      binom_res = list(stats::binom.test(
        x = base_count,
        n = position_depth,
        p = null_proportion,
        alternative = alternative_
      )),
      p_value = binom_res$p.value,
      conf_low = binom_res$conf.int[1],
      conf_high = binom_res$conf.int[2]
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(-binom_res)

  return(df_edit_evidence)
}
