#' Calculate Base Percentages
#'
#' Uses input from or `aln_to_base_and_position_tables()`  and returns a table with counts and percentages of each base at each position.
#'
#' @param df_aln_pos A data frame with base, position, and qname. Must already be sorted down to the chromosome level to be used here.
#' @param start Start minimum position. Desired positions must be greater than this value.
#' @param stop Stop minimum position. Desired positions must be less than this value.
#'
#' @details `Start` must be less than `stop.` Use `cut_site_protospacer_coords()` and `protospacer_coords_to_absolute_start_stop_coords()` to calculate `start` and `stop` if desired.
#'
#' @return A data frame with counts and percentages of each base at each position. The number of aligned reads (depth) for each position is also reported.
#'
#' @export
#'
#' @import dplyr
calculate_base_percentages <- function(df_aln_pos, start, stop) {

  testthat::expect_true(start < stop)

  df_aln_pos_filtered <- df_aln_pos %>%
    dplyr::filter(
      position >= start & position <= stop) %>%
    dplyr::group_by(position, base) %>%
    dplyr::mutate(base_count = dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::select(position, base, base_count) %>%
    base::unique() %>%
    dplyr::group_by(position) %>%
    dplyr::mutate(position_depth = sum(base_count)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(percentage = 100 * (base_count/position_depth))

  return(df_aln_pos_filtered)
}
