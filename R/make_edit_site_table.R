#' Make Edit Sites Table
#'
#' Description: Makes Edit sites dataframe. This is the experimental object that is used to identify
#' iGUIDE edit sites in alignment data
#'
#' @param ft_data_table A data table containing information about edit sites.
#' @param spec_info_combo_overview_table A data table containing information about specimen and annotations.
#' @param abundance_cutoff Numeric value specifying the abundance cutoff for filtering.
#'
#' @return A data frame containing information about edit sites.
#'
#' @import dplyr stringr
#'
#' @export
make_edit_site_table <- function(ft_data_table, spec_info_combo_overview_table, abundance_cutoff = 3) {
  # Further modify the edit site dataframe
  df_edit_sites <- ft_data_table %>%
    dplyr::mutate(pass_abundance_filter = abund > abundance_cutoff) %>%
    dplyr::left_join(
      y = spec_info_combo_overview_table %>%
        dplyr::select(annotation, specimen),
      by = 'annotation') %>%
    split_edit_site_id() %>%
    dplyr::arrange(dplyr::desc(abund)) %>%
    dplyr::mutate(
      unique_edit_site_id = seq(1, dplyr::n())
    )

  return(df_edit_sites)
}
