#' Convert a cut site position to protospacer and PAM coordinates
#'
#' Given a genomic cut site and strand orientation, this function calculates the
#' coordinates of the protospacer start, protospacer end, and PAM start.
#'
#' @param cut_site Numeric. The cut site position (1-based).
#' @param strand Character. The strand orientation, either `'+'` or `'-'`.
#' @param end_distance_from_cut_site Numeric. Distance from the cut site to the end of the protospacer. Default is 20.
#' @param expected_cut_distance_from_pam Numeric. Distance between the cut site and the end of expected PAM site. Default is 3.
#' @param pam_length Numeric. Length of the PAM sequence. Default is 3.
#'
#' @description Nomenclature visual: 5'->protospacer_end--------cut_site---protospacer_startPAM_end-PAM_start->3'
#'
#' @return A named list with elements:
#' \describe{
#'   \item{PAM_start}{The calculated start position of the PAM.}
#'   \item{protospacer_start}{The calculated start position of the protospacer.}
#'   \item{protospacer_end}{The calculated end position of the protospacer.}
#' }
#'
#' @examples
#' \dontrun{
#' cut_site_to_protospacer_coordinates(100, '+')
#' cut_site_to_protospacer_coordinates(100, '-', end_distance_from_cut_site = 25)
#' }
#'
#' @export
cut_site_to_protospacer_coordinates <- function(cut_site, strand,
                                                end_distance_from_cut_site = 20,
                                                expected_cut_distance_from_pam = 3,
                                                pam_length = 3) {
  cut_site <- as.numeric(cut_site)

  if (strand == '+') {
    protospacer_start <- cut_site + expected_cut_distance_from_pam
    PAM_start <- protospacer_start + pam_length
    protospacer_end <- cut_site - end_distance_from_cut_site
  } else {
    protospacer_start <- cut_site - expected_cut_distance_from_pam
    PAM_start <- protospacer_start - pam_length
    protospacer_end <- cut_site + end_distance_from_cut_site
  }
  return(list(
    'PAM_start' = PAM_start,
    'protospacer_start' = protospacer_start,
    'protospacer_end' = protospacer_end
  ))
}
