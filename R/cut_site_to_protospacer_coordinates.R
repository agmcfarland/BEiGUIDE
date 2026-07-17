#' Map Cut Site to Protospacer and PAM Coordinates
#'
#' Calculates the genomic coordinates for the protospacer and the boundary 
#' PAM start coordinate in a 1-based inclusive coordinate system.
#'
#' @param cut_site Numeric. Genomic coordinate of the double-strand break.
#' @param strand Character. The strand orientation, either `'+'` or `'-'`.
#' @param protospacer_distance_from_cut_PAM_distal Numeric. Bases from cut site to PAM-distal end of spacer (default: `17`).
#' @param protospacer_distance_from_cut_PAM_proximal Numeric. Bases from cut site to PAM-proximal end of spacer (default: `3`).
#' @param pam_length Numeric. Length of the PAM sequence (default: `3`).
#'
#' @return A named list containing:
#' \describe{
#'   \item{PAM_start}{The boundary genomic coordinate of the PAM (furthest from the spacer).}
#'   \item{protospacer_start}{The genomic coordinate of the PAM-proximal end of the spacer.}
#'   \item{protospacer_end}{The genomic coordinate of the PAM-distal end of the spacer.}
#' }
#' @export
cut_site_to_protospacer_coordinates <- function(cut_site, strand,
                                                protospacer_distance_from_cut_PAM_distal = 17,
                                                protospacer_distance_from_cut_PAM_proximal = 3,
                                                pam_length = 3) {
  cut_site <- as.numeric(cut_site)

  if (strand == '+') {
    # Spacer PAM-proximal end (right of cut)
    protospacer_start <- cut_site + protospacer_distance_from_cut_PAM_proximal
    # Spacer PAM-distal end (left of cut)
    protospacer_end <- cut_site - (protospacer_distance_from_cut_PAM_distal - 1)

    # PAM is to the right; boundary is shifted right by the PAM length
    PAM_start <- protospacer_start + pam_length

  } else if (strand == '-') {
    # Spacer PAM-proximal end (left of cut)
    protospacer_start <- cut_site - protospacer_distance_from_cut_PAM_proximal
    # Spacer PAM-distal end (right of cut)
    protospacer_end <- cut_site + (protospacer_distance_from_cut_PAM_distal - 1)

    # PAM is to the left; boundary is shifted left by the PAM length
    PAM_start <- protospacer_start - pam_length

  } else {
    stop("Strand must be either '+' or '-'.")
  }

  return(list(
    'PAM_start' = PAM_start,
    'protospacer_start' = protospacer_start,
    'protospacer_end' = protospacer_end
  ))
}