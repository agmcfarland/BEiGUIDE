#' Standardize Cut Site Orientation
#'
#' Adjusts genomic positions and base identities relative to a specified genomic position
#' (relative zero) according to the provided strand orientation. All relative positions read in the
#' 5'->3' direction.
#' For the positive strand, relative positions are computed as `relative_zero - genomic_position`. Genomic positions upstream
#' of the relative zero will have positive relative position values and downstream will be negative.
#' For the negative strand, relative positions are computed as
#' `genomic_position - relative_zero`, and `bases are complemented`. Genomic positions downstream of the relative zero
#' will have positive relative position values and upstream will be negative.
#'
#' @param strand Character. The strand orientation of the cut site (`'+'` or `'-'`).
#' @param relative_zero Numeric. The genomic coordinate representing the zero relative position (position zero).
#' @param genomic_positions Numeric vector. The genomic positions of each base to standardize.
#' @param bases Character vector. The base calls corresponding to each genomic position.
#'
#' @return A `data.frame` with columns:
#' \describe{
#'   \item{genomic_position}{Original genomic position of each base.}
#'   \item{base}{Base at the given position (reverse-complemented if strand is `'-'`).}
#'   \item{relative_position}{Position relative to the cut site. Is a factor.}
#' }
#'
#' @examples
#' \dontrun{
#' genomic_positions <- 100:105
#' bases <- c("A", "T", "G", "C", "C", "A")
#' standardize_cut_site_direction(
#'   strand = '+',
#'   relative_zero = 103,
#'   genomic_positions = genomic_positions,
#'   bases = bases
#' )
#'
#' standardize_cut_site_direction(
#'   strand = '-',
#'   relative_zero = 103,
#'   genomic_positions = genomic_positions,
#'   bases = bases
#' )
#' }
#'
#' @export
standardize_cut_site_direction <- function(
    strand,
    relative_zero,
    genomic_positions,
    bases
) {
  df_standard_loci <- data.frame(
    'genomic_position' = genomic_positions,
    'base' = bases
  )

  if (strand == '+') {
    df_standard_loci <- df_standard_loci %>%
      dplyr::mutate(
        relative_position = relative_zero - genomic_position
      )
  } else {
    df_standard_loci <- df_standard_loci %>%
      dplyr::mutate(
        relative_position = genomic_position - relative_zero
      )
    df_standard_loci$base <- stringr::str_split(as.character(Biostrings::complement(Biostrings::DNAString(paste(bases, collapse = '')))), '')[[1]]
  }

  df_standard_loci$relative_position <- factor(df_standard_loci$relative_position, levels = sort(df_standard_loci$relative_position, decreasing = T))

  return(df_standard_loci)
}
