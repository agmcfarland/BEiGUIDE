#' Filter BAM Alignments
#'
#' Takes as input a data frame of BAM alignments and filters based on edit site position and strand.
#'
#' @details
#' This function filters BAM alignments based on the specified edit site position and strand.
#'
#' The function uses different filtering criteria based on the strand of the edit site:
#' - For the '+' strand, alignments are filtered with flag '83', indicating the cut site is to the left of the PAM site.
#' - For the '-' strand, alignments are filtered with flag '99', indicating the cut site is to the right of the PAM site.
#'
#' Note: `pos` represents the 1-based start position of the alignment on the genome.
#' We use the genomic alignment `width` (not `qwidth`) to calculate the exact genomic end position (`pos + width - 1`).
#'
#' Outputted columns `edit_site_min` and `edit_site_max` are always relative to the plus strand of the genome.
#'
#' @param df_bam A data frame containing BAM alignments (must contain `cigar` and `width` columns).
#' @param specimen The specimen ID.
#' @param chromosome The chromosomal location of the edit site.
#' @param edit_site_position The position of the edit site.
#' @param edit_site_strand The strand of the edit site ('+' or '-').
#' @param allowed_aln_start_cut_site_PAM_distal The allowed alignment start with respect to the cut site, PAM distal. Default is 3.
#' @param allowed_aln_start_cut_site_PAM_proximal The allowed alignment start with respect to the cut site, PAM proximal Default is 3.
#'
#' @return A filtered data frame of BAM alignments.
#'
#' @export
#'
#' @import dplyr
filter_bam_alignments <- function(df_bam,
                                  specimen,
                                  chromosome,
                                  edit_site_position,
                                  edit_site_strand,
                                  allowed_aln_start_cut_site_PAM_distal = 3,
                                  allowed_aln_start_cut_site_PAM_proximal = 3) {

  # Initial filtering for specimen, chromosome
  df_bam_filtered <- df_bam %>%
    dplyr::filter(
      specimen_id == specimen,
      rname == chromosome
      )

  if (edit_site_strand == '-') {
    df_bam_filtered <- df_bam_filtered %>%
      dplyr::mutate(
        edit_site_max = edit_site_position + allowed_aln_start_cut_site_PAM_distal,
        edit_site_min = edit_site_position - allowed_aln_start_cut_site_PAM_proximal
      ) %>%
      dplyr::filter(
        flag == '99', # - strand cut with no indels
        # For minus strand, the 5' end of the read (pos) aligns to the DSB/cut site/dsODN incorporation site.
        pos >= edit_site_min,
        pos <= edit_site_max
      )
  }

  if (edit_site_strand == '+') {
    df_bam_filtered <- df_bam_filtered %>%
      dplyr::mutate(
        # Subtract 1 to account one-based coordinate system
        aln_pos_end = pos + qwidth - 1,
        edit_site_max = edit_site_position + allowed_aln_start_cut_site_PAM_proximal,
        edit_site_min = edit_site_position - allowed_aln_start_cut_site_PAM_distal
      ) %>%
      dplyr::filter(
        flag == '83', # + strand cut with no indels
        # For plus strand, the 3' end of the alignment (aln_pos_end) aligns to the DSB/cut site/dsODN incorporation site.
        aln_pos_end >= edit_site_min,
        aln_pos_end <= edit_site_max
      ) %>%
      dplyr::select(-aln_pos_end)
  }

  return(df_bam_filtered)
}
