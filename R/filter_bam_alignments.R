#' Filter BAM Alignments
#'
#' Takes as input a data frame of BAM alignments and filters based on edit site position and strand.
#'
#' @details
#' This function filters BAM alignments based on the specified edit site position and strand.
#' It removes alignments containing soft clipping ('S') or deletions ('D') in the CIGAR string.
#'
#' The function uses different filtering criteria based on the strand of the edit site:
#' - For the '+' strand, alignments are filtered with flag '83', indicating the cut site is to the left of the PAM site.
#' - For the '-' strand, alignments are filtered with flag '99', indicating the cut site is to the right of the PAM site.
#'
#' Note: `pos` represents the start position of the alignment and is always on the left-hand side of the aligned read,
#' regardless of the strand it mapped to.
#'
#' Outputted columns `edit_site_min` and `edit_site_max` are always relative to the plus strand of the genome. Therefore min is always less than
#' max regardless of strand being edited and always relative to the plus strand.
#'
#' @param df_bam A data frame containing BAM alignments.
#' @param specimen The iGUIDE specimen ID
#' @param chromosome The chromosomal loocation of the edit site.
#' @param edit_site_position The position of the edit site.
#' @param edit_site_strand The strand of thee edit site ('+' or '-').
#' @param allowed_bases_within_gRNA The number of base pairs allowed within the gRNA region. Default is 3.
#' @param allowed_bases_outside_gRNA The number of base pairs allowed outside the gRNA region. Default is 3.
#'
#' @return A filtered data frame of BAM alignments.
#'
#' @export
#'
#' @import dplyr
filter_bam_alignments <- function(df_bam, specimen, chromosome, edit_site_position, edit_site_strand, allowed_bases_within_gRNA = 3, allowed_bases_outside_gRNA = 3) {

  df_bam_filtered <- df_bam %>%
    dplyr::filter(
      specimen_id == specimen,
      rname == chromosome
      )

  if (edit_site_strand == '-') {
    df_bam_filtered <- df_bam_filtered %>%
      dplyr::mutate(
        edit_site_max = edit_site_position + allowed_bases_within_gRNA,
        edit_site_min = edit_site_position - allowed_bases_outside_gRNA
      ) %>%
      dplyr::filter(
        flag == '99', # - strand cut
        pos >= edit_site_min & pos <= edit_site_max,
        pos + qwidth > edit_site_position
      )
  }

  if (edit_site_strand == '+') {
    df_bam_filtered <- df_bam_filtered %>%
      dplyr::mutate(
        aln_pos_end = pos + qwidth,
        edit_site_max = edit_site_position + allowed_bases_outside_gRNA,
        edit_site_min = edit_site_position - allowed_bases_within_gRNA
      ) %>%
      dplyr::filter(
        flag == '83', # + strand cut
        aln_pos_end >= edit_site_min & aln_pos_end <= edit_site_max,
        pos < edit_site_position
      ) %>%
      dplyr::select(-aln_pos_end)
  }

  return(df_bam_filtered)
}
