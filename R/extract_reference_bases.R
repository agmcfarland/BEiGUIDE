#' Extract Reference Bases for Edit Sites
#'
#' This function extracts the genomic sequence for a defined chromosome with start and stop coordinates. Does not consider strand, assumes coordinates are in a plus strand orientation relative to the provided genomic sequence.
#'
#' @param genomic_sequence A `BSgenome` object containing the genomic sequence from which reference bases will be extracted.
#' @param chromosome Chromosome of edit site
#' @param start Start position for edit site. Relative to + strand of genomic sequence.
#' @param stop Stop position edit site. Relative to + strand of genomic sequence.
#'
#' @return table with a column with the reference base and its genomic location
#'
#' @export
#'
#' @import Biostrings BSgenome dplyr
extract_reference_bases <- function(genomic_sequence, chromosome, start, stop) {

  gr <- GenomicRanges::GRanges(chromosome, IRanges::IRanges(seq(start, stop)), strand = '*')

  reference_bases <- Biostrings::getSeq(genomic_sequence, gr)

  df_positions <- data.frame(
    'position' = seq(start, stop, by = 1),
    'reference_base' = stringr::str_split(as.character(reference_bases), '')[[1]]
  )

  return(df_positions)
#
#   df_base_percentages_inputs <- df_base_percentages %>%
#     # dplyr::group_by(specimen, annotation, edit.site) %>%
#     dplyr::mutate(
#       start = min(position),
#       stop = max(position),
#       chromosome = edit_site_chromosome,
#     ) %>%
#     dplyr::ungroup() %>%
#     dplyr::select(specimen, annotation, edit_site_target.seq, start, stop, chromosome, edit_site_strand, edit.site) %>%
#     base::unique()
#
#
#   grange_base_percentages <- GenomicRanges::makeGRangesFromDataFrame(df_base_percentages_inputs, keep.extra.columns = TRUE, ignore.strand = TRUE)
#
#   reference_bases <- Biostrings::getSeq(genomic_sequence, grange_base_percentages)
#
#   df_all_reference_bases <- do.call(rbind, lapply(1:nrow(df_base_percentages_inputs), function (x) {
#
#     df_temp_reference_bases <- data.frame('reference_base' = unlist(reference_bases[x]))
#
#     df_temp_reference_bases$reference_position = seq(df_base_percentages_inputs[x, ]$start, df_base_percentages_inputs[x, ]$stop)
#
#     df_temp_reference_bases <- cbind(df_temp_reference_bases, df_base_percentages_inputs[x, ])
#
#     return(df_temp_reference_bases)
#
#   }))
#
#   df_base_percentages <- merge(
#     df_base_percentages,
#     df_all_reference_bases %>%
#       dplyr::select(-c(edit_site_strand, start, stop)),
#     by.x = c('specimen', 'annotation', 'edit_site_target.seq', 'edit.site', 'position'),
#     by.y = c('specimen', 'annotation', 'edit_site_target.seq', 'edit.site', 'reference_position')
#   )

  # return(df_base_percentages)

}

