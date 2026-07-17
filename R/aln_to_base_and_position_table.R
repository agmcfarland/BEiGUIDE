#' Alignment to Base and Position Table
#'
#' Takes row produced by filter_bam_alignments and returns a table with base, position, and qname.
#'
#' @details input is a single row produced by filter_bam_alignments
#' otherwise input must have seq, pos, qwidth, and qname columns
#'
#' @param aln A single row data frame.
#'
#' @return A data frame with base, position, and qname.
#'
#' @import data.table
#'
#' @export
# aln_to_base_and_position_table <- function(aln) {
#   return(data.frame(
#     'base' = stringr::str_split(aln$seq, '')[[1]],
#     'position' = base::seq(aln$pos, aln$pos + aln$qwidth - 1, by = 1), # stays the same for + or -
#     'qname' = aln$qname
#   ))
# }

aln_to_base_and_position_table <- function(df_bam_filtered) {
  # Make data.table object
  df_bam_filtered <- data.table::as.data.table(df_bam_filtered)

  results <- base::lapply(1:nrow(df_bam_filtered), function(idx) {

    # Extract the single row
    df_idx_chunk <- df_bam_filtered[idx, ]

    seq_char <- base::as.character(df_idx_chunk$seq)
    base_lists <- base::strsplit(seq_char, "", fixed = TRUE)[[1]]

    # Compute length for specific read
    len <- df_idx_chunk$qwidth

    data.table::data.table(
      base = base_lists,
      position = df_idx_chunk$pos + 0:(len - 1),
      qname = df_idx_chunk$qname
    )
  })

  return(data.table::rbindlist(results))
}

