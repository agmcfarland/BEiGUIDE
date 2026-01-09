#' BAM File to Table
#'
#' Take in a BAM file, read in with scanBam, and tag each unique BAM file with the filepath's basename.
#'
#' @param bam_file A character string specifying the path to the BAM file.
#' @param scan_param_what_list A character vector specifying the scan parameters.
#'
#' @return A data frame containing information from the BAM file.
#'
#' @import Rsamtools dplyr
#'
#' @export
bam_file_to_table <- function(bam_file, scan_param_what_list = c('qname', 'rname', 'strand', 'pos', 'qwidth', 'seq', 'cigar', 'flag')) {

  scan_param <- Rsamtools::ScanBamParam(what = scan_param_what_list)

  df_bam <- Rsamtools::BamFile(bam_file) %>%
    Rsamtools::scanBam(param = scan_param) %>%
    base::as.data.frame() %>%
    dplyr::mutate(bam_id = base::basename(bam_file))

  return(df_bam)
}
