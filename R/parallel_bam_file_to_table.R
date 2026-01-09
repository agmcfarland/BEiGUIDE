#' Parallel BAM File List to Table
#'
#' Take in a list of BAM files, read in with scanBam, and tag each unique BAM file with the basename.
#'
#' @param bam_file_list A character vector specifying the paths to BAM files.
#' @param scan_param_what_list A character vector specifying the scan parameters.
#' @param number_of_cpu An integer specifying the number of CPU cores to use for parallel processing.
#'
#' @return A data frame containing information from the BAM files.
#'
#' @import parallel
#' 
#' @export
parallel_bam_file_list_to_table <- function(bam_file_list, 
                                           scan_param_what_list = c('qname', 'rname', 'strand', 'pos', 'qwidth', 'seq', 'cigar', 'flag'), 
                                           number_of_cpu = 4) {
  results_list <- parallel::mclapply(
    bam_file_list,
    function(bam_file) {
      bam_file_to_table(bam_file, scan_param_what_list)
    },
    mc.cores = number_of_cpu
  )

  df_bam <- do.call(rbind, results_list)


  df_bam$specimen_id <- stringr::str_extract(df_bam$bam_id, "^[^-]+")

  return(df_bam)
}