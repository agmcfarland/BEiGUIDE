% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parallel_bam_file_to_table.R
\name{parallel_bam_file_list_to_table}
\alias{parallel_bam_file_list_to_table}
\title{Parallel BAM File List to Table}
\usage{
parallel_bam_file_list_to_table(
  bam_file_list,
  scan_param_what_list = c("qname", "rname", "strand", "pos", "qwidth", "seq", "cigar",
    "flag"),
  number_of_cpu = 4
)
}
\arguments{
\item{bam_file_list}{A character vector specifying the paths to BAM files.}

\item{scan_param_what_list}{A character vector specifying the scan parameters.}

\item{number_of_cpu}{An integer specifying the number of CPU cores to use for parallel processing.}
}
\value{
A data frame containing information from the BAM files.
}
\description{
Take in a list of BAM files, read in with scanBam, and tag each unique BAM file with the basename.
}
