test_that("standardize_cut_site_direction works for both strands", {
  genomic_positions <- 100:105
  bases <- c("A", "T", "G", "C", "C", "A")
  relative_zero <- 103

  # Positive strand
  result_plus <- standardize_cut_site_direction(
    strand = "+",
    relative_zero = relative_zero,
    genomic_positions = genomic_positions,
    bases = bases
  )

  testthat::expect_equal(nrow(result_plus), length(genomic_positions))
  testthat::expect_named(result_plus, c("genomic_position", "base", "relative_position"))

  # relative_position should be 103 - position
  testthat::expect_equal(as.character(result_plus$relative_position), as.character(relative_zero - genomic_positions))
  testthat::expect_equal(result_plus$base, bases)

  # Negative strand
  result_minus <- standardize_cut_site_direction(
    strand = "-",
    relative_zero = relative_zero,
    genomic_positions = genomic_positions,
    bases = bases
  )

  testthat::expect_equal(nrow(result_minus), length(genomic_positions))
  testthat::expect_named(result_minus, c("genomic_position", "base", "relative_position"))

  # relative_position should be position - 103
  testthat::expect_equal(as.character(result_minus$relative_position), as.character(genomic_positions - relative_zero))

  # bases should be reverse-complemented
  revcomp <- as.character(Biostrings::complement(
    Biostrings::DNAString(paste(bases, collapse = ""))
  ))
  revcomp <- strsplit(revcomp, "")[[1]]
  testthat::expect_equal(result_minus$base, revcomp)
})
