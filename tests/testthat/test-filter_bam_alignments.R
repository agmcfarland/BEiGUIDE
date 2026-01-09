test_that("filter_bam_alignments works positive strand", {

  df_bam <- readRDS(testthat::test_path('testdata', 'df_bam_example.rds'))

  testthat::expect_equal(nrow(df_bam), 200000)

  df_bam_filtered <- filter_bam_alignments(
    df_bam,
    specimen = 'GTSP5614',
    chromosome = 'chr1',
    edit_site_position = 198706754,
    edit_site_strand = '+',
    allowed_bases_within_gRNA = 3,
    allowed_bases_outside_gRNA = 3
  )

  df_expected_bam_filtered <- readRDS(testthat::test_path('testdata', 'df_bam_filtered_example.rds'))

  testthat::expect_equal(ncol(df_bam_filtered), ncol(df_expected_bam_filtered))

  testthat::expect_equal(nrow(df_bam_filtered), nrow(df_expected_bam_filtered))

  testthat::expect_equal(df_bam_filtered$seq, df_expected_bam_filtered$seq)

  testthat::expect_equal(base::unique(df_bam_filtered$edit_site_min < df_bam_filtered$edit_site_max), TRUE)

})

test_that("filter_bam_alignments works negative strand", {

  df_bam <- readRDS(testthat::test_path('testdata', 'df_bam_example.rds')) %>%
    dplyr::filter(
    !stringr::str_detect(cigar, 'I'),
    !stringr::str_detect(cigar, 'S'),
    !stringr::str_detect(cigar, 'D')
    )

  df_bam_filtered <- filter_bam_alignments(
    df_bam,
    specimen = 'GTSP5614',
    chromosome = 'chr2',
    edit_site_position = 143961596,
    edit_site_strand = '-',
    allowed_bases_within_gRNA = 3,
    allowed_bases_outside_gRNA = 3
  )

  testthat::expect_equal(base::unique(df_bam_filtered$edit_site_min < df_bam_filtered$edit_site_max), TRUE)

  testthat::expect_equal(ncol(df_bam_filtered), ncol(readRDS(testthat::test_path('testdata', 'df_bam_filtered_example.rds'))))

  testthat::expect_equal(nrow(df_bam_filtered), 8)

  testthat::expect_equal(base::unique(df_bam_filtered$edit_site_min < df_bam_filtered$edit_site_max), TRUE)

})

