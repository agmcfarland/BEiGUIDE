test_that("calculate_base_percentages works", {

  df_aln_pos_filtered_expected <- readRDS(testthat::test_path('testdata', 'df_aln_pos_filtered_example.rds'))

  start <- 198706754 - 19
  stop <- 198706754 + 6

  df_base_percentages <- calculate_base_percentages(
    df_aln_pos = readRDS(testthat::test_path('testdata', 'df_aln_pos_example.rds')),
    start = start,
    stop = stop
  )

  testthat::expect_equal(nrow(df_base_percentages), nrow(df_aln_pos_filtered_expected))

  testthat::expect_equal(ncol(df_base_percentages), ncol(df_aln_pos_filtered_expected))


})

test_that("calculate_base_percentages errors", {

  df_aln_pos_filtered_expected <- readRDS(testthat::test_path('testdata', 'df_aln_pos_filtered_example.rds'))

  stop <- 198706754 - 19
  start <- 198706754 + 6

  testthat::expect_error(calculate_base_percentages(
    df_aln_pos = readRDS(testthat::test_path('testdata', 'df_aln_pos_example.rds')),
    start = start,
    stop = stop
  ))


})
