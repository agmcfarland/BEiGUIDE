test_that("quantify_edits works", {
  testthat::skip()
  test_data <- testthat::test_path('testdata', 'integration_1')

  test_output <- file.path(testthat::test_path('integration_test_1'))

  dir.create(test_output, showWarnings = FALSE)

  quantify_edits(
    base_directory = test_data,
    output_directory = test_output,
    n_processors = 8
  )

  df_edit_sites <- readRDS(file.path(test_output, 'quantify_edits', 'edit_sites.rds'))

  df_base_percentages <- readRDS(file.path(test_output, 'quantify_edits', 'base_percentages.rds'))

  testthat::expect_equal(nrow(df_edit_sites), 11)

  testthat::expect_equal(ncol(df_edit_sites), 13)

  testthat::expect_equal(nrow(df_base_percentages), 147)

  testthat::expect_equal(ncol(df_base_percentages), 16)

})
