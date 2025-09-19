
test_that("manage_run_directory works correctly", {
  # Create a temporary directory for testing
  temp_dir <- tempdir()
  run_params <- list(
    output_directory = temp_dir,
    analysis_output = file.path(temp_dir, "analysis"),
    overwrite = FALSE
  )
  unlink(temp_dir, recursive = T)
  unlink(file.path(temp_dir, "analysis"), recursive = T)

  # Case 1: output_directory does not exist
  non_existent_dir <- file.path(temp_dir, "non_existent")
  run_params$output_directory <- non_existent_dir
  manage_run_directory(run_params)
  testthat::expect_true(dir.exists(run_params$output_directory))
  testthat::expect_true(dir.exists(run_params$analysis_output))

  # Case 2: analysis_output exists and overwrite is TRUE
  dir.create(run_params$analysis_output, showWarnings = FALSE)
  run_params$output_directory <- temp_dir
  run_params$overwrite <- TRUE
  testthat::expect_true(dir.exists(run_params$analysis_output))

  # Case 3: analysis_output exists and overwrite is FALSE
  dir.create(run_params$analysis_output, showWarnings = FALSE)
  run_params$overwrite <- FALSE
  testthat::expect_error(
    manage_run_directory(run_params),
    paste("analysis directory exists and overwrite is set to FALSE ", run_params$analysis_output)
  )

  # Case 4: analysis_output does not exist
  unlink(run_params$analysis_output, recursive = TRUE)
  manage_run_directory(run_params)
  testthat::expect_true(dir.exists(run_params$analysis_output))
})
