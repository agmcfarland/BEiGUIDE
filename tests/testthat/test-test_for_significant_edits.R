testthat::test_that("test_for_significant_edits correctly filters and flags edits", {
  df <- data.frame(
    position = 1:4,
    base = c('G', 'A', 'T', 'C'),
    reference_base = c('A', 'A', 'T', 'C'),
    position_depth = c(20, 20, 20, 20),
    base_count = c(10, 2, 1, 1),
    p_value = c(0.01, 0.2, 0.5, 0.001),
    conf_low = c(0.2, 0.01, 0.0, 0.0),
    conf_high = c(0.8, 0.3, 0.1, 0.1),
    percentage = c(50, 10, 5, 5),
    stringsAsFactors = FALSE
  )

  # Positive strand test
  result_pos <- test_for_significant_edits(df)
  testthat::expect_true(all(c("all", "significant") %in% names(result_pos)))

  # All output should not include bases matching the reference
  testthat::expect_true(all(result_pos$all$base != result_pos$all$reference_base))

  # On-target flag is correct (A→G)
  testthat::expect_true(any(result_pos$all$on_target_edit))

  # Significant should only have p_value <= 0.05
  testthat::expect_true(all(result_pos$significant$p_value <= 0.05))

})
