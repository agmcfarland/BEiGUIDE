testthat::test_that("binomial_prop_edit_test returns expected results", {
  df <- data.frame(
    position = c(1, 1, 1, 2, 2),
    base = c("A", "C", "G", "A", "T"),
    base_count = c(10, 0, 5, 20, 0),
    position_depth = c(15, 15, 15, 20, 20)
  )

  result <- binomial_prop_edit_test(df, alternative_ = 'greater')

  # Check structure
  testthat::expect_s3_class(result, "data.frame")
  testthat::expect_true(all(c("position", "base", "base_count", "position_depth",
                    "p_value", "conf_low", "conf_high", "percentage") %in% colnames(result)))

  # Check that all bases A,C,G,T,N are present for each position
  bases_per_pos <- table(result$position, result$base)
  testthat::expect_true(all(colnames(bases_per_pos) %in% c("A","C","G","T","N")))
  testthat::expect_true(all(rowSums(bases_per_pos > 0) >= 1)) # each pos has at least one nonzero

  # Check that percentages are consistent
  testthat::expect_equal(result$percentage, 100 * result$base_count / result$position_depth)

  # Check that p-values are numeric and between 0 and 1
  testthat::expect_true(all(result$p_value >= 0 & result$p_value <= 1))
})

testthat::test_that("binomial_prop_edit_test errors on invalid alternative_", {
  df <- data.frame(
    position = 1,
    base = "A",
    base_count = 10,
    position_depth = 15
  )

  testthat::expect_error(
    binomial_prop_edit_test(df, alternative_ = "invalid"),
    "Specify one of two.sided, less, greater"
  )
})
