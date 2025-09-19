test_that("cut_site_to_protospacer_coordinates works for + strand", {
  result <- cut_site_to_protospacer_coordinates(
    cut_site = 100,
    strand = '+',
    end_distance_from_cut_site = 20,
    expected_cut_distance_from_pam = 3,
    pam_length = 3
  )

  expect_type(result, "list")
  expect_named(result, c("PAM_start", "protospacer_start", "protospacer_end"))

  # + strand logic
  expect_equal(result$protospacer_start, 103)
  expect_equal(result$PAM_start, 106)
  expect_equal(result$protospacer_end, 80)
})

test_that("cut_site_to_protospacer_coordinates works for - strand", {
  result <- cut_site_to_protospacer_coordinates(
    cut_site = 100,
    strand = '-',
    end_distance_from_cut_site = 20,
    expected_cut_distance_from_pam = 3,
    pam_length = 3
  )

  expect_type(result, "list")
  expect_named(result, c("PAM_start", "protospacer_start", "protospacer_end"))

  # - strand logic
  expect_equal(result$protospacer_start, 97)
  expect_equal(result$PAM_start, 94)
  expect_equal(result$protospacer_end, 120)
})

test_that("cut_site is coerced to numeric", {
  result <- cut_site_to_protospacer_coordinates("100", "+")
  expect_equal(result$protospacer_start, 103)
})
