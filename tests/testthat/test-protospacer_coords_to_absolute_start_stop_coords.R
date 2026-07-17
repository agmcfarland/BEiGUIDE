test_that("protospacer_coords_to_absolute_start_stop_coords maps boundaries correctly", {

  # --- Test Case 1: Plus Strand (+) ---
  # Protospacer is left (starts at 84), PAM is right (ends at 106)
  coords_plus <- list(
    protospacer_end = 84,
    PAM_start = 106
  )

  res_plus <- protospacer_coords_to_absolute_start_stop_coords(
    protospacer_coords = coords_plus,
    strand = "+"
  )

  expect_equal(res_plus$start, 84)
  expect_equal(res_plus$stop, 106)


  # --- Test Case 2: Minus Strand (-) ---
  # PAM is left (ends at 80), Protospacer is right (starts at 116)
  coords_minus <- list(
    protospacer_end = 116,
    PAM_start = 80
  )

  res_minus <- protospacer_coords_to_absolute_start_stop_coords(
    protospacer_coords = coords_minus,
    strand = "-"
  )

  expect_equal(res_minus$start, 80)
  expect_equal(res_minus$stop, 116)


  # --- Test Case 3: Error Handling ---
  expect_error(
    protospacer_coords_to_absolute_start_stop_coords(coords_plus, strand = "invalid"),
    "Strand must be either"
  )
})


test_that("Integration: Coordinates pipe seamlessly from cut site to absolute start/stop", {

  # --- Plus Strand Integration ---
  # Cut site at 100, 20 bp spacer (17 distal, 3 proximal), 3 bp PAM
  plus_coords <- cut_site_to_protospacer_coordinates(
    cut_site = 100,
    strand = "+",
    protospacer_distance_from_cut_PAM_distal = 16,
    protospacer_distance_from_cut_PAM_proximal = 3,
    pam_length = 3
  )

  plus_absolute <- protospacer_coords_to_absolute_start_stop_coords(
    protospacer_coords = plus_coords,
    strand = "+"
  )

  expect_equal(plus_absolute$start, 85)  # 100 - 16
  expect_equal(plus_absolute$stop, 106)  # 100 + 3 + 3


  # --- Minus Strand Integration ---
  # Cut site at 100, 20 bp spacer (17 distal, 3 proximal), 3 bp PAM
  minus_coords <- cut_site_to_protospacer_coordinates(
    cut_site = 100,
    strand = "-",
    protospacer_distance_from_cut_PAM_distal = 16,
    protospacer_distance_from_cut_PAM_proximal = 3,
    pam_length = 3
  )

  minus_absolute <- protospacer_coords_to_absolute_start_stop_coords(
    protospacer_coords = minus_coords,
    strand = "-"
  )

  expect_equal(minus_absolute$start, 94)   # 100 - 3 - 3
  expect_equal(minus_absolute$stop, 115)  # 100 + 16
})
