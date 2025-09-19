test_that("make_edit_site_table works", {

  subsetted_test_data_path <- testthat::test_path('testdata', 'integration_1')

  df_ft_data <- pull_ft_data_tables(base_directory = subsetted_test_data_path)

  df_annotations <- pull_combo_overview_table(base_directory = subsetted_test_data_path)

  df_edit_sites <- make_edit_site_table(
    ft_data_table = df_ft_data,
    spec_info_combo_overview_table = df_annotations,
    abundance_cutoff = 3)

  testthat::expect_equal(5815, nrow(df_edit_sites))

  testthat::expect_equal(11, nrow(df_edit_sites %>% dplyr::filter(pass_abundance_filter == T)))

  testthat::expect_equal(15, ncol(df_edit_sites))


})
