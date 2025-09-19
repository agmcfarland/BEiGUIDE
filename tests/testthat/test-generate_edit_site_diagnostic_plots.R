test_that("test_diagnostic_plots", {


  for (edit_site_id in c('1', '3', '5', '12')) {

    edit_site <- readRDS(testthat::test_path('testdata', 'example_output_1', 'edit_sites', paste0('edit_site_id_', edit_site_id)))

    plot_results <- generate_edit_site_diagnostic_plots(edit_site = edit_site)

    testthat::expect_equal(length(names(plot_results)), 4)

  }


})
