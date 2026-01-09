test_that("aln_to_base_and_position_table works", {

  df_bam_filtered_expected <- readRDS(testthat::test_path('testdata', 'df_bam_filtered_example.rds'))

  df_aln_result <- aln_to_base_and_position_table(df_bam_filtered = df_bam_filtered_expected)

  testthat::expect_equal(dplyr::n_distinct(df_bam_filtered_expected$qname), dplyr::n_distinct(df_aln_result$qname))

  # Test logic of lapply loop that expands each row in df_bam_filtered
  qname_ <- 'MN01490:144:000H5KVFN:1:23104:18170:3323'
  df_subset_input <- dplyr::filter(df_bam_filtered_expected, qname == qname_)
  df_subset_result <- dplyr::filter(df_aln_result, qname == qname_)
  testthat::expect_equal(min(df_subset_result$position), df_subset_input$pos)
  testthat::expect_equal(max(df_subset_result$position), df_subset_input$pos + df_subset_input$qwidth)
  testthat::expect_equal(max(df_subset_result$position), df_subset_input$pos + df_subset_input$qwidth)

  # Test logic of rbinded lapply loop output that is the fully expanded table from the original df_bam_filtered
  testthat::expect_equal(nrow(df_aln_result), sum(df_bam_filtered_expected$qwidth))


})
