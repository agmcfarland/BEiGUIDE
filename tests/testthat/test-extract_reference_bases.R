test_that("extract_reference_bases works", {

  set.seed(3)

  dna_letters <- data.frame(
    base = c('A', 'T','G', 'C', 'N'),
    probability_1 = c(0.24, 0.24, 0.24, 0.24, 0.04),
    probability_2 = c(0.18, 0.18, 0.3, 0.3, 0.04)
  )

  multifasta <- Biostrings::DNAStringSet(
    c(
      paste(sample(x = dna_letters$base, size = 10000, replace = TRUE, prob = dna_letters$probability_1), collapse = ''),
      paste(sample(x = dna_letters$base, size = 10000, replace = TRUE, prob = dna_letters$probability_2), collapse = '')
    ))

  names(multifasta) <- c('chr1', 'chr2')


  df_base_percentages <- data.frame(
    'position' = seq(100, 120),
    'specimen' = rep('GTSP1234', 21),
    'annotation' = rep('annot 1 (A)', 21),
    'base_count' = rep(100, 21),
    'position_depth' = rep(100, 21),
    'percentage' = rep(100, 21),
    'base' = 'A',
    'edit_site_target.seq' = rep('target_gene1', 21),
    'edit.site' = rep('chr1:+:120', 21),
    'edit_site_chromosome' = rep('chr1', 21),
    'edit_site_strand' = rep('+', 21)
  )

  df_base_percentages <- rbind(df_base_percentages, data.frame(
    'position' = seq(200, 220),
    'specimen' = rep('GTSP9999', 21),
    'annotation' = rep('annot 2 (A)', 21),
    'base_count' = rep(100, 21),
    'position_depth' = rep(100, 21),
    'percentage' = rep(100, 21),
    'base' = 'A',
    'edit_site_target.seq' = rep('target_gene2', 21),
    'edit.site' = rep('chr2:+:220', 21),
    'edit_site_chromosome' = rep('chr2', 21),
    'edit_site_strand' = rep('+', 21)
  ))

  df_base_percentages <- rbind(df_base_percentages, data.frame(
    'position' = seq(200, 220),
    'specimen' = rep('GTSP1234', 21),
    'annotation' = rep('annot 2 (B)', 21),
    'base_count' = rep(100, 21),
    'position_depth' = rep(100, 21),
    'percentage' = rep(100, 21),
    'base' = 'A',
    'edit_site_target.seq' = rep('target_gene2', 21),
    'edit.site' = rep('chr2:+:220', 21),
    'edit_site_chromosome' = rep('chr2', 21),
    'edit_site_strand' = rep('+', 21)
  ))


  reference_bases <- extract_reference_bases(multifasta, df_base_percentages)

  testthat::expect_equal(nrow(reference_bases), nrow(df_base_percentages))



  #### Testing with real data ###
#   devtools::load_all()
#
#   genomic_sequence <- Biostrings::readDNAStringSet('/data/iGUIDE/genomes/hg38.fasta')
#
  # df_base_percentages <- readRDS('/data/BEiGUIDE/tests/quantify_edits/base_percentages.rds') %>%
  #   dplyr::rename(edit_site_chromosome = edit_site_chromsome)
#
#   reference_bases <- extract_reference_bases(genomic_sequence, readRDS('/data/BEiGUIDE/tests/quantify_edits/base_percentages.rds') %>%
#                                                dplyr::rename(edit_site_chromosome = edit_site_chromsome))


})
