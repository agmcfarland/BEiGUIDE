#' Generate diagnostic plots for an edit site
#'
#' This function takes a single `edit_site` object and produces a set of four
#' diagnostic plots to visualize base composition, editing rates, read depth,
#' and per-base composition across the targeted locus. It is useful for
#' visually assessing editing outcomes and coverage at a specific target site.
#'
#' @param edit_site A named list containing:
#' \describe{
#'   \item{description}{A data frame or list with annotation information about the site, including at least:
#'     \code{annotation}, \code{edit.site}, \code{gene_id}, \code{target.seq}, \code{abund}, and \code{position}}
#'   \item{standardized_locus}{A data frame describing the reference locus around the edit site, including
#'     \code{genomic_position}, \code{base}, and \code{relative_position}}
#'   \item{base_composition}{A data frame describing observed base composition at each genomic position,
#'     including \code{position}, \code{base}, \code{percentage}, and \code{position_depth}}
#'   \item{significant_edits}{A list containing a data frame \code{significant} with columns
#'     \code{position}, \code{base}, \code{p_value}, and \code{on_target_edit}}
#' }
#'
#' @return A named list of four \code{ggplot} objects:
#' \describe{
#'   \item{heatmap}{A heatmap showing base composition percentages per position and base, with significance annotations and expected/on/off-target labels}
#'   \item{on_target_only}{A stacked bar plot showing the percentage of reads that are on-target edits at each position}
#'   \item{depth}{A bar plot showing read depth at each genomic position}
#'   \item{stacked_percentage}{A stacked bar plot showing the percentage composition of each base at each position}
#' }
#'
#' @details
#' The function generates ggplots using \pkg{ggplot2} and custom color scales.
#' The first heatmap highlights significant edits and the expected cut site,
#' while the other plots summarize on-target editing frequency, sequencing depth,
#' and overall base composition.
#'
#' @examples
#' \dontrun{
#' plots <- generate_edit_site_diagnostic_plots(my_edit_site)
#' plots$heatmap
#' plots$depth
#' }
#'
#' @import ggplot2 dplyr
#'
#' @export
generate_edit_site_diagnostic_plots <- function(edit_site) {

  plot_title <- as.character(edit_site$description$annotation)
  plot_subtitle <- paste0(edit_site$description$edit.site, ' ', edit_site$description$gene_id, ' sgRNA: ', edit_site$description$target.seq, ' abund: ', edit_site$description$abund)

  base_order <- c('A', 'G', 'T', 'C', 'N')
  base_plot_colors <- c("#4285F4", "#34A853", "#FBBC05", "#EA4335", 'lightgrey')
  names(base_plot_colors) <- base_order

  df_standardized_locus <- edit_site$standardized_locus %>%
    dplyr::arrange(relative_position) %>%
    dplyr::mutate(genomic_position = factor(genomic_position, levels = genomic_position))

  df_base_composition <- edit_site$base_composition %>%
    dplyr::left_join(
      edit_site$significant_edits$significant %>%
        dplyr::select(position, base, p_value, on_target_edit),
      by = c('position', 'base')
    ) %>%
    dplyr::mutate(
      genomic_position = factor(position, levels = df_standardized_locus$genomic_position),
      percentage_fill = ifelse(percentage > 5, paste(round(percentage, 1)), NA_real_),
      signifiance = dplyr::case_when(
        p_value < 0.001 ~ '***',
        p_value < 0.01 ~ '**',
        p_value < 0.05 ~ '*',
        TRUE ~ NA_character_
      ),
      percentage_significance_fill = paste0(percentage_fill, '\n', signifiance),
      percentage_significance_fill = stringr::str_replace_all(percentage_significance_fill, 'NA', ''),
      expected_edit = dplyr::case_when(
        (on_target_edit == TRUE) ~ 'on_target_edit',
        (on_target_edit == FALSE) ~ 'off_target_edit',
        base == reference_base ~ 'expected_base',
        TRUE ~ 'other'
      ),
      base = factor(base, levels = base_order)
    )

  cut_site_position <- factor(edit_site$description$position, levels = df_standardized_locus$genomic_position)

  color_plot_values <- c('#984EA3', 'orange', 'green', 'transparent')
  names(color_plot_values) <- c('expected_base', 'on_target_edit', 'off_target_edit', 'other')
  color_plot_values <- color_plot_values[names(color_plot_values) %in% base::unique(df_base_composition$expected_edit)]

  x_axis_plot_name <- df_standardized_locus$genomic_position
  names(x_axis_plot_name) <- paste0(df_standardized_locus$base, '\n', df_standardized_locus$relative_position)

  # diagnostic plot 1
  p1 <- ggplot(df_standardized_locus, aes(x = genomic_position)) +
    geom_vline(xintercept = cut_site_position, color = 'grey', linetype = 'dashed') +
    geom_tile(
      data = df_base_composition,
      aes(x = genomic_position, y = base, fill = percentage, color = expected_edit),
      linewidth = 0.8,
      width = 0.8,
      height = 0.8
    ) +
    geom_text(
      data = df_base_composition,
      aes(label = percentage_significance_fill, y = base, x = genomic_position),
      size = 2,
      color = 'white',
    ) +
    theme_light() +
    theme(
      aspect.ratio = 0.25,
      axis.title = element_blank(),
      legend.position = 'bottom',
      legend.title.position = 'top',
      legend.title = element_text(face = 'bold', hjust = 0.5)
    ) +
    scale_y_discrete(labels = base_order, drop = F) +
    scale_x_discrete(labels = names(x_axis_plot_name), drop = F) +
    scale_fill_gradientn(
      colors = c("transparent", "lightblue", "darkblue"),
      values = scales::rescale(c(0, 0.01, max(df_base_composition$percentage))),
      na.value = "white",
      limits = c(0, 100)
    ) +
    scale_color_manual(values = color_plot_values, drop = F) +
    guides(
      color = guide_legend(
        override.aes = list(
          color = color_plot_values,
          fill = 'white'
        )
      ),
      fill = guide_colorbar(
        barwidth = unit(5, 'cm'),
        barheight = unit(0.5, 'cm')
      )
    ) +
    labs(title = plot_title, subtitle = plot_subtitle)


  # diagnostic plot 2
  fill_plot_values_2 <- c('orange', 'lightblue', 'transparent')
  names(fill_plot_values_2) <- c(TRUE, FALSE)

  if (df_base_composition %>%
      dplyr::filter(!is.na(on_target_edit)) %>% nrow() > 0 ) {

    p2 <- ggplot(df_standardized_locus, aes(x = genomic_position)) +
      geom_bar(
        data = df_base_composition %>%
          dplyr::filter(!is.na(on_target_edit)),
        stat = 'identity',
        position = position_stack(),
        aes(x = genomic_position, y = percentage, fill = on_target_edit)
      ) +
      scale_x_discrete(labels = names(x_axis_plot_name), drop = F) +
      theme_light() +
      theme(
        aspect.ratio = 0.25,
        axis.title.x = element_blank(),
        legend.position = 'bottom',
        legend.title.position = 'top',
        legend.title = element_text(face = 'bold', hjust = 0.5)
      ) +
      labs(y = 'percentage') +
      scale_fill_manual(values = fill_plot_values_2) +
      labs(title = plot_title, subtitle = plot_subtitle)
  } else {
    p2 <- ggplot(df_standardized_locus, aes(x = genomic_position)) +
      geom_blank() +
      scale_x_discrete(labels = names(x_axis_plot_name), drop = F) +
      theme_light() +
      theme(
        aspect.ratio = 0.25,
        axis.title.x = element_blank(),
        legend.position = 'bottom',
        legend.title.position = 'top',
        legend.title = element_text(face = 'bold', hjust = 0.5)
      ) +
      labs(y = 'percentage') +
      labs(title = plot_title, subtitle = plot_subtitle)
    }

  # diagnostic plot 3
  p3 <- ggplot(df_standardized_locus, aes(x = genomic_position)) +
    geom_bar(
      data = df_base_composition %>%
        dplyr::select(genomic_position, position_depth) %>%
        dplyr::distinct(),
      stat = 'identity',
      aes(x = genomic_position, y = position_depth)
    ) +
    scale_x_discrete(labels = names(x_axis_plot_name), drop = F) +
    theme_light() +
    theme(
      aspect.ratio = 0.25,
      axis.title.x = element_blank(),
      legend.position = 'bottom',
      legend.title.position = 'top',
      legend.title = element_text(face = 'bold', hjust = 0.5)
    ) +
    labs(title = plot_title, subtitle = plot_subtitle, y = 'depth')

  # diagnostic plot 4
  p4 <- ggplot(df_standardized_locus, aes(x = genomic_position)) +
    geom_bar(
      data = df_base_composition,
      stat = 'identity',
      position = position_stack(),
      aes(x = genomic_position, y = percentage, fill = base)
    ) +
    scale_fill_manual(values = base_plot_colors) +
    scale_x_discrete(labels = names(x_axis_plot_name), drop = F) +
    theme_light() +
    theme(
      aspect.ratio = 0.25,
      axis.title.x = element_blank(),
      legend.position = 'bottom',
      legend.title.position = 'top',
      legend.title = element_text(face = 'bold', hjust = 0.5)
    ) +
    labs(title = plot_title, subtitle = plot_subtitle, y = 'percentage')

  return(
    list(
      'heatmap' = p1,
      'on_target_only' = p2,
      'depth' = p3,
      'stacked_percentage' = p4
    )
  )
}

