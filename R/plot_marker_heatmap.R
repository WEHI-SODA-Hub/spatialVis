library(dplyr)

#' @title Plot marker heatmap
#'
#' @description Plot a heatmap of mean marker intensities or proportions
#' @param spe SingleCellExperiment object containing marker intensities
#' and marker positivity
#' @param markers Character vector of markers to plot and the order in which
#' they should be plotted (NOTE: for expression, check that the markers are
#' present in the counts, and for proportions, check that they are present in
#' the colData)
#' @param cell_types Character vector of cell types to plot and the order in
#' which they should be plotted
#' @param parent_types Character vector of parent cell types to plot and the
#' order in which they should be plotted
#' @param cell_type_colname Column name in colData containing cell type
#' information (default: "HierarchyLevel4")
#' @param parent_colname Column name in colData containing parent cell type
#' information (default: "HierarchyLevel2")
#' @param value The value to plot. Either "expression" (default) or
#' "proportion"
#' @return A ggplot2 object
#' @export
#' @importFrom dplyr %>%
plot_marker_heatmap <- function(spe, markers = NULL,
                                cell_types = NULL,
                                parent_types = NULL,
                                cell_type_colname = "HierarchyLevel4",
                                parent_colname = "HierarchyLevel2",
                                value = "expression") {
  stopifnot(value %in% c("expression", "proportion"))
  stopifnot(cell_type_colname %in% colnames(SingleCellExperiment::colData(spe)))
  stopifnot(parent_colname %in% colnames(SingleCellExperiment::colData(spe)))

  # get all cell types if not provided
  if (is.null(cell_types)) {
    cell_types <- unique(SingleCellExperiment::colData(spe)[[cell_type_colname]]) # nolint: line_length_linter.
  }
  if (is.null(parent_types)) {
    parent_types <- unique(SingleCellExperiment::colData(spe)[[parent_colname]])
  }
  if (is.null(markers)) {
    markers <- switch(value,
                      expression = rownames(spe),
                      proportion = get_marker_names(spe))
  }

  df <- switch(value,
               expression = get_mean_intensities(spe, cell_types,
                                                 cell_type_colname,
                                                 parent_colname),
               proportion = get_proportions(spe, markers, cell_types,
                                            cell_type_colname, parent_types,
                                            parent_colname))

  stopifnot(nrow(df) > 0)

  # convert to long format and scale for plotting
  long_df <- df %>%
    tidyr::pivot_longer(-dplyr::all_of(c(parent_colname,
                                         cell_type_colname, "count")),
                        names_to = "marker",
                        values_to = value) %>%
    dplyr::group_by(marker) # nolint: object_usage_linter.

  if (value == "expression") {
    long_df <- long_df %>%
      dplyr::mutate(standardised_mean = scale(expression)[, 1]) %>%
      dplyr::filter(marker %in% markers) # nolint: object_usage_linter.
    value <- "standardised_mean"
  } else {
    long_df <- long_df %>%
      dplyr::filter(marker %in% markers) # nolint: object_usage_linter
  }

  raster_plot <- create_raster_plot(long_df, markers, cell_types,
                                    parent_types, cell_type_colname,
                                    parent_colname, fill = value)
  bar_plot <- create_bar_plot(df, cell_types, parent_types,
                              cell_type_colname, parent_colname)

  aligned <- cowplot::align_plots(raster_plot, bar_plot, align = "h",
                                  axis = "tblr")

  p <- gridExtra::grid.arrange(aligned[[1]], aligned[[2]],
                               ncol = 2, nrow = 1, widths = c(4, 1))

  p
}

# Heatmap marker plot
create_raster_plot <- function(long_mean_intensities, markers, cell_types,
                               parent_types, cell_type_colname, parent_colname,
                               fill = "standardised_mean") {

  # order factors according to arguments
  lmi_parent_factors <- factor(long_mean_intensities %>%
                                 dplyr::pull(parent_colname),
                               levels = parent_types)
  lmi_cell_factors <- factor(long_mean_intensities %>%
                               dplyr::pull(cell_type_colname),
                             levels = rev(cell_types))
  lmi_marker_factors <- factor(long_mean_intensities %>% dplyr::pull("marker"),
                               levels = markers)
  long_mean_intensities[, parent_colname] <- lmi_parent_factors
  long_mean_intensities[, cell_type_colname] <- lmi_cell_factors
  long_mean_intensities[, "marker"] <- lmi_marker_factors

  raster_plot <- ggplot2::ggplot(long_mean_intensities,
                                 ggplot2::aes(x = marker, # nolint: object_usage_linter, line_length_linter.
                                              y = !!as.name(cell_type_colname),
                                              fill = !!as.name(fill))) +
    ggplot2::geom_raster() +
    ggplot2::geom_tile(ggplot2::aes(fill = !!as.name(fill)),
                       colour = "black", linewidth = 0.5) +
    ggplot2::scale_fill_distiller(palette = "YlGnBu") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::facet_grid(rows = dplyr::vars(!!as.name(parent_colname)),
                        scales = "free_y", space = "free") +
    ggplot2::theme(
      panel.background = ggplot2::element_blank(),
      strip.text.y = ggplot2::element_blank(),
      legend.position = "top",
      plot.margin = ggplot2::margin(t = 0, # Top margin
                                    r = -0.8, # Right margin
                                    b = 0, # Bottom margin
                                    l = 0),
      panel.spacing = ggplot2::unit(0.001, "cm")
    )

  raster_plot
}

# best attempt to get positive/negative marker names from metadata
get_marker_names <- function(spe) {
  markers <- SingleCellExperiment::colData(spe) %>%
    as.data.frame() %>%
    dplyr::select(-dplyr::starts_with("Hierarchy"),
                  -dplyr::starts_with("cluster"),
                  -sample_id, -Cell_ID, -In.Tumour, -Image) # nolint: object_usage_linter, line_length_linter.

  colnames(markers)
}

# Bar plot of cell type counts
create_bar_plot <- function(mean_intensities, cell_types, parent_types,
                            cell_type_colname, parent_colname) {
  # order factors according to arguments
  mi_parent_factors <- factor(mean_intensities %>% dplyr::pull(parent_colname),
                              levels = parent_types)
  mi_cell_factors <- factor(mean_intensities %>% dplyr::pull(cell_type_colname),
                            levels = rev(cell_types))
  mean_intensities[, parent_colname] <- mi_parent_factors
  mean_intensities[, cell_type_colname] <- mi_cell_factors

  bar_plot <- ggplot2::ggplot(mean_intensities,
                              ggplot2::aes(x = !!as.name(cell_type_colname),
                                           y = count)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::theme(axis.title.y = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(angle = 90, hjust = 1),
                   axis.title.x = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank(),
                   strip.text.y = ggplot2::element_blank(),
                   plot.margin = ggplot2::margin(t = 0,  # Top margin
                                                 r = 10,  # Right margin
                                                 b = 0,  # Bottom margin
                                                 l = -0.8),
                   panel.background = ggplot2::element_blank(),
                   panel.spacing = ggplot2::unit(0.001, "cm")) +
    ggplot2::coord_flip(ylim = c(0, max(mean_intensities$count))) +
    ggplot2::facet_grid(rows = dplyr::vars(!!as.name(parent_colname)),
                        scales = "free_y", space = "free")

  bar_plot
}

# calculate proportions of positive markers
get_proportions <- function(spe, markers, cell_types, cell_type_colname,
                            parent_types, parent_colname) {
  proportions <- SingleCellExperiment::colData(spe) %>%
    as.data.frame() %>%
    dplyr::filter(!!as.name(parent_colname) %in% parent_types) %>% # nolint: object_usage_linter, line_length_linter.
    dplyr::filter(!!as.name(cell_type_colname) %in% cell_types) %>%
    dplyr::mutate(dplyr::across(dplyr::all_of(markers),
                                ~ as.numeric(grepl("\\+", .)))) %>%
    dplyr::select(dplyr::all_of(c(cell_type_colname, parent_colname,
                                  markers))) %>%
    dplyr::group_by(pick(!!as.name(cell_type_colname),
                         !!as.name(parent_colname))) %>%
    dplyr::summarise(
      count = dplyr::n(),
      dplyr::across(dplyr::where(is.numeric), ~ sum(.) / length(.))
    ) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(!!as.name(parent_colname), count)  %>%
    as.data.frame()

  proportions
}

# calculate mean intensities of markers
get_mean_intensities <- function(spe, cell_types, cell_type_colname,
                                 parent_colname) {
  mean_intensities <- SingleCellExperiment::colData(spe) %>%
    as.data.frame() %>%
    dplyr::select(dplyr::all_of(c(cell_type_colname, parent_colname))) %>%
    cbind(t(SingleCellExperiment::counts(spe))) %>%
    dplyr::group_by(pick(!!as.name(cell_type_colname),
                         !!as.name(parent_colname))) %>%
    dplyr::summarise(
      count = dplyr::n(),
      dplyr::across(dplyr::where(is.numeric), mean)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(!!as.name(parent_colname), count)  %>%
    dplyr::filter(!!as.name(cell_type_colname) %in% cell_types) %>%
    as.data.frame()

  mean_intensities
}
