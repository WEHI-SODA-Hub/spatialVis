library(dplyr)

# Heatmap marker plot
create_raster_plot <- function(long_mean_intensities, markers, celltypes,
                               parenttypes, celltype_colname, parent_colname,
                               fill = "standardised_mean") {

  # order factors according to arguments
  lmi_parent_factors <- factor(long_mean_intensities %>% pull(parent_colname),
                               levels = parenttypes)
  lmi_cell_factors <- factor(long_mean_intensities %>% pull(celltype_colname),
                             levels = rev(celltypes))
  lmi_marker_factors <- factor(long_mean_intensities %>% pull("marker"),
                               levels = markers)
  long_mean_intensities[, parent_colname] <- lmi_parent_factors
  long_mean_intensities[, celltype_colname] <- lmi_cell_factors
  long_mean_intensities[, "marker"] <- lmi_marker_factors

  raster_plot <- ggplot2::ggplot(long_mean_intensities,
                                 ggplot2::aes(x = marker, # nolint: object_usage_linter, line_length_linter.
                                              y = !!as.name(celltype_colname),
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

  return(raster_plot)
}

# Bar plot of cell type counts
create_bar_plot <- function(mean_intensities, celltypes, parenttypes,
                            celltype_colname, parent_colname) {
  # order factors according to arguments
  mi_parent_factors <- factor(mean_intensities %>% pull(parent_colname),
                              levels = parenttypes)
  mi_cell_factors <- factor(mean_intensities %>% pull(celltype_colname),
                            levels = rev(celltypes))
  mean_intensities[, parent_colname] <- mi_parent_factors
  mean_intensities[, celltype_colname] <- mi_cell_factors

  bar_plot <- ggplot2::ggplot(mean_intensities,
                              ggplot2::aes(x = !!as.name(celltype_colname),
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
    ggplot2::facet_grid(rows = vars(!!as.name(parent_colname)),
                        scales = "free_y", space = "free")

  return(bar_plot)
}

# calculate proportions of positive markers
get_proportions <- function(spe, markers, celltypes, celltype_colname,
                            parent_colname) {
  proportions <- SingleCellExperiment::colData(spe) %>%
    as.data.frame() %>%
    dplyr::mutate(dplyr::across(dplyr::all_of(markers),
                                ~ as.numeric(grepl("\\+", .)))) %>%
    dplyr::select(dplyr::all_of(c(celltype_colname, parent_colname,
                                  markers))) %>%
    dplyr::group_by(pick(!!as.name(celltype_colname),
                         !!as.name(parent_colname))) %>%
    dplyr::summarise(
      count = dplyr::n(),
      dplyr::across(dplyr::where(is.numeric), ~ sum(.) / length(.))
    ) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(!!as.name(parent_colname), count)  %>%
    dplyr::filter(!!as.name(celltype_colname) %in% celltypes) %>%
    as.data.frame()

  return(proportions)
}

# calculate mean intensities of markers
get_mean_intensities <- function(spe, celltypes, celltype_colname,
                                 parent_colname) {
  mean_intensities <- SingleCellExperiment::colData(spe) %>%
    as.data.frame() %>%
    dplyr::select(dplyr::all_of(c(celltype_colname, parent_colname))) %>%
    cbind(t(SingleCellExperiment::counts(spe))) %>%
    dplyr::group_by(pick(!!as.name(celltype_colname),
                         !!as.name(parent_colname))) %>%
    dplyr::summarise(
      count = dplyr::n(),
      dplyr::across(dplyr::where(is.numeric), mean)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(!!as.name(parent_colname), count)  %>%
    dplyr::filter(!!as.name(celltype_colname) %in% celltypes) %>%
    as.data.frame()

  return(mean_intensities)
}

#' @title Plot marker heatmap
#'
#' @description Plot a heatmap of mean marker intensities or proportions
#' @param spe SingleCellExperiment object containing marker intensities
#' and marker positivity
#' @param markers Character vector of markers to plot and the order in which
#' they should be plotted
#' @param celltypes Character vector of cell types to plot and the order in
#' which they should be plotted
#' @param parenttypes Character vector of parent cell types to plot and the
#' order in which they should be plotted
#' @param celltype_colname Column name in colData containing cell type
#' information (default: "HierarchyLevel4")
#' @param parent_colname Column name in colData containing parent cell type
#' information (default: "HierarchyLevel2")
#' @param value The value to plot. Either "expression" (default) or
#' "proportion"
#' @return A ggplot2 object
#' @export
plot_marker_heatmap <- function(spe, markers = NULL,
                                celltypes = NULL,
                                parenttypes = NULL,
                                celltype_colname = "HierarchyLevel4",
                                parent_colname = "HierarchyLevel2",
                                value = "expression") {

  stopifnot(celltype_colname %in% colnames(SingleCellExperiment::colData(spe)))
  stopifnot(parent_colname %in% colnames(SingleCellExperiment::colData(spe)))

  # get all cell types if not provided
  if (is.null(celltypes)) {
    celltypes <- unique(SingleCellExperiment::colData(spe)[[celltype_colname]])
  }
  if (is.null(parenttypes)) {
    parenttypes <- unique(SingleCellExperiment::colData(spe)[[parent_colname]])
  }
  if (is.null(markers)) {
    markers <- rownames(spe)
  }

  df <- switch(value,
               expression = get_mean_intensities(spe, celltypes,
                                                 celltype_colname,
                                                 parent_colname),
               proportion = get_proportions(spe, markers, celltypes,
                                            celltype_colname, parent_colname))

  stopifnot(nrow(df) > 0)

  # convert to long format and scale for plotting
  long_df <- df %>%
    tidyr::pivot_longer(-c(parent_colname, celltype_colname, "count"),
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

  raster_plot <- create_raster_plot(long_df, markers, celltypes,
                                    parenttypes, celltype_colname,
                                    parent_colname, fill = value)
  bar_plot <- create_bar_plot(df, celltypes, parenttypes,
                              celltype_colname, parent_colname)

  aligned <- cowplot::align_plots(raster_plot, bar_plot, align = "h",
                                  axis = "tblr")

  p <- gridExtra::grid.arrange(aligned[[1]], aligned[[2]],
                               ncol = 2, nrow = 1, widths = c(4, 1))

  return(p)
}
