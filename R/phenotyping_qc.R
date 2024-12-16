library(dplyr)

#' @title Plot marker heatmap
#'
#' @description Plot a heatmap of mean marker intensities in cell types.
#' @param spe SingleCellExperiment object containing marker intensities
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
#' @return A ggplot2 object
#' @export
plot_marker_heatmap <- function(spe, markers, celltypes, parenttypes,
                                celltype_colname = "HierarchyLevel4",
                                parent_colname = "HierarchyLevel2") {

  stopifnot(celltype_colname %in% colnames(SingleCellExperiment::colData(spe)))
  stopifnot(parent_colname %in% colnames(SingleCellExperiment::colData(spe)))

  # create mean intensites for each marker in each cell type
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

  stopifnot(nrow(mean_intensities) > 0)

  # convert to long format and scale for plotting
  long_mean_intensities <- mean_intensities %>%
    tidyr::pivot_longer(-c(parent_colname, celltype_colname, "count"),
                        names_to = "marker",
                        values_to = "mean") %>%
    dplyr::group_by(marker) %>%
    dplyr::mutate(standardised_mean = scale(mean)[, 1]) %>%
    dplyr::filter(marker %in% markers)

  # order factors according to arguments
  parenttype_factors <- factor(long_mean_intensities %>% pull(parent_colname),
                               levels = parenttypes)
  celltype_factors <- factor(long_mean_intensities %>% pull(celltype_colname),
                             levels = rev(celltypes))
  marker_factors <- factor(long_mean_intensities %>% pull("marker"),
                           levels = markers)
  long_mean_intensities[, parent_colname] <- parenttype_factors
  long_mean_intensities[, celltype_colname] <- celltype_factors
  long_mean_intensities[, "marker"] <- marker_factors

  raster_plot <- ggplot2::ggplot(long_mean_intensities,
                                 ggplot2::aes(x = marker,
                                              y = !!as.name(celltype_colname),
                                              fill = standardised_mean)) +
    ggplot2::geom_raster() +
    ggplot2::geom_tile(ggplot2::aes(fill = standardised_mean),
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
