library(dplyr)

#' Plot cell type proportions in each cluster
#'
#' @description Plot the proportions of cell types in each cluster
#' @param spe SingleCellExperiment object containing cluster information
#' @param cell_types Character vector of cell types to plot and the order in
#' which they should be plotted (default: plot all cell types)
#' @param cell_type_colname Column name in colData containing cell type
#' information (default: "HierarchyLevel4")
#' @param parent_types Character vector of parent cell types to include
#' (default: all parent types)
#' @param exclude_parent_types Character vector of parent cell types to exclude
#' @param plot_type Type of plot to generate ("bar" (default) or "heatmap")
#' @return ggplot object
#' @export
#' @importFrom dplyr %>%
plot_cluster_cell_props <- function(spe,
                                    cell_types = NULL,
                                    cell_type_colname = "HierarchyLevel4",
                                    parent_types = NULL,
                                    parent_type_colname = "HierarchyLevel1",
                                    exclude_parent_types = NULL,
                                    plot_type = "bar") {
  stopifnot("cluster" %in% colnames(SingleCellExperiment::colData(spe)))
  stopifnot(cell_type_colname %in% colnames(SingleCellExperiment::colData(spe)))
  stopifnot(plot_type %in% c("bar", "heatmap"))

  # get all parent and cell types if not provided
  if (is.null(parent_types)) {
    parent_types <-
      unique(SingleCellExperiment::colData(spe)[[parent_type_colname]])
  }
  if (is.null(cell_types)) {
    cell_types <-
      unique(SingleCellExperiment::colData(spe)[[cell_type_colname]])
  }

  # summarise the cluster memberships
  cluster_membership_props <- SingleCellExperiment::colData(spe) %>%
    as.data.frame() %>%
    dplyr::filter(!!as.name(parent_type_colname) %in% parent_types) %>%
    dplyr::filter(!(!!as.name(parent_type_colname) %in%
                      exclude_parent_types)) %>%
    dplyr::group_by_at(c("cluster", cell_type_colname)) %>%
    dplyr::summarise(count = dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::group_by_at(c("cluster")) %>%
    dplyr::group_modify(~{
      proportion <- .x$count / sum(.x$count)
      cbind(.x, proportion)
    })

  # add the levels
  cluster_levels <- factor(cluster_membership_props %>%
                             dplyr::pull(cell_type_colname),
                           levels = cell_types)
  cluster_membership_props[, cell_type_colname] <- cluster_levels

  nlevels <- levels(cluster_levels) %>% length()
  pal <- colorRampPalette(RColorBrewer::brewer.pal(n = 11,
                                                   name = "Spectral"))(nlevels)

  # plot the results
  if (plot_type == "heatmap") {
    prop_plot <- ggplot2::ggplot(cluster_membership_props,
                                 ggplot2::aes(x = !!as.name(cell_type_colname), # nolint: object_usage_linter, line_length_linter.
                                              y = factor(cluster))) + # nolint: object_usage_linter, line_length_linter.
      ggplot2::geom_raster() +
      ggplot2::geom_tile(ggplot2::aes(fill = proportion), colour = "black", # nolint: object_usage_linter, line_length_linter.
                         linewidth = 0.5) +
      ggplot2::scale_fill_distiller(palette = "YlGnBu") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45,
                                                         hjust = 1)) +
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
  } else {
    prop_plot <- ggplot2::ggplot(data = cluster_membership_props,
                                 ggplot2::aes(x = as.factor(cluster), y = count, # nolint: object_usage_linter, line_length_linter.
                                              fill = !!as.name(cell_type_colname))) + # nolint: line_length_linter.
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::geom_text(ggplot2::aes(label = count), size = 3,
                         position = ggplot2::position_stack(vjust = 0.5)) +
      ggplot2::scale_fill_manual(values = pal) +
      ggplot2::theme(axis.title.y = ggplot2::element_blank(),
                     axis.text.x = ggplot2::element_text(angle = 90, hjust = 1),
                     axis.title.x = ggplot2::element_blank(),
                     strip.text.y = ggplot2::element_blank(),
                     panel.background = ggplot2::element_blank())
  }
  return(prop_plot)
}
