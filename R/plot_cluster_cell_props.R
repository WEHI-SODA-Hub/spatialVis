library(dplyr)

#' Plot cell type proportions in each cluster
#'
#' @description Plot the proportions of cell types in each cluster
#' @param spe SingleCellExperiment object containing cluster information
#' @param celltypes Character vector of cell types to plot and the order in
#' which they should be plotted (default: plot all cell types)
#' @param celltype_colname Column name in colData containing cell type
#' information (default: "HierarchyLevel4")
#' @return ggplot object
#' @export
#' @importFrom dplyr %>%
plot_cluster_cell_props <- function(spe,
                                    celltypes = NULL,
                                    celltype_colname = "HierarchyLevel4") {
  stopifnot("cluster" %in% colnames(SingleCellExperiment::colData(spe)))

  # get all cell types if not provided
  if (is.null(celltypes)) {
    celltypes <- unique(SingleCellExperiment::colData(spe)[[celltype_colname]])
  }
  # summarise the cluster memberships
  membership_props <- SingleCellExperiment::colData(spe) %>%
    as.data.frame() %>%
    dplyr::group_by_at(c("cluster", celltype_colname)) %>%
    dplyr::summarise(count = n()) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(cluster) # nolint: object_usage_linter.

  # calculate the proportions in each cluster
  cluster_membership_props <- membership_props %>%
    dplyr::group_by_at(c("cluster")) %>%
    dplyr::group_modify(~{
      proportion <- .x$count / sum(.x$count)
      cbind(.x, proportion)
    })

  # add the levels
  cluster_levels <- factor(cluster_membership_props %>%
                             pull(celltype_colname), levels = celltypes)
  cluster_membership_props[, celltype_colname] <- cluster_levels

  nlevels <- levels(cluster_levels) %>% length()
  pal <- colorRampPalette(RColorBrewer::brewer.pal(n = 11,
                                                   name = "Spectral"))(nlevels)

  # plot the results
  prop_plot <- ggplot2::ggplot(data = cluster_membership_props,
                               ggplot2::aes(x = as.factor(cluster), y = count, # nolint: object_usage_linter, line_length_linter.
                                            fill = !!as.name(celltype_colname))) + # nolint: line_length_linter.
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::geom_text(ggplot2::aes(label = count), size = 3,
                       position = ggplot2::position_stack(vjust = 0.5)) +
    ggplot2::scale_fill_manual(values = pal) +
    ggplot2::theme(axis.title.y = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(angle = 90, hjust = 1),
                   axis.title.x = ggplot2::element_blank(),
                   strip.text.y = ggplot2::element_blank(),
                   panel.background = ggplot2::element_blank())

  return(prop_plot)
}
