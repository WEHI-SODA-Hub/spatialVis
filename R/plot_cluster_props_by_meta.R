library(dplyr)

#' Plot sample proportions by some meta-value (e.g., sample_id)
#'
#' @description Plot the proportions of cell types in each cluster
#' @param spe SingleCellExperiment object containing cluster information
#' @param meta Column name in colData containing meta information (default:
#' "sample_id")
#' @return ggplot object
#' @export
#' @importFrom dplyr %>%
plot_cluster_props_by_meta <- function(spe, meta_col = "sample_id") {
  stopifnot("cluster" %in% colnames(SingleCellExperiment::colData(spe)))
  stopifnot(meta_col %in% colnames(SingleCellExperiment::colData(spe)))

  # summarise the cluster memberships
  cluster_membership_props <- SingleCellExperiment::colData(spe) %>%
    as.data.frame() %>%
    dplyr::group_by_at(c("cluster", meta_col)) %>%
    dplyr::summarise(count = dplyr::n()) %>%
    dplyr::ungroup()

  cluster_levels <- factor(cluster_membership_props %>%
                             dplyr::pull(cluster)) # nolint: object_usage_linter, line_length_linter.
  cluster_membership_props[, "cluster"] <- cluster_levels

  nlevels <- levels(cluster_levels) %>% length()
  pal <- colorRampPalette(RColorBrewer::brewer.pal(n = 11,
                                                   name = "Spectral"))(nlevels)

  prop_plot <- ggplot2::ggplot(data = cluster_membership_props,
                               ggplot2::aes(x = count, y = !!as.name(meta_col),
                                            fill = factor(cluster))) + # nolint: object_usage_linter, line_length_linter.
    ggplot2::geom_bar(position = "fill", stat = "identity") +
    ggplot2::scale_fill_manual(values = pal) +
    ggplot2::theme(axis.title.y = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(angle = 90, hjust = 1),
                   axis.title.x = ggplot2::element_blank(),
                   strip.text.y = ggplot2::element_blank(),
                   panel.background = ggplot2::element_blank())

  prop_plot
}
