library(dplyr)

#' Plot cell type proportions
#'
#' @description Plot the proportions of cell types in each sample
#' @param spe SingleCellExperiment object containing cell type information
#' @param celltypes Character vector of cell types to plot and the order in
#' which they should be plotted (default: plot all cell types)
#' @param celltype_colname Column name in colData containing cell type
#' information (default: "HierarchyLevel4")
#' @param stack If set to TRUE will create a stacked proportion plot per sample,
#' otherwise will create a standard bar plot per sample (default: TRUE).
#' @return ggplot object
#' @export
#' @importFrom dplyr %>%
plot_cell_props <- function(spe,
                            parenttypes = NULL,
                            parent_colname = "HierarchyLevel1",
                            celltypes = NULL,
                            celltype_colname = "HierarchyLevel4",
                            stack = TRUE) {
  stopifnot(parent_colname %in% colnames(SingleCellExperiment::colData(spe)))
  stopifnot(celltype_colname %in% colnames(SingleCellExperiment::colData(spe)))

  # get all parent and cell types if not provided
  if (is.null(parenttypes)) {
    parenttypes <- unique(SingleCellExperiment::colData(spe)[[parent_colname]])
  }
  if (is.null(celltypes)) {
    celltypes <- unique(SingleCellExperiment::colData(spe)[[celltype_colname]])
  }

  membership_props <- SingleCellExperiment::colData(spe) %>%
    as.data.frame() %>%
    dplyr::filter(!!as.name(celltype_colname) %in% celltypes) %>%
    dplyr::filter(!!as.name(parent_colname) %in% parenttypes) %>%
    dplyr::group_by_at(c("sample_id", celltype_colname)) %>%
    dplyr::summarise(count = n()) %>%
    dplyr::group_modify(~{
      proportion <- .x$count / sum(.x$count)
      cbind(.x, proportion)
    })

  nlevels <- length(celltypes)
  pal <- colorRampPalette(RColorBrewer::brewer.pal(n = 11,
                                                   name = "Spectral"))(nlevels)

  if (stack) {
    prop_plot <- ggplot2::ggplot(data = membership_props,
                                 ggplot2::aes(x = proportion, # nolint: object_usage_linter, line_length_linter.
                                              y = sample_id, # nolint: object_usage_linter, line_length_linter.
                                              fill = !!as.name(celltype_colname))) + # nolint: line_length_linter.
      ggplot2::geom_bar(position = "stack", stat = "identity") +
      ggplot2::scale_fill_manual(values = pal) +
      ggplot2::geom_text(ggplot2::aes(label = round(proportion, 2)), size = 3,
                         position = ggplot2::position_stack(vjust = 0.5)) +
      ggplot2::theme(axis.title.y = ggplot2::element_blank(),
                     axis.text.x = ggplot2::element_text(angle = 90, hjust = 1),
                     axis.title.x = ggplot2::element_blank(),
                     strip.text.y = ggplot2::element_blank(),
                     panel.background = ggplot2::element_blank())
  } else {
    prop_plot <- ggplot2::ggplot(data = membership_props,
                                 ggplot2::aes(x = !!as.name(celltype_colname), # nolint: line_length_linter.
                                              y = proportion)) + # nolint: object_usage_linter, line_length_linter.
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::facet_wrap(~sample_id) +
      ggplot2::theme(axis.title.y = ggplot2::element_text(),
                     axis.text.x = ggplot2::element_text(angle = 90, hjust = 1),
                     axis.title.x = ggplot2::element_blank(),
                     strip.text.y = ggplot2::element_blank(),
                     panel.background = ggplot2::element_blank())
  }
  return(prop_plot)
}
