library(dplyr)

#' Plot cell type proportions
#'
#' @description Plot the proportions of cell types in each sample
#' @param spe SingleCellExperiment object containing cell type information
#' @param parent_types Character vector of parent cell types to include
#' (default: NULL = all parent types)
#' @param cell_types Character vector of cell types to plot and the order in
#' which they should be plotted (default: NULL = plot all cell types)
#' @param cell_type_colname Column name in colData containing cell type
#' information (default: "HierarchyLevel4")
#' @param stack If set to TRUE will create a stacked proportion plot per sample,
#' otherwise will create a standard bar plot per sample (default: TRUE).
#' @param facet_by Column name in colData to facet the plot by (default:
#' "sample_id")
#' @return ggplot object
#' @export
#' @importFrom dplyr %>%
plot_cell_props <- function(spe,
                            parent_types = NULL,
                            parent_colname = "HierarchyLevel1",
                            cell_types = NULL,
                            cell_type_colname = "HierarchyLevel4",
                            stack = TRUE,
                            facet_by = "sample_id") {
  stopifnot(parent_colname %in% colnames(SingleCellExperiment::colData(spe)))
  stopifnot(cell_type_colname %in% colnames(SingleCellExperiment::colData(spe)))
  stopifnot(facet_by %in% colnames(SingleCellExperiment::colData(spe)))

  # get all parent and cell types if not provided
  if (is.null(parent_types)) {
    parent_types <- unique(SingleCellExperiment::colData(spe)[[parent_colname]])
  }
  if (is.null(cell_types)) {
    cell_types <- unique(SingleCellExperiment::colData(spe)[[cell_type_colname]]) # nolint: line_length_linter.
  }

  membership_props <- SingleCellExperiment::colData(spe) %>%
    as.data.frame() %>%
    dplyr::filter(!!as.name(cell_type_colname) %in% cell_types) %>%
    dplyr::filter(!!as.name(parent_colname) %in% parent_types) %>%
    dplyr::group_by_at(c(facet_by, cell_type_colname)) %>%
    dplyr::summarise(count = dplyr::n()) %>%
    dplyr::group_modify(~{
      proportion <- .x$count / sum(.x$count)
      cbind(.x, proportion)
    })

  nlevels <- length(cell_types)
  pal <- colorRampPalette(RColorBrewer::brewer.pal(n = 11,
                                                   name = "Spectral"))(nlevels)

  if (stack) {
    prop_plot <- ggplot2::ggplot(data = membership_props,
                                 ggplot2::aes(x = proportion, # nolint: object_usage_linter, line_length_linter.
                                              y = !!rlang::sym(facet_by),
                                              fill = !!as.name(cell_type_colname))) + # nolint: line_length_linter.
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
                                 ggplot2::aes(x = !!as.name(cell_type_colname), # nolint: line_length_linter.
                                              y = proportion)) + # nolint: object_usage_linter, line_length_linter.
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::facet_wrap(dplyr::vars(!!rlang::sym(facet_by))) +
      ggplot2::theme(axis.title.y = ggplot2::element_text(),
                     axis.text.x = ggplot2::element_text(angle = 90, hjust = 1),
                     axis.title.x = ggplot2::element_blank(),
                     strip.text.y = ggplot2::element_blank(),
                     panel.background = ggplot2::element_blank())
  }
  return(prop_plot)
}
