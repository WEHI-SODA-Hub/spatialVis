library(dplyr)

#' Plot cell types on spatial coordinates
#'
#' @param spe SingleCellExperiment object
#' @param cell_type_colname Column name in colData of spe that contains cell
#' types (default: "HierarchyLevel4")
#' @param parent_types Character vector of parent cell types to include
#' (default: all parent types)
#' @param parent_type_colname Column name in colData of spe that contains parent
#' cell types (default: "HierarchyLevel1")
#' @param centroid_x_col Column name in spatial coordinates that contains x
#' coordinates (default: "Centroid X")
#' @param centroid_y_col Column name in spatial coordinates that contains y
#' coordinates (default: "Centroid Y")
#' @param pointsize Size of points in the plot (default: 0.5)
#' @return ggplot2 object
#' @export
#' @importFrom dplyr %>%
plot_cell_types <- function(spe,
                            cell_type_colname = "HierarchyLevel4",
                            parent_types = NULL,
                            parent_type_colname = "HierarchyLevel1",
                            centroid_x_col = "Cell.X.Position",
                            centroid_y_col = "Cell.Y.Position",
                            pointsize = 0.5) {
  stopifnot(cell_type_colname %in% colnames(SingleCellExperiment::colData(spe)))
  stopifnot(parent_type_colname %in% colnames(SingleCellExperiment::colData(spe))) # nolint: line_length_linter.

  if (is.null(parent_types)) {
    parent_types <- unique(SingleCellExperiment::colData(spe)[[parent_type_colname]]) # nolint: line_length_linter.
  }

  coords <- SpatialExperiment::spatialCoords(spe) %>% as.data.frame()
  centroid_x_col <- grep(centroid_x_col, colnames(coords), value = TRUE)
  centroid_y_col <- grep(centroid_y_col, colnames(coords), value = TRUE)

  stopifnot(length(centroid_x_col) == 1 && length(centroid_y_col) == 1)

  spe_cell_types <- SingleCellExperiment::colData(spe)[[cell_type_colname]]
  spe_parent_types <- SingleCellExperiment::colData(spe)[[parent_type_colname]]
  df <- dplyr::mutate(coords, cell_type = spe_cell_types,
                      parent_type = spe_parent_types) %>%
    dplyr::filter(parent_type %in% parent_types) # nolint: object_usage_linter.

  ncells <- df$cell_type %>% unique() %>% length()
  pal <- colorRampPalette(RColorBrewer::brewer.pal(n = 11,
                                                   name = "Spectral"))(ncells)

  cellplot <- ggplot2::ggplot(df,
                              ggplot2::aes(x = !!as.name(centroid_x_col),
                                           y = !!as.name(centroid_y_col))) +
    ggplot2::geom_point(ggplot2::aes(color = cell_type), size = pointsize) + # nolint: object_usage_linter, line_length_linter.
    ggplot2::theme_minimal() +
    ggplot2::scale_colour_manual(values = pal)

  cellplot
}
