library(dplyr)

#' Plot cell types on spatial coordinates
#'
#' @param spe SingleCellExperiment object
#' @param celltype_colname Column name in colData of spe that contains cell
#' types (default: "HierarchyLevel4")
#' @param centroid_x_col Column name in spatial coordinates that contains x
#' coordinates (default: "Centroid X")
#' @param centroid_y_col Column name in spatial coordinates that contains y
#' coordinates (default: "Centroid Y")
#' @param pointsize Size of points in the plot (default: 0.5)
#' @return ggplot2 object
#' @export
#' @importFrom dplyr %>%
plot_celltypes <- function(spe, celltype_colname = "HierarchyLevel4",
                           centroid_x_col = "Cell.X.Position",
                           centroid_y_col = "Cell.Y.Position",
                           pointsize = 0.5) {
  stopifnot(celltype_colname %in% colnames(SingleCellExperiment::colData(spe)))

  coords <- SpatialExperiment::spatialCoords(spe) %>% as.data.frame()
  centroid_x_col <- grep(centroid_x_col, colnames(coords), value = TRUE)
  centroid_y_col <- grep(centroid_y_col, colnames(coords), value = TRUE)

  stopifnot(length(centroid_x_col) == 1 && length(centroid_y_col) == 1)

  celltypes <- SingleCellExperiment::colData(spe)[[celltype_colname]]
  df <- dplyr::mutate(coords, celltype = celltypes)

  ncells <- df$celltype %>% unique() %>% length()
  pal <- colorRampPalette(RColorBrewer::brewer.pal(n = 11,
                                                   name = "Spectral"))(ncells)

  cellplot <- ggplot2::ggplot(df,
                              ggplot2::aes(x = !!as.name(centroid_x_col),
                                           y = !!as.name(centroid_y_col))) +
    ggplot2::geom_point(ggplot2::aes(color = celltype), size = pointsize) + # nolint: object_usage_linter, line_length_linter.
    ggplot2::theme_minimal() +
    ggplot2::scale_colour_manual(values = pal)

  return(cellplot)
}
