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
#' @param colour_by Column name in colData of spe to colour points by (options:
#' "cell_type" (default), "cluster" or "tumour")
#' @return ggplot2 object
#' @export
#' @importFrom dplyr %>%
plot_cells_spatially <- function(spe,
                                 cell_type_colname = "HierarchyLevel4",
                                 parent_types = NULL,
                                 parent_type_colname = "HierarchyLevel1",
                                 centroid_x_col = "Cell.X.Position",
                                 centroid_y_col = "Cell.Y.Position",
                                 pointsize = 0.5,
                                 colour_by = "cell_type") {
  stopifnot(cell_type_colname %in% colnames(SingleCellExperiment::colData(spe)))
  stopifnot(parent_type_colname %in% colnames(SingleCellExperiment::colData(spe))) # nolint: line_length_linter.
  if (colour_by == "cluster") {
    stopifnot("cluster" %in% colnames(SingleCellExperiment::colData(spe)))
  }
  if (colour_by == "tumour") {
    stopifnot("In.Tumour" %in% colnames(SingleCellExperiment::colData(spe)))
  }

  if (is.null(parent_types)) {
    parent_types <- unique(SingleCellExperiment::colData(spe)[[parent_type_colname]]) # nolint: line_length_linter.
  }

  coords <- SpatialExperiment::spatialCoords(spe) %>% as.data.frame()
  centroid_x_col <- grep(centroid_x_col, colnames(coords), value = TRUE)
  centroid_y_col <- grep(centroid_y_col, colnames(coords), value = TRUE)

  stopifnot(length(centroid_x_col) == 1 && length(centroid_y_col) == 1)

  spe_parent_types <- SingleCellExperiment::colData(spe)[[parent_type_colname]]
  if (colour_by == "cluster") {
    spe_clusters <- SingleCellExperiment::colData(spe)[["cluster"]]
    df <- dplyr::mutate(coords, cluster = as.factor(spe_clusters),
                        parent_type = spe_parent_types)

    ncols <- df$cluster %>% unique() %>% length()
  } else if (colour_by == "tumour") {
    spe_tumour <- SingleCellExperiment::colData(spe)[["In.Tumour"]]
    df <- dplyr::mutate(coords, tumour = as.factor(spe_tumour),
                        parent_type = spe_parent_types)

    ncols <- df$tumour %>% unique() %>% length()
  } else {
    spe_cell_types <- SingleCellExperiment::colData(spe)[[cell_type_colname]]
    df <- dplyr::mutate(coords, cell_type = spe_cell_types,
                        parent_type = spe_parent_types)
    ncols <- df$cell_type %>% unique() %>% length()
  }

  df <- dplyr::filter(df, parent_type %in% parent_types) # nolint: object_usage_linter, line_length_linter.
  pal <- colorRampPalette(RColorBrewer::brewer.pal(n = 11,
                                                   name = "Spectral"))(ncols)

  cellplot <- ggplot2::ggplot(df,
                              ggplot2::aes(x = !!as.name(centroid_x_col),
                                           y = !!as.name(centroid_y_col))) +
    ggplot2::geom_point(ggplot2::aes(color = !!as.name(colour_by)),
                        size = pointsize) + # nolint: object_usage_linter.
    ggplot2::theme_minimal() +
    ggplot2::scale_colour_manual(values = pal)

  cellplot
}
