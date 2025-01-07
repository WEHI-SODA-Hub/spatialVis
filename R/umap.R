library(dplyr)

#' @title Plot UMAP
#'
#' @description Plot a UMAP of cell types
#' @param spe SingleCellExperiment object containing marker intensities
#' @param markers Character vector of markers to plot and the order in which
#' they should be plotted
#' @param celltypes Character vector of cell types to plot and the order in
#' which they should be plotted
#' @param celltype_colname Column name in colData containing cell type
#' information (default: "HierarchyLevel4")
#' @param parent_colname Column name in colData containing parent cell type
#' information (default: "HierarchyLevel2")
#' @return A ggplot2 object
#' @export
#' @importFrom dplyr %>%
plot_umap <- function(spe, markers = NULL, celltypes = NULL,
                      celltype_colname = "HierarchyLevel4",
                      parent_colname = "HierarchyLevel2") {

  stopifnot(celltype_colname %in% colnames(SingleCellExperiment::colData(spe)))

  if (is.null(celltypes)) {
    celltypes <- unique(SingleCellExperiment::colData(spe)[[celltype_colname]])
  }
  if (is.null(markers)) {
    markers <- rownames(spe)
  }

  celltype_intensities <- SingleCellExperiment::colData(spe) %>%
    as.data.frame() %>%
    dplyr::select(dplyr::all_of(c(celltype_colname, parent_colname))) %>%
    cbind(t(SingleCellExperiment::counts(spe))) %>%
    dplyr::filter(!!as.name(celltype_colname) %in% celltypes) %>%
    dplyr::mutate(dplyr::across(dplyr::where(is.numeric),
                                ~tidyr::replace_na(., 0)))

  # run umap
  system.time({
    cell_umap <- celltype_intensities %>%
      dplyr::select(-dplyr::where(is.character)) %>%
      umap::umap(.) # nolint: object_usage_linter.
  })

  # Extract UMAP coordinates
  umap_coordinates <- cell_umap$layout

  # Convert to data frame
  umap_df <- as.data.frame(umap_coordinates)

  umap_df$cell_type <- factor(celltype_intensities[[celltype_colname]],
                              levels = celltypes)

  p <- ggplot2::ggplot(as.data.frame(umap_df),
                       ggplot2::aes(x = V1, y = V2, color = cell_type)) + # nolint: object_usage_linter, line_length_linter.
    ggplot2::geom_point(size = 0.5, stroke = 0) +
    ggplot2::theme(panel.background = ggplot2::element_blank()) +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 4))) # nolint: line_length_linter.

  return(p)
}
