library(dplyr)

#' Filter a SingleCellExperiment object by metadata
#'
#' @description Filter a SingleCellExperiment object by metadata
#' @param spe SingleCellExperiment object
#' @param metadata_col Column name in colData of spe to filter by
#' @param metadata_values List of the metadata column values to perform
#' filtering on
#' @param invert Logical indicating whether to invert the filter (default:
#' FALSE)
#' @return Filtered SingleCellExperiment object
#' @export
#' @importFrom dplyr %>%

filter_spe_by_metadata <- function(spe, metadata_col, metadata_values,
                                   invert = FALSE) {
  stopifnot(metadata_col %in% colnames(SingleCellExperiment::colData(spe)))

  # Extract col data (metadata)
  col_data <- SingleCellExperiment::colData(spe) %>%
    as.data.frame()

  # Get filter indices
  if (invert) {
    filter_indices <- !col_data[[metadata_col]] %in% metadata_values
  } else {
    filter_indices <- col_data[[metadata_col]] %in% metadata_values
  }

  # Expression data
  counts <- SingleCellExperiment::counts(spe)[, filter_indices]

  # Spatial coordinates
  coords <- SpatialExperiment::spatialCoords(spe)[filter_indices, ]

  fspe <- SpatialExperiment::SpatialExperiment(
    assays = counts,
    spatialCoords = coords,
    colData = col_data[filter_indices, ],
    rowData = rownames(spe)
  )
  fspe
}
