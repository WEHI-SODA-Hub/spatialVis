library(dplyr)

#' @title Get segmentation geometry from GeoJSON file
#'
#' @description This function reads a GeoJSON file containing cell segmentation
#' geometry and returns a list with the coordinates of the cell and nucleus
#' geometries. If `only_keep_cells` is TRUE, it will filter the data to only
#' keep features with objectType "cell".
#'
#' @param geojson_file Path to the GeoJSON file containing segmentation data
#' @param only_keep_cells Logical indicating whether to filter the data to
#' only keep cells (default: TRUE). If TRUE, only features with objectType
#' "cell" will be kept.
#'
#' @return A list with two elements: `cell` and `nucleus`, each containing a
#' list of geometries.
#' @export
#' @importFrom dplyr %>%
get_segmentation_geometry <- function(geojson_file,
                                      only_keep_cells = TRUE) {
  stopifnot(file.exists(geojson_file))

  seg <- jsonlite::fromJSON(geojson_file)

  stopifnot(
    "features" %in% names(seg),
    "properties" %in% names(seg$features),
    "objectType" %in% names(seg$features$properties),
    "geometry" %in% names(seg$features),
    "nucleusGeometry" %in% names(seg$features),
    "coordinates" %in% names(seg$features$geometry),
    "coordinates" %in% names(seg$features$nucleusGeometry)
  )

  geom_data <- list()
  geom_data$cell <- seg$features$geometry$coordinates
  geom_data$nucleus <- seg$features$nucleusGeometry$coordinates

  if (only_keep_cells) {
    is_cell <- seg$features$properties$objectType == "cell"
    geom_data$cell <- geom_data$cell[is_cell]
    geom_data$nucleus <- geom_data$nucleus[is_cell]
  }

  # Clear seg as we don't need other properties
  rm(seg)

  geom_data
}
