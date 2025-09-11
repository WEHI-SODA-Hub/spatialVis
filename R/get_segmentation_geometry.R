#' @title Get segmentation geometry from GeoJSON file
#'
#' @description This function reads a GeoJSON file containing cell segmentation
#' geometry and returns a list with the coordinates of the cell and nucleus
#' geometries. If `only_keep_cells` is TRUE, it will filter the data to only
#' keep features with objectType "cell". It also filters out multipolygon
#' geometries, as they are not supported in the current plotting function.
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

  is_polygon <- seg$features$geometry$type == "Polygon" &
    seg$features$nucleusGeometry$type == "Polygon"

  if (only_keep_cells) {
    # Keep only object type cells with polygon geometries
    is_cell <- seg$features$properties$objectType == "cell"
    keep <- is_polygon & is_cell

    geom_data$cell <- geom_data$cell[keep]
    geom_data$nucleus <- geom_data$nucleus[keep]
  } else {
    # Just filter out non-polygon geometries
    geom_data$cell <- geom_data$cell[is_polygon]
    geom_data$nucleus <- geom_data$nucleus[is_polygon]
  }

  message(paste0(
    "Filtered out ", sum(!is_polygon, na.rm = TRUE),
    " non-polygon geometries.\n",
    "Final number of cell geometries: ", length(geom_data$cell)
  ))

  # Clear seg as we don't need other properties
  rm(seg)

  geom_data
}
