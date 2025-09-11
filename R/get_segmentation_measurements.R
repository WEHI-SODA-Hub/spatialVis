#' @title Get segmentation data from GeoJSON file
#'
#' @description This function reads a GeoJSON file containing cell segmentation
#' measurements and returns a data.table with the properties of the features.
#' If no measurements are present, it will return a data.table with only the
#' object types.
#'
#' @param geojson_file Path to the GeoJSON file containing segmentation data
#' @param only_keep_cells Logical indicating whether to filter the data to
#' only keep cells (default: TRUE). If TRUE, only features with objectType
#' "cell" will be kept.
#'
#' @return A data.table object
#' @export
#' @importFrom dplyr %>%
get_segmentation_measurements <- function(geojson_file,
                                          only_keep_cells = TRUE) {
  stopifnot(file.exists(geojson_file))

  seg <- jsonlite::fromJSON(geojson_file)

  stopifnot(
    "features" %in% names(seg),
    "properties" %in% names(seg$features),
    "objectType" %in% names(seg$features$properties)
  )

  seg_properties <- seg$features$properties

  # Clear seg as we don't need other properties
  rm(seg)

  if ("measurements" %in% names(seg_properties)) {
    dt <- data.table::data.table(
      objectType = seg_properties$objectType,
      seg_properties$measurements
    )
  } else {
    # If no measurements are available, we'll just get a data.table containing
    # the object types -- this is only useful for counting cells
    dt <- data.table::as.data.table(seg_properties)
  }

  if (only_keep_cells) {
    dt <- dplyr::filter(dt, objectType == "cell") # nolint
  }
  dt
}
