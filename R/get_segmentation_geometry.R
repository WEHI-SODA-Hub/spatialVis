#' @title Get segmentation geometry from GeoJSON file
#'
#' @description This function reads a GeoJSON file containing cell segmentation
#' geometry and returns a list with the coordinates of the cell and nucleus
#' geometries. If nuclear geometries are not present, it will only return cell
#' geometries. If `only_keep_cells` is TRUE, it will filter the data to only
#' keep features with objectType "cell". It also filters out multipolygon
#' geometries, as they are not supported in the current plotting function.
#'
#' @param geojson_file Path to the GeoJSON file containing segmentation data
#' @param only_keep_cells Logical indicating whether to filter the data to
#' only keep cells (default: TRUE). If TRUE, only features with objectType
#' "cell" will be kept.
#' @param nuclei_only Logical indicating whether to treat the primary geometry
#' as nuclear geometry (default: FALSE). When TRUE, `geometry` coordinates are
#' returned in the `nucleus` slot instead of `cell`.
#'
#' @return A list containing primary geometries in `cell` (default mode) or
#' `nucleus` (`nuclei_only = TRUE`). In default mode, a secondary `nucleus`
#' slot is included when `nucleusGeometry` is present.
#' @export
#' @importFrom dplyr %>%
get_segmentation_geometry <- function(geojson_file,
                                      only_keep_cells = TRUE,
                                      nuclei_only = FALSE) {
  stopifnot(file.exists(geojson_file))

  seg <- jsonlite::fromJSON(geojson_file)

  stopifnot(
    "features" %in% names(seg),
    "properties" %in% names(seg$features),
    "objectType" %in% names(seg$features$properties),
    "geometry" %in% names(seg$features),
    "coordinates" %in% names(seg$features$geometry)
  )

  geom_data <- list()
  primary_slot <- if (nuclei_only) "nucleus" else "cell"
  geom_data[[primary_slot]] <- seg$features$geometry$coordinates
  is_polygon <- seg$features$geometry$type == "Polygon"

  whole_cell_only <- FALSE
  if (nuclei_only) {
    # In nuclei-only mode, the primary geometry is already treated as nucleus.
    whole_cell_only <- TRUE
  } else if (!"nucleusGeometry" %in% names(seg$features)) {
    message("No nucleus geometry found in GeoJSON file.")
    whole_cell_only <- TRUE
  } else {
    stopifnot("coordinates" %in% names(seg$features$nucleusGeometry))
    geom_data$nucleus <- seg$features$nucleusGeometry$coordinates
    is_polygon <- is_polygon & seg$features$nucleusGeometry$type == "Polygon"
  }

  has_dims <- sapply(geom_data[[primary_slot]], function(x) {
    !is.null(dim(x))
  })

  n_geom <- length(geom_data[[primary_slot]])
  if (only_keep_cells) {
    # Keep only object type cells with polygon geometries
    is_cell <- seg$features$properties$objectType == "cell"
    keep <- is_polygon & is_cell & has_dims

    geom_data[[primary_slot]] <- geom_data[[primary_slot]][keep]
    if (!whole_cell_only) {
      geom_data$nucleus <- geom_data$nucleus[keep]
    }
  } else {
    # Just filter out non-polygon geometries
    geom_data[[primary_slot]] <- geom_data[[primary_slot]][is_polygon & has_dims] # nolint
    if (!whole_cell_only) {
      geom_data$nucleus <- geom_data$nucleus[is_polygon & has_dims]
    }
  }

  message(paste0(
    "Filtered out ", n_geom - length(geom_data[[primary_slot]]),
    " geometries.\n",
    "Final number of ", primary_slot, " geometries: ",
    length(geom_data[[primary_slot]])
  ))

  # Clear seg as we don't need other properties
  rm(seg)

  geom_data
}
