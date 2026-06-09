write_seg_geometry_geojson <- function(path, include_nucleus_geometry = TRUE) {
  cell_polygons <- list(
    list(c(0, 0), c(2, 0), c(2, 2), c(0, 2), c(0, 0)),
    list(c(3, 3), c(5, 3), c(5, 5), c(3, 5), c(3, 3))
  )

  nucleus_polygons <- list(
    list(c(0.5, 0.5), c(1.5, 0.5), c(1.5, 1.5), c(0.5, 1.5), c(0.5, 0.5)),
    list(c(3.5, 3.5), c(4.5, 3.5), c(4.5, 4.5), c(3.5, 4.5), c(3.5, 3.5))
  )

  make_feature <- function(i) {
    feature <- list(
      type = "Feature",
      id = paste0("cell-", i),
      geometry = list(type = "Polygon", coordinates = list(cell_polygons[[i]])),
      properties = list(objectType = "cell")
    )

    if (include_nucleus_geometry) {
      feature$nucleusGeometry <- list(
        type = "Polygon",
        coordinates = list(nucleus_polygons[[i]])
      )
    }

    feature
  }

  geojson <- list(
    type = "FeatureCollection",
    features = list(
      make_feature(1),
      make_feature(2)
    )
  )

  jsonlite::write_json(geojson, path = path, auto_unbox = TRUE, pretty = TRUE)
}
