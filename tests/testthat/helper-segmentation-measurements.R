write_seg_measurements_geojson <- function(path) {
  geojson <- list(
    type = "FeatureCollection",
    features = list(
      list(
        type = "Feature",
        id = "annotation-1",
        geometry = list(type = "Polygon", coordinates = list(list(
          c(0, 0), c(1, 0), c(1, 1), c(0, 1), c(0, 0)
        ))),
        properties = list(objectType = "annotation")
      ),
      list(
        type = "Feature",
        id = "cell-1",
        geometry = list(type = "Polygon", coordinates = list(list(
          c(0, 0), c(2, 0), c(2, 2), c(0, 2), c(0, 0)
        ))),
        properties = list(
          objectType = "cell",
          measurements = list(
            "Cell: Area \u00b5m^2" = 10,
            "Cell: Length \u00b5m" = 20,
            "Cell: Circularity" = 0.5,
            "Cell: Solidity" = 0.8,
            "Cell: Max diameter \u00b5m" = 6,
            "Cell: Min diameter \u00b5m" = 3,
            "Nucleus: Area \u00b5m^2" = 4,
            "Nucleus: Length \u00b5m" = 10,
            "Nucleus: Circularity" = 0.6,
            "Nucleus: Solidity" = 0.85,
            "Nucleus: Max diameter \u00b5m" = 4,
            "Nucleus: Min diameter \u00b5m" = 2,
            "Nucleus/Cell area ratio" = 0.4,
            "DAPI: Cell: Mean" = 110,
            "DAPI: Nucleus: Mean" = 90,
            "Nucleus: DAPI: Mean" = 999,
            "DAPI: Cell: Percentile: 70.0" = 100,
            "DAPI: Nucleus: Percentile: 70.0" = 80,
            "Cell: ErosionBin_1: Area_px" = 40,
            "Cell: ErosionBin_1: Area_Fraction" = 0.4,
            "Cell: ExpansionBin_2: Depth_px" = 5,
            "Nucleus: ErosionBin_1: Area_px" = 20,
            "Neighbours: Mean: Cell: ExpansionBin_3: Area_px" = 4446
          )
        )
      ),
      list(
        type = "Feature",
        id = "cell-2",
        geometry = list(type = "Polygon", coordinates = list(list(
          c(2, 2), c(4, 2), c(4, 4), c(2, 4), c(2, 2)
        ))),
        properties = list(
          objectType = "cell",
          measurements = list(
            "Cell: Area \u00b5m^2" = 14,
            "Cell: Length \u00b5m" = 24,
            "Cell: Circularity" = 0.7,
            "Cell: Solidity" = 0.9,
            "Cell: Max diameter \u00b5m" = 7,
            "Cell: Min diameter \u00b5m" = 4,
            "Nucleus: Area \u00b5m^2" = 6,
            "Nucleus: Length \u00b5m" = 12,
            "Nucleus: Circularity" = 0.65,
            "Nucleus: Solidity" = 0.87,
            "Nucleus: Max diameter \u00b5m" = 4.5,
            "Nucleus: Min diameter \u00b5m" = 2.5,
            "Nucleus/Cell area ratio" = 0.43,
            "DAPI: Cell: Mean" = 130,
            "DAPI: Nucleus: Mean" = 100,
            "Nucleus: DAPI: Mean" = 999,
            "DAPI: Cell: Percentile: 70.0" = 120,
            "DAPI: Nucleus: Percentile: 70.0" = 95,
            "Cell: ErosionBin_1: Area_px" = 50,
            "Cell: ErosionBin_1: Area_Fraction" = 0.5,
            "Cell: ExpansionBin_2: Depth_px" = 7,
            "Nucleus: ErosionBin_1: Area_px" = 24,
            "Neighbours: Mean: Cell: ExpansionBin_3: Area_px" = 5555
          )
        )
      )
    )
  )

  jsonlite::write_json(geojson, path = path, auto_unbox = TRUE, pretty = TRUE)
}
