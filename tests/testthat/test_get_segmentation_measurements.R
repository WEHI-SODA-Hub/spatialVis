test_that("get_segmentation_measurements() works", {
  measurement_data <- get_segmentation_measurements(
    "../../inst/extdata/segmentation.geojson"
  )
  # Check that the data is a data frame
  expect_is(measurement_data, "data.table")
})

test_that("get_segmentation_measurements() drops extra geojson metric columns by default", { # nolint
  geojson_file <- tempfile(fileext = ".geojson")
  on.exit(unlink(geojson_file), add = TRUE)
  write_seg_measurements_geojson(geojson_file)

  measurement_data <- get_segmentation_measurements(geojson_file)

  expect_false("DAPI: Cell: Percentile: 70.0" %in% colnames(measurement_data))
  expect_false("Cell: ErosionBin_1: Area_px" %in% colnames(measurement_data))
  expect_false(
    "Neighbours: Mean: Cell: ExpansionBin_3: Area_px" %in%
      colnames(measurement_data)
  )
})

test_that("get_segmentation_measurements() keeps extra geojson metric columns when requested", { # nolint
  geojson_file <- tempfile(fileext = ".geojson")
  on.exit(unlink(geojson_file), add = TRUE)
  write_seg_measurements_geojson(geojson_file)

  measurement_data <- get_segmentation_measurements(
    geojson_file,
    keep_extra_measurements = TRUE
  )

  expect_true("DAPI: Cell: Percentile: 70.0" %in% colnames(measurement_data))
  expect_true("Cell: ErosionBin_1: Area_px" %in% colnames(measurement_data))
  expect_true("Cell: ExpansionBin_2: Depth_px" %in% colnames(measurement_data))
  expect_true(
    "Neighbours: Mean: Cell: ExpansionBin_3: Area_px" %in%
      colnames(measurement_data)
  )
})
