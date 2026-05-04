test_that("summarise_segmentation_stats() works", {
  measurement_data <- get_segmentation_measurements(
    "../../inst/extdata/segmentation.geojson"
  )

  seg_stats <- summarise_segmentation_stats(measurement_data)

  # Check that the data is a data frame
  expect_is(seg_stats, "data.frame")

  # Check that the data has the correct number of rows
  expect_equal(nrow(seg_stats), 13)
})

test_that("summarise_segmentation_stats() includes new metrics", {
  geojson_file <- tempfile(fileext = ".geojson")
  on.exit(unlink(geojson_file), add = TRUE)
  write_seg_measurements_geojson(geojson_file)

  measurement_data <- get_segmentation_measurements(
    geojson_file,
    keep_extra_measurements = TRUE
  )
  seg_stats <- summarise_segmentation_stats(
    measurement_data,
    keep_extra_measurements = TRUE
  )

  expect_true("DAPI: Cell: Percentile: 70.0" %in% seg_stats$Measurement)
  expect_true("DAPI: Nucleus: Percentile: 70.0" %in% seg_stats$Measurement)
  expect_true("Cell: ErosionBin_1: Area_px" %in% seg_stats$Measurement)
  expect_true("Cell: ErosionBin_1: Area_Fraction" %in% seg_stats$Measurement)
  expect_true("Cell: ExpansionBin_2: Depth_px" %in% seg_stats$Measurement)
  expect_true("Nucleus: ErosionBin_1: Area_px" %in% seg_stats$Measurement)
  expect_false(any(startsWith(seg_stats$Measurement, "Neighbours: ")))
  expect_equal(
    seg_stats$Mean[seg_stats$Measurement == "DAPI: Cell: Percentile: 70.0"],
    110
  )
  expect_equal(
    seg_stats$Median[seg_stats$Measurement == "Cell: ErosionBin_1: Area_px"],
    45
  )
})

test_that("summarise_segmentation_stats() drops new metrics by default", {
  geojson_file <- tempfile(fileext = ".geojson")
  on.exit(unlink(geojson_file), add = TRUE)
  write_seg_measurements_geojson(geojson_file)

  measurement_data <- get_segmentation_measurements(
    geojson_file,
    keep_extra_measurements = TRUE
  )
  seg_stats <- summarise_segmentation_stats(measurement_data)

  expect_false("DAPI: Cell: Percentile: 70.0" %in% seg_stats$Measurement)
  expect_false("Cell: ErosionBin_1: Area_px" %in% seg_stats$Measurement)
  expect_false(any(startsWith(seg_stats$Measurement, "Neighbours: ")))
})
