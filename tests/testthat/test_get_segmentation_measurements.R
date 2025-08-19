test_that("get_segmentation_measurements() works", {
  measurement_data <- get_segmentation_measurements(
    "../../inst/extdata/test_segmentations.geojson"
  )
  # Check that the data is a data frame
  expect_is(measurement_data, "data.table")
})
