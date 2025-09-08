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
