test_that("plot_intesity_measurements() works", {
  measurement_data <- get_segmentation_measurements(
    "../../inst/extdata/test_segmentations.geojson"
  )

  seg_intensity_plot <- plot_intensity_measurements(measurement_data)

  # Check that the data is a data frame
  expect_is(seg_intensity_plot, "ggplot")
})
