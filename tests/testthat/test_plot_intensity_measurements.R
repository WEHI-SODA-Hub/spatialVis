test_that("plot_intensity_measurements() works", {
  measurement_data <- get_segmentation_measurements(
    "../../inst/extdata/segmentation.geojson"
  )

  seg_intensity_plot <- plot_intensity_measurements(measurement_data)

  expect_is(seg_intensity_plot, "ggplot")
})
