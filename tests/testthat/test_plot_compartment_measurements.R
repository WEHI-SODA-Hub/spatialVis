test_that("plot_compartment_measurements() works", {
  measurement_data <- get_segmentation_measurements(
    "../../inst/extdata/segmentation.geojson"
  )

  measurement_plot <- plot_intensity_measurements(measurement_data)

  # Check that the data is a data frame
  expect_is(measurement_plot, "ggplot")
})
