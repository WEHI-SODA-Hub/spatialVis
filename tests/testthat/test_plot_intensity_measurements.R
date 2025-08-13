test_that("plot_intesity_measurements() works", {
  # Load the data
  seg_intensity_plot <- plot_intensity_measurements(
    "../../inst/extdata/test_segmentations.geojson"
  )

  # Check that the data is a data frame
  expect_is(seg_intensity_plot, "ggplot")
})
