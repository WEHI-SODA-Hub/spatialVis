test_that("plot_compartment_measurements() works", {
  measurement_data <- get_segmentation_measurements(
    "../../inst/extdata/segmentation.geojson"
  )

  measurement_plot <- plot_compartment_measurements(measurement_data)

  # Check that the data is a data frame
  expect_is(measurement_plot, "ggplot")
})

test_that("plot_compartment_measurements() ignores extra metric columns", {
  geojson_file <- tempfile(fileext = ".geojson")
  on.exit(unlink(geojson_file), add = TRUE)
  write_seg_measurements_geojson(geojson_file)

  measurement_data <- get_segmentation_measurements(
    geojson_file,
    keep_extra_measurements = TRUE
  )
  measurement_plot <- plot_compartment_measurements(measurement_data)

  expect_false(any(grepl("ErosionBin|ExpansionBin|Percentile",
                         measurement_plot$data$measurement)))
  expect_true(any(grepl("^Area", measurement_plot$data$measurement)))
})
