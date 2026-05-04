test_that("plot_intensity_measurements() works", {
  measurement_data <- get_segmentation_measurements(
    "../../inst/extdata/segmentation.geojson"
  )

  seg_intensity_plot <- plot_intensity_measurements(measurement_data)

  expect_is(seg_intensity_plot, "ggplot")
})

test_that("plot_intensity_measurements() supports legacy intensity layout", {
  geojson_file <- tempfile(fileext = ".geojson")
  on.exit(unlink(geojson_file), add = TRUE)
  write_seg_measurements_geojson(geojson_file, intensity_layout = "legacy")

  measurement_data <- get_segmentation_measurements(
    geojson_file,
    keep_extra_measurements = TRUE
  )
  seg_intensity_plot <- plot_intensity_measurements(
    measurement_data,
    compartments = c("Cell", "Nucleus"),
    calculation = "Mean"
  )

  expect_equal(sort(unique(seg_intensity_plot$data$channel)), "DAPI")
  expect_equal(
    sort(unique(seg_intensity_plot$data$compartment)),
    c("Cell", "Nucleus")
  )
  expect_equal(nrow(seg_intensity_plot$data), 4)
  expect_false(any(seg_intensity_plot$data$channel == "Nucleus"))
})

test_that("plot_intensity_measurements() supports QuPath >0.6 intensity layout", { # nolint
  geojson_file <- tempfile(fileext = ".geojson")
  on.exit(unlink(geojson_file), add = TRUE)
  write_seg_measurements_geojson(geojson_file, intensity_layout = "qpath06")

  measurement_data <- get_segmentation_measurements(
    geojson_file,
    keep_extra_measurements = TRUE
  )
  seg_intensity_plot <- plot_intensity_measurements(
    measurement_data,
    compartments = c("Cell", "Nucleus"),
    calculation = "Mean"
  )

  expect_equal(sort(unique(seg_intensity_plot$data$channel)), "DAPI")
  expect_equal(
    sort(unique(seg_intensity_plot$data$compartment)),
    c("Cell", "Nucleus")
  )
  expect_equal(nrow(seg_intensity_plot$data), 4)
  expect_false(any(seg_intensity_plot$data$channel == "Cell"))
})
