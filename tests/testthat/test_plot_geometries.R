test_that("plot_geometries() works", {
  geom_data <- get_segmentation_geometry(
    "../../inst/extdata/segmentation.geojson"
  )

  geometry_plot <- plot_geometries(geom_data)

  expect_is(geometry_plot, "ggplot")
})

test_that("plot_geometries() works for whole-cell-only geometry", {
  geojson_file <- withr::local_tempfile(fileext = ".geojson")
  write_seg_geometry_geojson(geojson_file, include_nucleus_geometry = FALSE)
  geom_data <- get_segmentation_geometry(geojson_file)

  expect_message(
    geometry_plot <- plot_geometries(
      geom_data,
      area_width = 2,
      area_height = 2,
      min_objects = 1
    ),
    "No nucleus geometry found in GeoJSON file."
  )

  expect_is(geometry_plot, "ggplot")
})

test_that("plot_geometries() works for nucleus-only geometry lists", {
  geom_data <- list(
    nucleus = list(array(c(0, 0, 2, 0, 2, 2, 0, 2, 0, 0), dim = c(5, 1, 2)))
  )

  expect_message(
    geometry_plot <- plot_geometries(
      geom_data,
      area_width = 2,
      area_height = 2,
      min_objects = 1
    ),
    "No cell geometry found in GeoJSON file."
  )

  expect_is(geometry_plot, "ggplot")
})

test_that("plot_geometries() works with nuclei_only output from parser", {
  geojson_file <- withr::local_tempfile(fileext = ".geojson")
  write_seg_geometry_geojson(geojson_file, include_nucleus_geometry = FALSE)
  geom_data <- get_segmentation_geometry(geojson_file, nuclei_only = TRUE)

  expect_message(
    geometry_plot <- plot_geometries(
      geom_data,
      area_width = 2,
      area_height = 2,
      min_objects = 1
    ),
    "No cell geometry found in GeoJSON file."
  )

  expect_is(geometry_plot, "ggplot")
})

test_that("plot_geometries() works with a background image", {
  geom_data <- get_segmentation_geometry(
    "../../inst/extdata/segmentation.geojson"
  )

  # Download tiff temporarily
  image_file <- withr::local_tempfile(fileext = ".tiff")
  download.file(
    "https://github.com/WEHI-SODA-Hub/spatialproteomics/raw/refs/heads/main/tests/data/mesmer/test_data.tiff", # nolint
    destfile = image_file,
    mode = "wb"
  )

  geometry_plot <- plot_geometries(geom_data, image = image_file)

  expect_is(geometry_plot, "ggplot")
})
