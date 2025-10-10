test_that("plot_geometries() works", {
  geom_data <- get_segmentation_geometry(
    "../../inst/extdata/segmentation.geojson"
  )

  geometry_plot <- plot_geometries(geom_data)

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
