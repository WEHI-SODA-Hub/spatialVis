test_that("plot_geometries() works", {
  geom_data <- get_segmentation_geometry(
    "../../inst/extdata/segmentation.geojson"
  )

  geometry_plot <- plot_geometries(geom_data)

  expect_is(geometry_plot, "ggplot")
})
