test_that("plot_geometries() works", {
  geometry_plot <- plot_geometries(
    "../../inst/extdata/segmentation.geojson"
  )

  expect_is(geometry_plot, "ggplot")
})
