test_that("get_segmentation_geometry() works", {
  geom_data <- get_segmentation_geometry(
    "../../inst/extdata/segmentation.geojson"
  )

  expect_is(geom_data, "list")
  expect_is(geom_data$cell, "list")
  expect_is(geom_data$nucleus, "list")
})
