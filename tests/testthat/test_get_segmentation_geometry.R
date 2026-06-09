test_that("get_segmentation_geometry() works", {
  geom_data <- get_segmentation_geometry(
    "../../inst/extdata/segmentation.geojson"
  )

  expect_is(geom_data, "list")
  expect_is(geom_data$cell, "list")
  expect_is(geom_data$nucleus, "list")
})

test_that("get_segmentation_geometry() supports whole-cell-only segmentations", {
  geojson_file <- withr::local_tempfile(fileext = ".geojson")
  write_seg_geometry_geojson(geojson_file, include_nucleus_geometry = FALSE)

  expect_message(
    geom_data <- get_segmentation_geometry(geojson_file),
    "No nucleus geometry found in GeoJSON file."
  )

  expect_is(geom_data, "list")
  expect_is(geom_data$cell, "list")
  expect_false("nucleus" %in% names(geom_data))
  expect_gt(length(geom_data$cell), 0)
})

test_that("get_segmentation_geometry() supports nuclei_only mode", {
  geojson_file <- withr::local_tempfile(fileext = ".geojson")
  write_seg_geometry_geojson(geojson_file, include_nucleus_geometry = FALSE)

  geom_data <- get_segmentation_geometry(geojson_file, nuclei_only = TRUE)

  expect_is(geom_data, "list")
  expect_false("cell" %in% names(geom_data))
  expect_is(geom_data$nucleus, "list")
  expect_gt(length(geom_data$nucleus), 0)
})
