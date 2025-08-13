test_that("summarise_segmentation_stats() works", {
  # Load the data
  seg_stats <- summarise_segmentation_stats(
    "../../inst/extdata/test_segmentations.geojson"
  )

  # Check that the data is a data frame
  expect_is(seg_stats, "data.frame")

  # Check that the data has the correct number of rows
  expect_equal(nrow(seg_stats), 14)
})
