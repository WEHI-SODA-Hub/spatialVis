test_that("load_hierarchies() works", {
  # Load the data
  hierarchy <- load_hierarchies("../../inst/extdata/hierarchy.yaml")

  # Check that the data is a data frame
  expect_is(hierarchy, "data.frame")

  # Check that the data has the correct number of rows
  expect_equal(nrow(hierarchy), 13)

  # Check that the data has the correct number of columns
  expect_equal(ncol(hierarchy), 4)
})
