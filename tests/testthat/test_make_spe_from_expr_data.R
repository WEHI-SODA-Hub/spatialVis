test_that("make_spe_from_expr_data() works", {
  # Load the data
  spe <- make_spe_from_expr_data("../../data-raw/simulated.csv",
                                 "../../inst/extdata/hierarchy.yaml")
  col_data <- SingleCellExperiment::colData(spe)
  spatial_coords <- SpatialExperiment::spatialCoords(spe)

  # Check that the data is a data frame
  expect_is(spe, "SpatialExperiment")

  # Check that the data has the correct number of rows
  expect_equal(nrow(col_data), 5000)

  # Check that the data has the correct number of columns
  expect_equal(ncol(col_data), 18)

  # Check that the coord data has the correct number of rows
  expect_equal(nrow(spatial_coords), 5000)

  # Check that the coord data has the correct number of columns
  expect_equal(ncol(spatial_coords), 2)
})
