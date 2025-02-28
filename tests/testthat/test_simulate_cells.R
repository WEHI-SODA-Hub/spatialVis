test_that("simulate_cells() returns a data frame", {
  cells <- simulate_cells()

  # test that the output is a data frame
  expect_is(cells, "data.frame")

  # test that the data frame has the correct dimensions
  expect_equal(ncol(cells), 16)
  expect_equal(nrow(cells), 5000)
})
