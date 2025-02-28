test_that("create_spatial_clusters() works", {
  # Load the data
  data("spe")

  # Create spatial clusters
  spe <- create_spatial_clusters(spe)

  # Check that the data is a data frame
  expect_type(spe$cluster, "integer")

  # Check that the data has the correct number of values
  expect_equal(length(spe$cluster), 5000)
})
