test_that("plot_cluster_props_by_meta() works", {
  # Load the data
  data("spe")

  # Create spatial clusters
  spe <- create_spatial_clusters(spe)

  # test that the output is a ggplot object
  expect_is(plot_cluster_props_by_meta(spe), "ggplot")
})
