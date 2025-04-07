test_that("plot_cells_spatially", {
  # load the data
  data("spe")

  # test that the output is a ggplot object
  expect_is(plot_cells_spatially(spe), "ggplot")

  # test with cluster colouring
  spe <- create_spatial_clusters(spe)
  expect_is(plot_cells_spatially(spe, colour_by = "cluster"), "ggplot")

  # test with tumour colouring
  expect_is(plot_cells_spatially(spe, colour_by = "In.Tumour"), "ggplot")
})
