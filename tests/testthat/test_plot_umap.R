test_that("plot_umap", {
  # load the data
  data("spe")

  # test that the output is a ggplot object
  expect_is(plot_umap(spe), "ggplot")
})
