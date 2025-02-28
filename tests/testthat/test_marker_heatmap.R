test_that("plot_marker_heatmap", {
  # load the data
  data("spe")

  # test that the output is a gtable object
  expect_is(plot_marker_heatmap(spe), "gtable")
})
