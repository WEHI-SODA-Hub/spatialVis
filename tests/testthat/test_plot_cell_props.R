test_that("plot_cell_props", {
  # load the data
  data("spe")

  # test that the output is a ggplot object
  expect_is(plot_cell_props(spe), "ggplot")
  expect_is(plot_cell_props(spe, stack = FALSE), "ggplot")
})
