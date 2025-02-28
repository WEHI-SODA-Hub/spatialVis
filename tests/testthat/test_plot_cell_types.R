test_that("plot_cell_types", {
  # load the data
  data("spe")

  # test that the output is a ggplot object
  expect_is(plot_cell_types(spe), "ggplot")
})
