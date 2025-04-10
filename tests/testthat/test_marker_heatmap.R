test_that("plot_marker_heatmap", {
  # load the data
  data("spe")

  # the function suppresses the tablegrob output so we need to check the
  # character string output
  expect_true(
    grepl("TableGrob", plot_marker_heatmap(spe)[[1]])
  )

  expect_true(
    grepl("TableGrob", plot_marker_heatmap(spe, value = "proportion")[[1]])
  )
})
