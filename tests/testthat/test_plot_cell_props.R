test_that("plot_cell_props", {
  # load the data
  data("spe")

  # test that the output is a ggplot object
  expect_is(plot_cell_props(spe), "ggplot")
  expect_is(plot_cell_props(spe, stack = FALSE), "ggplot")
})

test_that("plot_cell_props with custom cell types", {
  # load the data
  data("spe")

  # test that the output is a ggplot object
  expect_is(plot_cell_props(spe, cell_types = c("CD4", "CD8")), "ggplot")

  # test with no stacking
  expect_is(plot_cell_props(spe, cell_types = c("CD4", "CD8"), stack = FALSE),
            "ggplot")

  # test by filtering out parent types
  expect_is(plot_cell_props(spe, cell_types = c("CD4", "CD8"),
                            parent_types = c("T_cells")), "ggplot")

  # test by excluding parent types
  expect_is(plot_cell_props(spe, cell_types = c("CD4", "CD8"),
                            exclude_parent_types = c("T_cells")), "ggplot")
})

test_that("plot_cell_props with facet_by variable", {
  # load the data
  data("spe")

  # add a custom facet variable -- disease type
  spe$diagnosis <- sample(
    c("Healthy", "Disease"),
    size = ncol(spe),
    replace = TRUE
  )

  # test that the output is a ggplot object
  expect_is(plot_cell_props(spe, facet_by = "sample_id"), "ggplot")
  expect_is(plot_cell_props(spe, facet_by = "diagnosis", stack = FALSE),
            "ggplot")
})
