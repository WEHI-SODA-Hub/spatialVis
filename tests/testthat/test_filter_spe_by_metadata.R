test_that("filter_spe_by_metadata() works", {
  data(spe)

  # Test filtering by metadata
  filtered_spe <- filter_spe_by_metadata(
    spe,
    "HierarchyLevel1",
    c("Immune_cells", "Stromal_cells")
  )

  # Orthogonal filtering
  filtered_df <- SingleCellExperiment::colData(spe) %>%
    as.data.frame() %>%
    dplyr::filter(HierarchyLevel1 %in% c("Immune_cells", "Stromal_cells"))

  expect_equal(
    nrow(SingleCellExperiment::colData(filtered_spe)),
    nrow(filtered_df)
  )

  expect_equal(
    nrow(SpatialExperiment::spatialCoords(filtered_spe)),
    nrow(filtered_df)
  )

  expect_equal(
    length(colnames(filtered_spe)),
    nrow(filtered_df)
  )
})
