test_that("make_spe_from_expr_data() works", {
  # Load the data
  spe <- make_spe_from_expr_data("../../data-raw/simulated.csv",
                                 "../../inst/extdata/hierarchy.yaml")
  col_data <- SingleCellExperiment::colData(spe)
  spatial_coords <- SpatialExperiment::spatialCoords(spe)

  # Check that the data is a data frame
  expect_is(spe, "SpatialExperiment")

  # Check that the data has the correct number of rows
  expect_equal(nrow(col_data), 5000)

  # Check that the data has the correct number of columns
  expect_equal(ncol(col_data), 18)

  # Check that the coord data has the correct number of rows
  expect_equal(nrow(spatial_coords), 5000)

  # Check that the coord data has the correct number of columns
  expect_equal(ncol(spatial_coords), 2)
})

test_that("make_spe_from_expr_data() works for split marker input", {
  # split the marker columns
  simu <- read.delim("../../data-raw/simulated.csv", sep = ",", header = TRUE)

  # Move Class column to the end
  simu <- simu %>%
    dplyr::select(-Class, dplyr::everything(), Class)

  # Get marker columns
  marker_cols <- simu[["Class"]][1] %>%
    stringr::str_split(":") %>%
    unlist() %>%
    stringr::str_trim() %>%
    stringr::str_replace_all("[-+]$", "") %>%
    stringr::str_replace_all("[\\-\\+_() ]", "_")
  marker_cols[1] <- "CellType"

  # Split the Class column into separate marker cols
  simu <- simu %>%
    tidyr::separate("Class", into = marker_cols, sep = ":")

  # Replace dots with spaces in column names to mimic the original data
  colnames(simu) <- gsub("\\.\\.Mean", ": Mean", colnames(simu))
  colnames(simu) <- gsub("\\.", " ", colnames(simu))

  # Write the data to a temporary file
  tempfile <- tempfile(fileext = ".csv")
  write.csv(simu, tempfile, row.names = FALSE)

  # Load the data
  spe <- make_spe_from_expr_data(tempfile,
                                 "../../inst/extdata/hierarchy.yaml",
                                 marker_col = "CellType",
                                 are_markers_split = TRUE)
  col_data <- SingleCellExperiment::colData(spe)
  spatial_coords <- SpatialExperiment::spatialCoords(spe)

  # Check that the data is a data frame
  expect_is(spe, "SpatialExperiment")

  # Check that the data has the correct number of rows
  expect_equal(nrow(col_data), 5000)

  # Check that the data has the correct number of columns
  expect_equal(ncol(col_data), 18)

  # Check that the coord data has the correct number of rows
  expect_equal(nrow(spatial_coords), 5000)

  # Check that the coord data has the correct number of columns
  expect_equal(ncol(spatial_coords), 2)

})
