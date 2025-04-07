library(dplyr)

#' @title Create a SpatialExperiment object from expression data
#'
#' @description Takes an expression file and a data frame containing the
#' hierarchy levels and constructs a SpatialExperiment object.
#' @param expression_file Path to the expression data file
#' @param hierarchy_file Path to the YAML file containing cell hierarchies
#' @param metadata_cols Character vector of metadata columns (default: c("Image"
#' , "Class", "In Tumour"))
#' @param marker_col Column name for the marker class (default: "Class")
#' @param centroid_x_col Column name for the x-coordinate of the cell centroid
#' (default: "Centroid X")
#' @param centroid_y_col Column name for the y-coordinate of the cell centroid
#' (default: "Centroid Y")
#' @param are_markers_split Logical indicating if the marker column is already
#' split into multiple columns (default: FALSE)
#' @return A SpatialExperiment object
#' @export
#' @importFrom dplyr %>%
make_spe_from_expr_data <- function(expression_file, hierarchy_file,
                                    metadata_cols = c("Image", "In Tumour"),
                                    marker_col = "Class",
                                    centroid_x_col = "Centroid X",
                                    centroid_y_col = "Centroid Y",
                                    are_markers_split = FALSE) {
  hierarchy_df <- spatialVis::load_hierarchies(hierarchy_file)
  exp_data <- data.table::fread(expression_file, sep = ",", check.names = FALSE)

  stopifnot(
    all(marker_col %in% colnames(exp_data)),
    all(metadata_cols %in% colnames(exp_data)),
    nrow(hierarchy_df) > 0
  )

  # extract expression columns and construct matrix with markers for rows
  exp_cols <- grep(": Mean", colnames(exp_data), value = TRUE)
  exp_matrix <- exp_data %>%
    dplyr::select(dplyr::all_of(exp_cols)) %>%
    t() %>%
    as.matrix()

  rownames(exp_matrix) <- rownames(exp_matrix) %>%
    stringr::str_replace_all("(: Cell: Mean)|(: Nucleus: Mean)", "") %>%
    stringr::str_replace_all("[\\-\\+_() ]", "_")

  hierarchy_level <- paste0("HierarchyLevel", ncol(hierarchy_df))
  if (are_markers_split) {
    # if markers are already split, set marker cols to everything after the
    # marker col
    marker_idx <- which(colnames(exp_data) == marker_col)
    marker_cols <- colnames(exp_data)[marker_idx:ncol(exp_data)]

    stopifnot(length(marker_cols) > 0)
    marker_cols[1] <- hierarchy_level

    cell_metadata <- exp_data %>%
      dplyr::rename(!!hierarchy_level := !!rlang::sym(marker_col)) %>% # nolint: object_usage_linter, line_length_linter.
      dplyr::select(dplyr::all_of(c(marker_cols, metadata_cols)))
  } else {
    # split marker column into constituent parts
    marker_cols <- exp_data[[marker_col]][1] %>%
      stringr::str_split(":") %>%
      unlist() %>%
      stringr::str_trim() %>%
      stringr::str_replace_all("[-+]$", "") %>%
      stringr::str_replace_all("[\\-\\+_() ]", "_")
    marker_cols[1] <- hierarchy_level
    cell_metadata <- dplyr::select(exp_data,
                                   dplyr::all_of(c(marker_col,
                                                   metadata_cols))) %>%
      tidyr::separate(marker_col, into = marker_cols, sep = ":")
  }

  # extract metadata and split markers out into separate columns
  cell_metadata <- cell_metadata %>%
    dplyr::mutate(dplyr::across(dplyr::any_of(marker_cols[-1]),
                                ~ stringr::str_sub(., -1))) %>%
    dplyr::mutate(Cell_ID = paste0("Cell", seq(1, nrow(.)))) %>% # nolint: object_usage_linter, line_length_linter.
    dplyr::left_join(hierarchy_df, by = hierarchy_level) %>%
    dplyr::mutate(dplyr::across(dplyr::starts_with("HierarchyLevel"),
                                ~ stringr::str_replace_all(., " ", "_")))

  centroid_x_col <- grep(centroid_x_col, colnames(exp_data), value = TRUE)
  centroid_y_col <- grep(centroid_y_col, colnames(exp_data), value = TRUE)
  cell_coords <- exp_data %>%
    dplyr::select(dplyr::all_of(c(centroid_x_col, centroid_y_col))) %>%
    dplyr::rename(`Cell.Y.Position` = !!as.name(centroid_y_col), # nolint: object_usage_linter, line_length_linter.
                  `Cell.X.Position` = !!as.name(centroid_x_col)) %>% # nolint: object_usage_linter, line_length_linter.
    as.matrix()

  # row_data contains the marker names
  row_data <- data.frame(rownames(exp_matrix))

  # col_data contains the metadata
  col_data <- data.frame(cell_metadata, row.names = cell_metadata$Cell_ID)

  spe <- SpatialExperiment::SpatialExperiment(
    assays = list(counts = exp_matrix),
    spatialCoords = cell_coords,
    colData = col_data,
    rowData = row_data
  )
  spe
}
