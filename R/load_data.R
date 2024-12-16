library(dplyr)

#' @title Load cell hierarchies from YAML file
#'
#' @description Load cell hierarchies from a YAML file and return a data frame
#' with the hierarchy levels as columns.
#' @param hierarchy_file Path to the YAML file containing cell hierarchies
#' @return A data frame with the hierarchy levels as columns
#' @export
load_hierarchies <- function(hierarchy_file) {
  # load YAML file containing cell hierarchies
  hierarchy <- yaml::yaml.load_file(hierarchy_file) %>%
    get_hierarchy_list %>%
    make_hierarchies_table

  return(hierarchy)
}

# recursive helper function to build list obtained list containing cell
# hierarchies from terminal to top level nodes
get_hierarchy_list <- function(tree, path = NULL) {
  # return path if terminal node (first element null) is reached
  if (is.null(tree[[1]])) {
    return(list(path))
  }

  # otherwise, keep traversing the tree
  result <- list()
  if (is.list(tree)) {
    for (name in names(tree)) {
      result <- append(result, get_hierarchy_list(tree[[name]], c(name, path)))
    }
  }

  return(result)
}

# helper function to create a data frame from a list of hierarchies
make_hierarchies_table <- function(hierarchy_list) {
  # helper function to pad hierarchy to a fixed depth
  pad_hierarchy <- function(hierarchy, depth) {
    last_element <- hierarchy[length(hierarchy)]
    length(hierarchy) <- depth
    hierarchy[is.na(hierarchy)] <- last_element
    return(hierarchy)
  }

  # find depth of hierarchy
  depth <- max(purrr::map_int(hierarchy_list, length))

  # pad and reverse so higher levels come first
  padded_hierarchies <- lapply(hierarchy_list, rev) %>%
    lapply(pad_hierarchy, depth)

  # convert to data frame
  hierarchy_df <- purrr::map_dfr(padded_hierarchies,
                                 ~data.frame(t(.), stringsAsFactors = FALSE))

  # rename columns
  colnames(hierarchy_df) <- paste0("HierarchyLevel", seq(depth))

  return(hierarchy_df)
}

#' @title Create a SpatialExperiment object from expression data
#'
#' @description Takes an expression file and a data frame containing the
#' hierarchy levels and constructs a SpatialExperiment object.
#' @param expression_file Path to the expression data file
#' @param hierarchy_df Data frame containing the hierarchy levels
#' @param metadata_cols Character vector of metadata columns (default: c("Image"
#' , "Class", "In Tumour"))
#' @param centroid_x_col Column name for the x-coordinate of the cell centroid
#' (default: "Centroid X")
#' @param centroid_y_col Column name for the y-coordinate of the cell centroid
#' (default: "Centroid Y")
#' @return A SpatialExperiment object
#' @export
make_spe_from_expr_data <- function(expression_file, hierarchy_df,
                                    metadata_cols = c("Image",
                                                      "Class",
                                                      "In Tumour"),
                                    centroid_x_col = "Centroid X",
                                    centroid_y_col = "Centroid Y") {
  exp_data <- data.table::fread(expression_file, sep = ",", check.names = FALSE)

  # extract expression columns and construct matrix with markers for rows
  exp_cols <- grep(": Mean", colnames(exp_data), value = TRUE)
  exp_matrix <- exp_data %>%
    dplyr::select(dplyr::all_of(exp_cols)) %>%
    t() %>%
    as.matrix()

  rownames(exp_matrix) <- rownames(exp_matrix) %>%
    stringr::str_replace_all("(: Cell: Mean)|(: Nucleus: Mean)", "")

  colnames(exp_data)   # split marker column into constituent parts
  marker_cols <- exp_data$Class[1] %>%
    stringr::str_split(":") %>%
    unlist() %>%
    stringr::str_trim() %>%
    stringr::str_replace_all("[-+]$", "")
  hierarchy_level <- paste0("HierarchyLevel", ncol(hierarchy_df))
  marker_cols[1] <- hierarchy_level

  # extract metadata and split markers out into separate columns
  centroid_x_col <- grep(centroid_x_col, colnames(exp_data), value = TRUE)
  centroid_y_col <- grep(centroid_y_col, colnames(exp_data), value = TRUE)
  cell_metadata <- dplyr::select(exp_data, dplyr::all_of(metadata_cols)) %>%
    tidyr::separate("Class", into = marker_cols, sep = ":") %>%
    dplyr::left_join(hierarchy_df, by = hierarchy_level)
  cell_coords <- exp_data %>%
    dplyr::select(dplyr::all_of(c(centroid_x_col, centroid_y_col))) %>%
    as.matrix()

  # row_data contains the marker names
  row_data <- data.frame(rownames(exp_matrix))

  # col_data contains the metadata
  col_data <- data.frame(cell_metadata)

  spe <- SpatialExperiment::SpatialExperiment(
    assays = list(counts = exp_matrix),
    spatialCoords = cell_coords,
    colData = col_data,
    rowData = row_data
  )
  return(spe)
}
