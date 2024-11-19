library(yaml)
library(data.table)
library(dplyr)
library(tidyverse)
library(tidyr)
library(purrr)
library(stringr)
library(data.table)
library(ggplot2)
library(SpatialExperiment)

load_hierarchies <- function(hierarchy_file) {
  # load YAML file containing cell hierarchies
  hierarchy <- yaml.load_file(hierarchy_file) %>%
    get_hierarchy_list %>%
    make_hierarchies_table

  return(hierarchy)
}

get_hierarchy_list <- function(tree, path = NULL) {
  # build list of hierarchies obtained from YAML file containing
  # cell hierarchies from terminal to top level nodes

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

make_hierarchies_table <- function(hierarchy_list) {

  pad_hierarchy <- function(hierarchy, depth) {
    length(hierarchy) <- depth
    hierarchy[is.na(hierarchy)] <- hierarchy[1]
    return(hierarchy)
  }

  # find depth of hierarchy
  depth <- max(map_int(hierarchy_list, length))

  # pad and reverse so higher levels come first
  padded_hierarchies <- lapply(hierarchy_list, pad_hierarchy, depth) %>%
    lapply(rev)

  # convert to data frame
  hierarchy_df <- map_dfr(padded_hierarchies,
                          ~data.frame(t(.), stringsAsFactors = FALSE))

  # rename columns
  colnames(hierarchy_df) <- paste0("HierarchyLevel", seq(depth))

  return(hierarchy_df)
}

make_spe_from_expr_data <- function(expression_file, hierarchy_df) {
  exp_data <- fread(expression_file, sep = ",", check.names = FALSE)

  # extract expression columns and construct matrix with markers for rows
  exp_cols <- grep(": Mean", colnames(exp_data), value = TRUE)
  exp_matrix <- exp_data %>% select(all_of(exp_cols)) %>% t() %>% as.matrix()

  # split marker column into constituent parts
  marker_cols <- exp_data$Class[1] %>%
    str_split(":") %>%
    unlist() %>%
    str_trim() %>%
    str_replace_all("[-+]$", "")
  hierarchy_level <- paste0("HierarchyLevel", ncol(hierarchy_df))
  marker_cols[1] <- hierarchy_level

  # extract metadata and split markers out into separate columns
  # TODO: parametrise column names
  metadata_cols <- c("Image", "Class", "In Tumour")
  cell_metadata <- select(exp_data, all_of(metadata_cols)) %>%
    separate("Class", into = marker_cols, sep = ":") %>%
    left_join(hierarchy_df, by = hierarchy_level)
  cell_coords <- exp_data %>%
    select(c("Centroid X um", "Centroid Y um")) %>%
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
