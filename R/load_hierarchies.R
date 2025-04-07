library(dplyr)

#' @title Load cell hierarchies from YAML file
#'
#' @description Load cell hierarchies from a YAML file and return a data frame
#' with the hierarchy levels as columns.
#' @param hierarchy_file Path to the YAML file containing cell hierarchies
#' @return A data frame with the hierarchy levels as columns
#' @export
#' @importFrom dplyr %>%
load_hierarchies <- function(hierarchy_file) {
  # load YAML file containing cell hierarchies
  hierarchy <- yaml::yaml.load_file(hierarchy_file) %>%
    get_hierarchy_list %>%
    make_hierarchies_table

  hierarchy
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

  result
}

# helper function to create a data frame from a list of hierarchies
make_hierarchies_table <- function(hierarchy_list) {
  # helper function to pad hierarchy to a fixed depth
  pad_hierarchy <- function(hierarchy, depth) {
    last_element <- hierarchy[length(hierarchy)]
    length(hierarchy) <- depth
    hierarchy[is.na(hierarchy)] <- last_element
    hierarchy
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

  hierarchy_df
}
