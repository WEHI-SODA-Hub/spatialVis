library(dplyr)

#' @title Create spatial clusters
#'
#' @description Create spatial clusters based on the proportions of nearest
#' neighbours of cell types
#' @param cell_df A data frame containing cell data
#' @param image_col Column name in cell_df that contains image names
#' @param phenotype_col Column name in cell_df that contains cell type
#' information
#' @param centroid_x_col Column name in cell_df that contains x-coordinates
#' @param centroid_y_col Column name in cell_df that contains y-coordinates
#' @param k_closest The number of nearest neighbours to consider
#' @param max_dist The maximum distance to consider a neighbour
#' @param number_clusters_kmeans The number of clusters to create with kmeans
#' @param keep_proportions Whether to keep the proportions of nearest neighbours
#' @param seed The random seed to use
#' @return A data frame with spatial clusters
#' @export
#' @importFrom dplyr %>%
create_spatial_clusters <- function(spe,
                                    image_col = "Image",
                                    phenotype_col = "HierarchyLevel4",
                                    centroid_x_col = "Cell.X.Position",
                                    centroid_y_col = "Cell.Y.Position",
                                    k_closest = 10,
                                    max_dist = 50,
                                    number_clusters_kmeans = 10,
                                    seed = 123) {

  # Get the spatial coordinates and grep out whole column names
  coords <- SpatialExperiment::spatialCoords(spe) %>% as.data.frame()
  centroid_x_col <- grep(centroid_x_col, colnames(coords), value = TRUE)
  centroid_y_col <- grep(centroid_y_col, colnames(coords), value = TRUE)

  # Add the spatial coordinates to the cell data frame
  cell_df <- SingleCellExperiment::colData(spe) %>%
    as.data.frame() %>%
    dplyr::mutate(!!as.name(centroid_x_col) := coords[[centroid_x_col]], # nolint: object_usage_linter, line_length_linter.
                  !!as.name(centroid_y_col) := coords[[centroid_y_col]])

  proportions_df <- find_nn_proportions(cell_df,
                                        image_col,
                                        phenotype_col,
                                        centroid_x_col,
                                        centroid_y_col,
                                        k_closest,
                                        max_dist)

  set.seed(seed)
  cluster_objs <- kmeans(proportions_df %>%
                           dplyr::select(-!!dplyr::sym(image_col)),
                         number_clusters_kmeans)
  spe$cluster <- cluster_objs$cluster

  return(spe)
}

# Find proportions of nearest neighbours
find_nn_proportions <- function(cell_df,
                                image_col,
                                phenotype_col,
                                centroid_x_col,
                                centroid_y_col,
                                k_closest = 10,
                                max_dist = 50) {
  cell_df$temp_cell_id <- rownames(cell_df)

  proportions_df <- cell_df %>%
    dplyr::group_by(!!as.name(image_col)) %>%
    dplyr::group_modify(~ {
      calc_nn_proportions(cell_df = .x,
                          phenotype_col = phenotype_col,
                          centroid_x_col = centroid_x_col,
                          centroid_y_col = centroid_y_col,
                          k_closest = k_closest,
                          max_dist = max_dist)
    }, .keep = TRUE) %>%
    dplyr::ungroup()


  proportions_df <- proportions_df[match(cell_df$temp_cell_id,
                                         proportions_df$temp_cell_id), ]
  proportions_df <- proportions_df %>% replace(is.na(.), 0) # nolint: object_usage_linter, line_length_linter.
  proportions_df <- proportions_df %>% subset(select = -c(temp_cell_id)) # nolint: object_usage_linter, line_length_linter.

  return(proportions_df)
}

# Calculate the proportions of nearest neighbours
calc_nn_proportions <- function(cell_df,
                                phenotype_col,
                                centroid_x_col,
                                centroid_y_col,
                                k_closest = 10,
                                max_dist = 50) {
  if (nrow(cell_df) == 1) {
    return(cell_df[, c("temp_cell_id")])
  }

  # Adjust k_closest to include the cell itself
  k_closest <- k_closest + 1
  if (k_closest > nrow(cell_df)) {
    warning(paste0(cell_df$Image[1],
                   " has less rows than the set nearest neighbours. NN will be calculated with nrow()")) # nolint: line_length_linter.
    k_closest <- nrow(cell_df)
  }


  cell_df <- data.frame(cell_df, check.names = FALSE)

  # Find nearest neighbors using RANN package
  closest <- RANN::nn2(
    data = cell_df[, c(centroid_x_col, centroid_y_col)],
    query = cell_df[, c(centroid_x_col, centroid_y_col)],
    k = k_closest
  )

  # Data frame of: rows = query cells and columns = top k closest cell indexes
  dists <- data.frame(closest$nn.dists)
  colnames(dists) <- paste("Dist", 1:k_closest - 1, sep = ":")

  # Data frame of index of the closest cells
  nn_idx <- data.frame(closest$nn.idx)
  colnames(nn_idx) <- paste("NN:idx", 1:k_closest - 1, sep = ":")

  # Remove the first column since it is the cell itself
  dists <- dists %>% dplyr::select(-1)
  nn_idx <- nn_idx %>% dplyr::select(-1)

  # Initialise empty dataframe to store the string nearest neighbours
  nn <- nn_idx
  colnames(nn) <- paste("NN", 1:(k_closest - 1), sep = ":")
  nn[] <- cell_df[match(as.integer(unlist(nn_idx)),
                        as.integer(rownames(cell_df))), phenotype_col]

  # Rename cell type if too far away
  null_cell_type <- "TOO_FAR"

  # Remove neighbours which are too far away
  nn[dists > max_dist] <- null_cell_type

  # calculate the proportions. For each cell, count the number of a certain cell
  # type being its neighbour.
  nn$temp_id <- as.numeric(rownames(nn))
  proportions <- nn %>%
    tidyr::pivot_longer(cols = -temp_id, names_to = "Hierarchy:Level", # nolint: object_usage_linter, line_length_linter.
                        values_to = "Phenotype") %>%
    dplyr::group_by(temp_id) %>%
    dplyr::count(Phenotype) %>% # nolint: object_usage_linter.
    tidyr::pivot_wider(names_from = Phenotype, values_from = n) %>%
    replace(is.na(.), 0) # nolint: object_usage_linter.
  proportions <- proportions %>% subset(select = -c(temp_id)) # nolint: object_usage_linter, line_length_linter.

  proportions$temp_cell_id <- cell_df$temp_cell_id


  f <- function(x) x / (k_closest - 1)

  if ("TOO_FAR" %in% colnames(proportions)) {
    proportions <- proportions %>% subset(select = -c(TOO_FAR)) # nolint: object_usage_linter, line_length_linter.
  }
  proportions <- proportions %>%
    dplyr::mutate_if(is.numeric, f)

  return(proportions)
}
