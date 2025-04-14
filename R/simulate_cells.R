library(dplyr)

# List of example markers
example_markers <- c("CD3", "CD8", "CD4", "CD20", "CD68",
                     "PD1", "PDL1", "Ki67", "FoxP3", "Cytokeratin")
#' @title Simulate cells
#' @description Simulate cells for testing package functions
#' @param markers A character vector of marker names
#' @param n_cells The number of cells to simulate
#' @param width The width of the image
#' @param height The height of the image
#' @param seed The random seed to use
#' @return A data frame with simulated cell data
#' @export
#' @importFrom dplyr %>%
simulate_cells <- function(markers = example_markers, n_cells = 5000,
                           width = 2000, height = 2000, seed = 123,
                           plot_image = FALSE) {
  set.seed(seed)

  # use spaSim to simulate background cells
  bg <- spaSim::simulate_background_cells(n_cells = n_cells,
                                          width = width,
                                          height = height,
                                          method = "Hardcore",
                                          min_d = 10,
                                          oversampling_rate = 1.6,
                                          Cell.Type = "Other",
                                          plot_image = FALSE)

  idents <- c("Stromal cells", "Epithelial cells", "Other")
  mix_bg <- spaSim::simulate_mixing(bg_sample = bg,
                                    idents = idents,
                                    props = c(0.2, 0.3, 0.5),
                                    plot_image = FALSE)

  cluster1 <- list(name_of_cluster_cell = "Treg cells", size = 500,
                   shape = "Oval", centre_loc = data.frame(x = 600, y = 600),
                   infiltration_types = c("Treg cells", "Other"),
                   infiltration_proportions = c(0.1, 0.05))
  cluster2 <- list(name_of_cluster_cell = "CD4 T cells", size = 600,
                   shape = "Irregular",
                   centre_loc = data.frame(x = 1500, y = 500),
                   infiltration_types = c("CD4 T cells", "Other"),
                   infiltration_proportions = c(0.1, 0.05))
  cluster_properties <- list(cluster1, cluster2)
  clusters <- spaSim::simulate_clusters(bg_sample = mix_bg,
                                        n_clusters = 2,
                                        bg_type = "Other",
                                        cluster_properties = cluster_properties,
                                        plot_image = plot_image)

  # create empty expression values
  expr <- matrix(nrow = n_cells, ncol = length(markers), 0)
  colnames(expr) <- markers

  # create a data frame with expression values for each marker
  marker_data <- cbind(clusters, expr) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(dplyr::across(dplyr::all_of(markers),
                                ~ simulate_expression(1, mean_log = 0.5,
                                                      sd_log = 0.5))) %>%
    dplyr::ungroup()

  # add randomly positive marker to each cell
  marker_data$Class <- sapply(marker_data$Cell.Type,
    function(x) {
      paste(x, get_random_positive_marker(markers), sep = ": ")
    }
  ) %>%
    as.character()

  # 'fake' static values
  marker_data$Image <- "image.tiff"
  marker_data$In_Tumour <- 0

  # rename cols
  marker_data <- marker_data %>%
    dplyr::rename(`Centroid Y` = Cell.Y.Position, # nolint: object_usage_linter.
                  `Centroid X` = Cell.X.Position) %>% # nolint: object_usage_linter, line_length_linter.
    dplyr::rename_with(~ paste0(.x, ": Cell: Mean"), dplyr::all_of(markers))

  marker_data
}

# Return one positive marker, all others are negative in format:
# CD3-: CD8-: CD4+ ... etc.
get_random_positive_marker <- function(markers) {
  positive_marker <- sample(markers, 1)
  result <- paste0(
    sapply(markers, function(marker) {
      if (marker == positive_marker) {
        paste0(marker, "+:")
      } else {
        paste0(marker, "-:")
      }
    }),
    collapse = " "
  )
  gsub(":$", "", result)
}

# Function to simulate marker expression with a log-normal distribution
simulate_expression <- function(n, mean_log = 0, sd_log = 1) {
  rlnorm(n, meanlog = mean_log, sdlog = sd_log)
}
