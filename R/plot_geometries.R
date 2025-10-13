#' @title Plot cell and nucleus geometries from segmentation GeoJSON file
#'
#' @description This function takes cell and nucleus geometries from a list and
#' creates a ggplot showing the outlines of the cells and nuclei within a
#' randomly selected area containing at least a specified minimum number of
#' cells.
#'
#' @param geom_data A list containing cell and nucleus geometries, as
#' returned by `get_segmentation_geometry()`
#' @param area_width Width of the area to plot (default: 500)
#' @param area_height Height of the area to plot (default: 500)
#' @param min_objects Minimum number of cell objects required in the area
#' (default: 10)
#' @param image Optional raster image to use as background (default: NULL). The
#' image should have at least 2 layers: nuclear and membrane in that order.
#'
#' @return A ggplot object containing the cell and nucleus geometries
#' @export
#' @importFrom dplyr %>%
plot_geometries <- function(geom_data, area_width = 500,
                            area_height = 500, min_objects = 10,
                            image = NULL) {
  stopifnot(
    "cell" %in% names(geom_data),
    "nucleus" %in% names(geom_data)
  )

  bbox <- get_random_bbox(geom_data, area_width, area_height, min_objects)
  cell_ids_in_area <- get_objects_in_area(geom_data, bbox, type = "cell")

  cells <- geom_data$cell[cell_ids_in_area]
  nuclei <- geom_data$nucleus[cell_ids_in_area]

  plot_cell_geometries(cells, nuclei, bbox, image)
}

# Function to plot cell and nucleus geometries
plot_cell_geometries <- function(cells, nuclei, bbox, image = NULL) {
  plot_data <- data.frame()

  # Helper function to convert polygon coordinates to data frame
  polygon_to_df <- function(points, id, type = "cell") {
    data.frame(
      x = points[, , 1],
      y = points[, , 2],
      id = id,
      type = type
    )
  }

  for (i in seq_along(cells)) {
    points <- cells[[i]]
    cell_df <- polygon_to_df(points, paste0("cell_", i), "cell")
    plot_data <- rbind(plot_data, cell_df)
  }

  for (i in seq_along(nuclei)) {
    points <- nuclei[[i]]
    nuc_df <- polygon_to_df(points, paste0("nucleus_", i), "nucleus")
    plot_data <- rbind(plot_data, nuc_df)
  }

  # Get plot limits which are likely larger than the bounding box
  xmin <- min(c(plot_data$x, bbox$xmin))
  xmax <- max(c(plot_data$x, bbox$xmax))
  ymin <- min(c(plot_data$y, bbox$ymin))
  ymax <- max(c(plot_data$y, bbox$ymax))

  if (!is.null(image)) {
    img <- terra::rast(image) |> suppressWarnings()

    n_layers <- terra::nlyr(img)
    if (n_layers < 2) {
      warning("Image must have at least 2 layers to plot nuclear and membrane ",
              "channels. Skipping background image.")
      image <- NULL
    } else {
      nuc <- terra::crop(img[[1]], terra::ext(xmin, xmax, ymin, ymax)) |>
        terra::flip("vertical")
      mem <- terra::crop(img[[2]], terra::ext(xmin, xmax, ymin, ymax)) |>
        terra::flip("vertical")

      # Normalize values to 0-1 range
      nuc_minmax <- terra::minmax(nuc)
      mem_minmax <- terra::minmax(mem)
      nuc_norm <- (nuc - nuc_minmax[1]) / diff(nuc_minmax)
      mem_norm <- (mem - mem_minmax[1]) / diff(mem_minmax)

      # Apply gamma correction for better visibility
      gamma <- 0.5
      nuc_norm <- nuc_norm ^ gamma
      mem_norm <- mem_norm ^ gamma

      # Create RGB values with nuclei and membrane in green and blue channels
      raster_df <- expand.grid(
        x = seq(round(xmin), round(xmax), length.out = ncol(nuc)),
        y = seq(round(ymin), round(ymax), length.out = nrow(nuc))
      )
      raster_df$red <- 0
      raster_df$green <- as.vector(nuc_norm)
      raster_df$blue <- as.vector(mem_norm)

      # Use geom_raster with rgb()
      raster_df$colour <- grDevices::rgb(
        raster_df$red, raster_df$green, raster_df$blue
      )
    }
  }

  # Colours for plotting when no image is provided
  cell_colour <- "#4daf4a"
  nuc_colour <- "#377eb8"
  plot_title <- paste0("Objects overlapping area\n",
                       "X: ", round(bbox$xmin), "-", round(bbox$xmax), ", ",
                       "Y: ", round(bbox$ymin), "-", round(bbox$ymax))

  geom_plot <- ggplot2::ggplot(plot_data,
                  ggplot2::aes(x = x, y = y, group = id, colour = type)) # nolint
  # If we have a background image, add it to the plot under the geometries
  if (!is.null(image)) {
    # Tweak colours for visibility against image
    cell_colour <- "#00BFFF"
    nuc_colour <- "#FFFF00"

    geom_plot <- geom_plot +
      ggplot2::geom_raster(data = raster_df,
                           ggplot2::aes(x = x, y = y, fill = I(colour)), # nolint
                           inherit.aes = FALSE, alpha = 0.8)
  } else {
    geom_plot <- geom_plot + ggplot2::scale_fill_gradient(
      low = "white", high = "black", guide = "none"
    )
  }
  geom_plot +
    ggplot2::geom_polygon(fill = NA, linewidth = 0.5) +
    ggplot2::annotate("rect",
                      xmin = bbox$xmin, xmax = bbox$xmax,
                      ymin = bbox$ymin, ymax = bbox$ymax,
                      alpha = 0.3, fill = NA, linewidth = 1,
                      colour = "pink") +
    ggplot2::scale_color_manual(
      values = c("cell" = cell_colour, "nucleus" = nuc_colour)
    ) +
    ggplot2::coord_equal() +
    ggplot2::scale_y_reverse() +
    ggplot2::theme(strip.text.y = ggplot2::element_blank(),
                   panel.background = ggplot2::element_blank(),
                   legend.position = "bottom") +
    ggplot2::labs(title = plot_title, x = "X coordinate", y = "Y coordinate")
}

# Get bounding box from coordinates
get_bbox <- function(coords) {
  x_coords <- coords[, , 1]
  y_coords <- coords[, , 2]
  list(
    xmin = min(x_coords),
    xmax = max(x_coords),
    ymin = min(y_coords),
    ymax = max(y_coords)
  )
}

# Select random bounding box with enough objects
get_random_bbox <- function(geom_data,
                            area_width = 500,
                            area_height = 500,
                            min_objects = 10) {
  all_bboxes <- lapply(geom_data$cell, function(coords) {
    get_bbox(coords)
  })

  # Find global bounds
  global_xmin <- min(sapply(all_bboxes, function(bb) bb$xmin))
  global_xmax <- max(sapply(all_bboxes, function(bb) bb$xmax))
  global_ymin <- min(sapply(all_bboxes, function(bb) bb$ymin))
  global_ymax <- max(sapply(all_bboxes, function(bb) bb$ymax))

  # Try random areas until we find one with enough objects
  max_attempts <- 100
  attempt <- 1

  while (attempt <= max_attempts) {
    # Random starting point
    start_x <- runif(1, global_xmin, global_xmax - area_width)
    start_y <- runif(1, global_ymin, global_ymax - area_height)

    area_bbox <- list(
      xmin = start_x,
      xmax = start_x + area_width,
      ymin = start_y,
      ymax = start_y + area_height
    )

    # Count objects that overlap with this area
    overlapping_objects <- 0
    for (i in seq_along(all_bboxes)) {
      bb <- all_bboxes[[i]]
      if (bb$xmax >= area_bbox$xmin &&
            bb$xmin <= area_bbox$xmax &&
            bb$ymax >= area_bbox$ymin &&
            bb$ymin <= area_bbox$ymax) {
        overlapping_objects <- overlapping_objects + 1
      }
    }

    if (overlapping_objects >= min_objects) {
      return(area_bbox)
    }

    attempt <- attempt + 1
  }

  stop("Could not find area with sufficient objects after ",
       max_attempts, " attempts")
}

# Get object IDs within area bounding box
get_objects_in_area <- function(geom_data, area_bbox, type = "cell") {
  object_ids <- c()
  coords <- list()

  if (type == "cell") {
    coords <- geom_data$cell
  } else if (type == "nucleus") {
    coords <- geom_data$nucleus
  } else {
    stop("Invalid type: ", type)
  }

  for (i in seq_along(coords)) {
    points <- coords[[i]]
    bb <- get_bbox(points)

    # Check if object overlaps with area
    if (bb$xmax >= area_bbox$xmin &&
          bb$xmin <= area_bbox$xmax &&
          bb$ymax >= area_bbox$ymin &&
          bb$ymin <= area_bbox$ymax) {
      object_ids <- c(object_ids, i)
    }
  }
  object_ids
}
