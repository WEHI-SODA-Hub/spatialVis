library(dplyr)

#' @title Plot intensity measurements for all channels from cell segmentation
#' GeoJSON file
#'
#' @description Take a GeoJSON file containing cell segmentation output and
#' plot intensity measurements for the specified compartment.
#' @param geojson_file Path to the GeoJSON file containing segmentation data
#' @param compartment Compartments to plot intensity measurements for (default:
#' c("Nucleus", "Cell", "Cytoplasm", "Membrane")). Must be at least one of the
#' following: "Cell", "Nucleus", "Cytoplasm", "Membrane".
#' @param calculation Calculation to use for intensity measurements (default:
#' "Median"). Options are "Mean", "Median", "Min", "Max" and "Std.Dev".
#' @export
#' @importFrom dplyr %>%
plot_intensity_measurements <- function(geojson_file,
                                        compartments =
                                          c("Nucleus", "Cell",
                                            "Cytoplasm", "Membrane"),
                                        calculation = "Median") {
  seg <- jsonlite::fromJSON(geojson_file)

  stopifnot(
    "features" %in% names(seg),
    "properties" %in% names(seg$features),
    "objectType" %in% names(seg$features$properties),
    "measurements" %in% names(seg$features$properties),
    all(compartments %in% c("Cell", "Nucleus", "Cytoplasm", "Membrane")),
    calculation %in% c("Mean", "Median", "Min", "Max", "Std.Dev")
  )

  properties <- seg$features$properties
  measurement_data <- properties$measurements

  # Find measurement columns for all compartments and channels
  measurements <- paste(compartments, calculation, sep = ": ")
  col_idxs <- lapply(measurements, function(measurement) {
    col_idx <- which(stringr::str_detect(colnames(measurement_data),
                                         stringr::fixed(measurement)))
    if (length(col_idx) == 0) {
      warning(paste("No measurements found for", measurement))
    }
    col_idx
  }) %>%
    unlist()

  if (length(col_idxs) == 0) {
    stop("No measurements found in data")
  }

  # Format dataframe for plotting
  df <- measurement_data[, col_idxs] %>%
    dplyr::filter(!if_all(everything(), is.na)) %>%
    tidyr::pivot_longer(everything(),
                        names_to = "measurement",
                        values_to = "intensity") %>%
    tidyr::separate(measurement, into = c("channel", "compartment", #nolint
                                          "calculation"),
                    sep = ": ", extra = "merge", fill = "right")

  pal <- colorRampPalette(
    RColorBrewer::brewer.pal(n = 11, name = "Spectral")
  )(length(unique(df$channel)))

  ggplot2::ggplot(data = df, ggplot2::aes(x = channel, y = intensity, #nolint
                                          fill = channel)) +
    ggplot2::geom_violin(adjust = 2) +
    ggplot2::geom_jitter(width = 0.1, alpha = 0.6) +
    ggplot2::facet_wrap(~ compartment, scales = "free_y") +
    ggplot2::ggtitle(paste(calculation, "intensity measurements")) +
    ggplot2::scale_fill_manual(values = pal) +
    ggplot2::theme(axis.title.y = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(angle = 90, hjust = 1),
                   axis.title.x = ggplot2::element_blank(),
                   strip.text.y = ggplot2::element_blank(),
                   panel.background = ggplot2::element_blank())
}
