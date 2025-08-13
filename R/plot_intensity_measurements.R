library(dplyr)

#' @title Plot intensity measurements for all channels from cell segmentation
#' GeoJSON file
#'
#' @description Take a GeoJSON file containing cell segmentation output and
#' plot intensity measurements for the specified compartment.
#' @param geojson_file Path to the GeoJSON file containing segmentation data
#' @param compartment Compartment to plot intensity measurements for (default:
#' "Cell"). Options are "Cell", "Nucleus", "Cytoplasm", and "Membrane".
#' @param calculation Calculation to use for intensity measurements (default:
#' "Median"). Options are "Mean", "Median", "Min", "Max" and "Std.Dev".
#' @export
#' @importFrom dplyr %>%
plot_intensity_measurements <- function(geojson_file,
                                        compartment = "Cell",
                                        calculation = "Median") {
  seg <- jsonlite::fromJSON(geojson_file)

  stopifnot(
    "features" %in% names(seg),
    "properties" %in% names(seg$features),
    "objectType" %in% names(seg$features$properties),
    "measurements" %in% names(seg$features$properties),
    compartment %in% c("Cell", "Nucleus", "Cytoplasm", "Membrane"),
    calculation %in% c("Mean", "Median", "Min", "Max", "Std.Dev")
  )

  properties <- seg$features$properties
  measurement_data <- properties$measurements

  measurement <- paste(compartment, calculation, sep = ": ")
  col_idx <- which(stringr::str_detect(colnames(measurement_data),
                                       stringr::fixed(measurement)))

  if (length(col_idx) == 0) {
    stop(paste("Measurement", measurement, "not found in data"))
  }

  channels <- colnames(measurement_data)[col_idx] %>%
    stringr::str_split(":") %>%
    lapply(., first) %>% #nolint
    unlist() %>%
    stringr::str_trim()

  if (length(channels) == 0) {
    stop(paste("No channels found for measurement", measurement))
  }

  df <- measurement_data[, col_idx] %>%
    dplyr::filter(!if_all(everything(), is.na)) %>%
    tidyr::pivot_longer(everything(),
                        names_to = "channel",
                        values_to = "intensity") %>%
    dplyr::mutate(channel = stringr::str_split(channel, ":") %>% #nolint
                      lapply(., first) %>% #nolint
                      unlist() %>%
                      stringr::str_trim())

  pal <- colorRampPalette(
    RColorBrewer::brewer.pal(n = 11, name = "Spectral")
  )(length(channels))

  ggplot2::ggplot(data = df, ggplot2::aes(x = channel, y = intensity, #nolint
                                          fill = channel)) +
    ggplot2::geom_violin(adjust = 2) +
    ggplot2::geom_jitter(width = 0.1, alpha = 0.2) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
    ggplot2::ggtitle(paste(compartment, calculation,
                           "intensity measurements")) +
    ggplot2::scale_fill_manual(values = pal) +
    ggplot2::theme(axis.title.y = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(angle = 90, hjust = 1),
                   axis.title.x = ggplot2::element_blank(),
                   strip.text.y = ggplot2::element_blank(),
                   panel.background = ggplot2::element_blank())
}
