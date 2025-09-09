library(dplyr)

#' @title Plot intensity measurements for all channels from cell segmentation
#' GeoJSON file
#'
#' @description Take segmentation measurements and plot intensity measurements
#' for the specified compartments.
#'
#' @param geojson_file Path to the GeoJSON file containing segmentation data
#' @param compartment Compartments to plot intensity measurements for (default:
#' c("Nucleus", "Cell", "Cytoplasm", "Membrane")). Must be at least one of the
#' following: "Cell", "Nucleus", "Cytoplasm", "Membrane".
#' @param calculation Calculation to use for intensity measurements (default:
#' "Median"). Options are "Mean", "Median", "Min", "Max" and "Std.Dev".
#' @param percentiles Numeric vector of length 2 specifying the lower and upper
#' percentiles to use for x-axis limits (default: c(0.05, 0.95)). Must be
#' between 0 and 1, and the first value must be less than the second value.
#'
#' @export
#' @importFrom dplyr %>%
plot_intensity_measurements <- function(measurement_data,
                                        compartments =
                                          c("Nucleus", "Cell",
                                            "Cytoplasm", "Membrane"),
                                        calculation = "Median",
                                        percentiles = c(0.05, 0.95)) {
  stopifnot(
    all(compartments %in% c("Cell", "Nucleus", "Cytoplasm", "Membrane")),
    calculation %in% c("Mean", "Median", "Min", "Max", "Std.Dev"),
    length(percentiles) == 2,
    all(percentiles >= 0 & percentiles <= 1),
    percentiles[1] < percentiles[2]
  )

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
  df <- measurement_data[, col_idxs] %>% # nolint
    dplyr::filter(!dplyr::if_all(dplyr::everything(), is.na)) %>%
    tidyr::pivot_longer(dplyr::everything(),
                        names_to = "measurement",
                        values_to = "intensity") %>%
    tidyr::separate(measurement, into = c("channel", "compartment", #nolint
                                          "calculation"),
                    sep = ": ", extra = "merge", fill = "right")

  # Calculate percentile to set x-axis limits
  lower_limit <- quantile(df$intensity, percentiles[1], na.rm = TRUE)
  upper_limit <- quantile(df$intensity, percentiles[2], na.rm = TRUE)

  # Palette based on number of channels
  pal <- colorRampPalette(
    RColorBrewer::brewer.pal(n = 11, name = "Spectral")
  )(length(compartments))

  ggplot2::ggplot(data = df, ggplot2::aes(x = intensity, fill = compartment)) + #nolint
    ggplot2::geom_density(alpha = 0.4) +
    ggplot2::facet_wrap(~ channel, scales = "free_y") +
    ggplot2::ggtitle(paste(calculation, "intensity measurements")) +
    ggplot2::scale_fill_manual(values = pal) +
    ggplot2::xlim(lower_limit, upper_limit) +
    ggplot2::theme(axis.title.y = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(angle = 90, hjust = 1),
                   axis.title.x = ggplot2::element_blank(),
                   strip.text.y = ggplot2::element_blank(),
                   panel.background = ggplot2::element_blank(),
                   legend.position = "bottom")
}
