library(dplyr)

#' @title Plot compartment measurements
#'
#' @description This function takes segmentation measurements and creates violin
#' plots for each compartment.
#'
#' @param measurement_data A data.table containing segmentation measurements
#' @param compartments Character vector of compartments to plot (default:
#' c("Cell", "Nucleus"))
#' @param measurements Character vector of measurements to plot for each
#' compartment (default: c("Area", "Length", "Circularity", "Solidity",
#' "Max diameter", "Min diameter"))
#'
#' @return A ggplot object containing the violin plots for each compartment
#' @export
#' @importFrom dplyr %>%
plot_compartment_measurements <- function(measurement_data,
                                          compartments = c("Cell", "Nucleus"),
                                          measurements =
                                            c("Area", "Length",
                                              "Circularity",
                                              "Solidity",
                                              "Max diameter",
                                              "Min diameter")) {
  stopifnot(
    all(compartments %in% c("Cell", "Nucleus")),
    all(measurements %in% c("Area", "Length", "Circularity", "Solidity",
                            "Max diameter", "Min diameter"))
  )

  # Create list of measurements
  measurements <- tidyr::expand_grid(
    compartment = compartments,
    measurement = measurements
  ) %>%
    dplyr::mutate(measurement_colname = paste0(compartment, ": ", measurement) # nolint
    ) %>%
    dplyr::select(measurement_colname) %>% # nolint
    unlist()

  # Extract all measurements above from the data
  col_found <- sapply(names(measurement_data),
                      function(x) any(startsWith(x, measurements)))

  # Format dataframe for plotting
  df <- measurement_data %>%
    dplyr::select(names(.)[col_found]) %>% #nolint
    tidyr::pivot_longer(dplyr::everything(),
                        names_to = "measurement",
                        values_to = "value") %>%
    tidyr::separate(measurement, into = c("compartment", "measurement"), #nolint
                    sep = ": ", extra = "merge", fill = "right")

  pal <- colorRampPalette(
    RColorBrewer::brewer.pal(n = 11, name = "Spectral")
  )(length(compartments))

  ggplot2::ggplot(data = df, ggplot2::aes(x = compartment, y = value, #nolint
                                          fill = compartment)) + #nolint
    ggplot2::geom_violin(adjust = 2, alpha = 0.8) +
    ggplot2::geom_jitter(width = 0.1, alpha = 0.1, size = 0.2) +
    ggplot2::facet_wrap(~ measurement, scales = "free_y") +
    ggplot2::ggtitle("Compartment measurements") +
    ggplot2::scale_fill_manual(values = pal) +
    ggplot2::theme(axis.title.y = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(angle = 90, hjust = 1),
                   axis.title.x = ggplot2::element_blank(),
                   strip.text.y = ggplot2::element_blank(),
                   panel.background = ggplot2::element_blank())
}
