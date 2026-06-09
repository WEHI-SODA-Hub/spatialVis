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

  has_nuclei <- grepl("Nucleus", colnames(measurement_data)) |> any()
  has_cells <- grepl("Cell", colnames(measurement_data)) |> any()
  stopifnot(has_nuclei || has_cells)

  if (!has_nuclei) {
    compartments <- c("Cell")
    message(paste("No nucleus measurements found in GeoJSON file.",
                  "Plotting whole cell measurements only."))
  } else if (!has_cells) {
    compartments <- c("Nucleus")
    message(paste("No cell measurements found in GeoJSON file.",
                  "Plotting whole nucleus measurements only."))
  }

  # Create list of measurements
  measurements <- tidyr::expand_grid(
    compartment = compartments,
    measurement = measurements
  ) |>
    dplyr::mutate(measurement_colname = paste0(compartment, ": ", measurement) # nolint
    ) |>
    dplyr::select(measurement_colname) %>% # nolint
    unlist()

  measurement_regexes <- vapply(
    measurements,
    function(measurement) {
      escaped_measurement <- measurement
      for (special_character in c("\\", ".", "|", "(", ")", "[", "]",
                                  "{", "}", "^", "$", "*", "+", "?")) {
        escaped_measurement <- gsub(
          special_character,
          paste0("\\\\", special_character),
          escaped_measurement,
          fixed = TRUE
        )
      }
      paste0("^", escaped_measurement, "(?:$|\\s)")
    },
    character(1)
  )
  col_found <- vapply(names(measurement_data), function(measurement_name) {
    any(stringr::str_detect(
      measurement_name,
      stringr::regex(measurement_regexes)
    ))
  }, logical(1))

  # Format dataframe for plotting
  df <- measurement_data %>%  # nolint
    dplyr::select(names(.)[col_found]) %>% #nolint
    tidyr::pivot_longer(dplyr::everything(),
                        names_to = "measurement",
                        values_to = "value") |>
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
