library(dplyr)

#' @title Summarise segmentation statistics from GeoJSON file
#'
#' @description Take segmentation measurements and summarise statistics,
#' returning an output table. If measurements are not available, the function
#' will return only the cell counts.
#'
#' @param geojson_file Path to the GeoJSON file containing segmentation data
#' @param compartments Character vector of compartments to summarise (default:
#' c("Cell", "Nucleus"))
#' @param compartment_measurements Character vector of measurements to summarise
#' for each compartment (default: c("Area", "Length", "Circularity", "Solidity",
#' "Max diameter", "Min diameter"))
#' @param other_measurements Character vector of additional measurements to
#' summarise (default: c("Nucleus/Cell area ratio"))
#' @param summary_funcs List of functions to summarise the measurements
#' (default: list(Mean = mean, Median = median))
#'
#' @return A data.frame object
#' @export
#' @importFrom dplyr %>%
summarise_segmentation_stats <- function(measurement_data,
                                         compartments = c("Cell", "Nucleus"),
                                         compartment_measurements =
                                           c("Area", "Length",
                                             "Circularity",
                                             "Solidity",
                                             "Max diameter",
                                             "Min diameter"),
                                         other_measurements =
                                           c("Nucleus/Cell area ratio"),
                                         summary_funcs =
                                           list(Mean = mean, Median = median)) {
  # Create list of measurements
  measurements <- tidyr::expand_grid(
    compartment = compartments,
    measurement = compartment_measurements
  ) %>%
    dplyr::mutate(measurement_colname = paste0(compartment, ": ", measurement) # nolint
    ) %>%
    dplyr::select(measurement_colname) %>% # nolint
    unlist()
  measurements <- c(measurements, other_measurements) %>% as.character()

  # Calculate stats on data
  summary_list <- list()
  for (measurement in measurements) {
    col_idx <- which(stringr::str_detect(colnames(measurement_data),
                                         stringr::fixed(measurement)))
    if (length(col_idx) == 1) {
      tmp <- measurement_data[[col_idx]]
      summary_list[[measurement]] <- sapply(
        summary_funcs, function(func) {
          func(tmp[!is.na(tmp)])
        }
      )
    } else {
      warning(paste("Measurement", measurement, "not found in data"))
      summary_list[[measurement]] <- rep(NA_real_,
                                         length(summary_funcs))
    }
  }

  summary_df <- do.call(rbind, summary_list) %>%
    as.data.frame() %>%
    dplyr::mutate(Measurement = rownames(.)) %>% # nolint
    dplyr::select(Measurement, dplyr::everything()) # nolint
  rownames(summary_df) <- NULL

  summary_df
}
