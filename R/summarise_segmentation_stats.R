library(dplyr)

.esc_meas_regex <- function(x) {
  escaped <- x
  regex_special_characters <- c("\\", ".", "|", "(", ")", "[", "]",
                                "{", "}", "^", "$", "*", "+", "?")

  for (special_character in regex_special_characters) {
    escaped <- gsub(
      special_character,
      paste0("\\\\", special_character),
      escaped,
      fixed = TRUE
    )
  }

  escaped
}

.find_meas_col <- function(measurement_names, measurement) {
  measurement_regex <- paste0(
    "^",
    .esc_meas_regex(measurement),
    "(?:$|\\s)"
  )

  which(
    stringr::str_detect(
      measurement_names,
      stringr::regex(measurement_regex)
    )
  )
}

#' @title Summarise segmentation statistics from GeoJSON file
#'
#' @description Take segmentation measurements and summarise statistics,
#' returning an output table. If measurements are not available, the function
#' will return only the cell counts.
#'
#' @param compartments Character vector of compartments to summarise (default:
#' c("Cell", "Nucleus"))
#' @param compartment_measurements Character vector of measurements to summarise
#' for each compartment (default: c("Area", "Length", "Circularity", "Solidity",
#' "Max diameter", "Min diameter"))
#' @param other_measurements Character vector of additional measurements to
#' summarise (default: c("Nucleus/Cell area ratio"))
#' @param summary_funcs List of functions to summarise the measurements
#' (default: list(Mean = mean, Median = median))
#' @param keep_extra_measurements Logical indicating whether to summarise
#' additional percentile and erosion/expansion-bin measurements when present
#' (default: FALSE). Neighbour-derived measurements are always excluded.
#' @details Percentile and erosion/expansion-bin measurements are excluded by
#' default.
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
                                           list(Mean = mean, Median = median),
                                         keep_extra_measurements = FALSE) {
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
    col_idx <- .find_meas_col(colnames(measurement_data), measurement)
    if (length(col_idx) == 1) {
      tmp <- measurement_data[[col_idx]]
      measurement_colname <- colnames(measurement_data)[col_idx]
      summary_list[[measurement_colname]] <- sapply(
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

  if (keep_extra_measurements) {
    extra_measurements <- names(measurement_data)[
      vapply(
        names(measurement_data),
        .is_extra_col,
        logical(1),
        compartments = compartments
      )
    ]
    extra_measurements <- extra_measurements[
      !startsWith(extra_measurements, "Neighbours: ")
    ]

    for (measurement in extra_measurements) {
      tmp <- measurement_data[[measurement]]
      summary_list[[measurement]] <- sapply(
        summary_funcs, function(func) {
          func(tmp[!is.na(tmp)])
        }
      )
    }
  }

  summary_df <- do.call(rbind, summary_list) %>%
    as.data.frame() %>%
    dplyr::mutate(Measurement = rownames(.)) %>% # nolint
    dplyr::select(Measurement, dplyr::everything()) # nolint
  rownames(summary_df) <- NULL

  summary_df
}
