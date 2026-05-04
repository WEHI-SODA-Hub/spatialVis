.is_extra_col <- function(measurement_name,
                          compartments = c("Cell", "Nucleus",
                                           "Cytoplasm", "Membrane")) {
  if (identical(measurement_name, "objectType")) {
    return(FALSE)
  }
  if (startsWith(measurement_name, "Neighbours: ")) {
    return(TRUE)
  }

  measurement_parts <- strsplit(measurement_name, ": ", fixed = TRUE)[[1]]
  has_compartment <- any(measurement_parts %in% compartments)

  is_percentile_measurement <- has_compartment &&
    "Percentile" %in% measurement_parts

  is_bin_measurement <- length(measurement_parts) == 3 &&
    measurement_parts[1] %in% compartments &&
    stringr::str_detect(measurement_parts[2],
                        "^(ErosionBin|ExpansionBin)_[0-9]+$") &&
    stringr::str_detect(
      measurement_parts[3],
      stringr::regex("^(Area_px|Area_Fraction|Depth_px)$",
                     ignore_case = TRUE)
    )

  is_percentile_measurement || is_bin_measurement
}
