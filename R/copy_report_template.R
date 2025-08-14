#' Copy Quarto Report Template
#'
#' This function copies the Quarto (.qmd) report template from the package to
#' the current directory.
#' @param template_name Name of the template file (default:
#' "phenotyping_report_template.qmd"). Options are:
#' phenotyping_report_template.qmd or segmentation_report_template.qmd.
#' @export
copy_report_template <- function(template_name =
                                   "phenotyping_report_template.qmd",
                                 output_dir = ".", overwrite = FALSE) {

  # Check that the output directory exists
  if (!dir.exists(output_dir)) {
    stop("Output directory does not exist.")
  }

  # Get the path to the template file within the package
  template_path <- system.file(
    "templates",
    template_name,
    package = "spatialVis"
  )

  if (template_path == "") {
    stop("Template not found in the package.")
  }

  # Copy the template to the current directory
  file.copy(template_path, output_dir, overwrite = overwrite)
  message("Template copied to the current directory: ", template_name)
}
