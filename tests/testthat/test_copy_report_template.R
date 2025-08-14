test_that("copy_report_template() works for phenotyping report", {
  # Create a temporary directory
  temp_dir <- tempdir()

  # Copy the report template to the temporary directory
  template_name <- "phenotyping_report_template.qmd"
  copy_report_template(
    template_name = template_name,
    output_dir = temp_dir
  )

  # Check that the report template was copied successfully
  expect_true(file.exists(file.path(temp_dir, template_name)))
})

test_that("copy_report_template() works for segmentation report", {
  # Create a temporary directory
  temp_dir <- tempdir()

  # Copy the report template to the temporary directory
  template_name <- "segmentation_report_template.qmd"
  copy_report_template(
    template_name = template_name,
    output_dir = temp_dir
  )

  # Check that the report template was copied successfully
  expect_true(file.exists(file.path(temp_dir, template_name)))
})
