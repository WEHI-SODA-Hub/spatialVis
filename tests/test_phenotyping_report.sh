#!/usr/bin/env bash -C -e -u -o pipefail

Rscript -e "devtools::install_local('.', force = TRUE, upgrade = 'never')"
Rscript -e "spatialVis::copy_report_template(
    template_name = 'phenotyping_report_template.qmd',
    output_dir = '.',
    overwrite = TRUE
)"
sample=Test
quarto render phenotyping_report_template.qmd \
        --to html \
        --no-cache \
        --output phenotyping_report.html \
        -P hierarchy_file:inst/extdata/hierarchy.yaml \
        -P expression_file:data-raw/simulated.csv \
        -P sample_name:$sample \
        -P save_rdata:false
