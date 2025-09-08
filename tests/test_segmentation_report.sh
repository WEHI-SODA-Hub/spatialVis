#!/usr/bin/env bash -C -e -u -o pipefail

Rscript -e "devtools::install_local('.', force = TRUE, upgrade = 'never')"
Rscript -e "spatialVis::copy_report_template(
    template_name = 'segmentation_report_template.qmd',
    output_dir = '.',
    overwrite = TRUE
)"
sample=Test
quarto render segmentation_report_template.qmd \
    --to html \
    --no-cache \
    --output cellpose_report.html \
    -P sample_name:$sample \
    -P geojson_file:inst/extdata/segmentation.geojson \
    -P nuclear_channel:"Channel 1" \
    -P membrane_channels:"Channel 2" \
    -P run_cellpose:true \
    -P cellpose_diameter:30 \
    -P cellpose_min_area:200 \
    -P cellpose_flow_threshold:0.4 \
    -P cellpose_cellprob_threshold:0.0 \
    -P cellpose_model_type:null \
    -P cellpose_pretrained_model:null
    
quarto render segmentation_report_template.qmd \
    --to html \
    --no-cache \
    --output mesmer_report.html \
    -P sample_name:$sample \
    -P geojson_file:inst/extdata/segmentation.geojson \
    -P nuclear_channel:"Channel 1" \
    -P membrane_channels:"Channel 2" \
    -P run_mesmer:true \
    -P mesmer_segmentation_level:null \
    -P mesmer_maxima_threshold:0.1 \
    -P mesmer_interior_threshold:0.3 \
    -P mesmer_maxima_smooth:0 \
    -P mesmer_min_nuclei_area:40 \
    -P mesmer_remove_border_cells:true \
    -P mesmer_pixel_expansion:3 \
    -P mesmer_padding:0
