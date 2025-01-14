# spatialVis

SpatialVis is an R package for downstream QC, visualisation and analysis of MIBI
data. Currently, the following analyses are supported:

- marker heatmaps (expression and proportion positive)
- UMAP
- spatial plotting of cell types

The package can also simulate simple spatial data for testing.

The package loads data exported from QuPath post-phenotyping and constructs a 
`SpatialExperiment` object. A file containing the cell type hierarchy is also
required. An example of these files can be foundunder the `data` directory.

## Installation

The package can be installed as follows (tested on R version 4.4.1).

Dependencies:

- [R](https://www.r-project.org/)
- [ImageMagick](https://imagemagick.org/)

Installation:

```R
if (!requireNamespace("devtools", quietly = TRUE)) {
      install.packages("devtools")}
if (!requireNamespace("magick", quietly = TRUE)) {
      install.packages("magick")}
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("SpatialExperiment", "spaSim"), ask = FALSE)
devtools::install_github("WEHI-SODA-Hub/spatialVis")
```

## Basic usage

Loading your data:

```R
library(spatialVis)

hierarchy_df <- load_hierarchies("data/hierarchy.yaml")
spe <- make_spe_from_expr_data("data/simulated.csv", hierarchy_df)
```

The `spe` object can now be used to generate plots and analyses.

Marker heatmap:

```R
plot_marker_heatmap(spe)
```

Marker positivity proportion heatmap:

```R
plot_marker_heatmap(spe, value = "proportion")
```

UMAP:

```R
plot_umap(spe)
```

Spatial plot of cells:

```R
plot_celltypes(spe)
```

All of these functions have configurable parameters. Check the function help for
more info.
