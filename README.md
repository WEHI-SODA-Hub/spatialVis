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

spe <- make_spe_from_expr_data("expression.csv", "hierarchy.yaml")
```

See `raw-data/simulated.csv` for an example of the expression data format, and
`inst/extdata/hierarchy.yaml` for an example of a hierarchy file.

You can test the package by loading a the simulated data set:

```R
data(spe)
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
plot_cell_types(spe)
```

Plot proportion of cells:

```R
plot_cell_props(spe)
```

Create spatial clusters and plot cell_type proportions by cluster:

```R
spe <- create_spatial_clusters(spe)
plot_cluster_cell_props(spe)
```

All of these functions have configurable parameters. Check the function help for
more info.
