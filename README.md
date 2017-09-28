ezSeurat: Interactive single-cell data analysis via the Seurat pipeline
============

**NOTE**: <u>This is the development version of ezSeurat</u>

### ezSeurat

User friendly GUI implementation of Seurat analysis pipeline, from input of expression sparse matrices to DEG analysis.
This package was designed to allow analysis to be done with easier changing of visualisation parameters.

### Installation

As this is an R package, it is assumed you have R installed. However, do ensure that you're using the latest version.
It is also recommended to run R via RStudio.

To install ezSeurat:

``` r
library("devtools")
install_github("MattMyint/ezSeurat")
```

### Usage

After successfully installing ezSeurat, run the following:

``` r
library("ezSeurat")
ezSeurat()
```

### Acknowledgements

Many thanks to Satija Lab for developing Seurat

Rahul Satija, Andrew Butler and Paul Hoffman (2017). Seurat: Tools for Single Cell Genomics. R package
  version 2.0.1. https://CRAN.R-project.org/package=Seurat


