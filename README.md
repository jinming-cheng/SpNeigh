
<!-- README.md is generated from README.Rmd. Please edit that file -->

# SpNeigh <img src="man/figures/SpNeigh_logo.png" align="right" alt="" width="120" />

<!-- badges: start -->

[![R-CMD-check](https://github.com/jinming-cheng/SpNeigh/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/jinming-cheng/SpNeigh/actions)
[![Codecov test
coverage](https://codecov.io/gh/jinming-cheng/SpNeigh/graph/badge.svg)](https://app.codecov.io/gh/jinming-cheng/SpNeigh)
<!-- badges: end -->

SpNeigh provides methods for neighborhood-aware analysis of spatial
transcriptomics data. It supports boundary detection, spatial weighting
(centroid- and boundary-based), spatially informed differential
expression using spline-based models, and spatial enrichment analysis
via the Spatial Enrichment Index (SEI). Designed for compatibility with
Seurat objects, SpatialExperiment objects and spatial data frames,
SpNeigh enables interpretable, publication-ready analysis of spatial
gene expression patterns.

Quick start guide can be found
[here](https://jinming-cheng.github.io/SpNeigh/articles/SpNeigh.html).

## Installation

The *SpNeigh* package can be installed from GitHub by using:

``` r
devtools::install_github("jinming-cheng/SpNeigh")
```

## Citation

Please cite this article if you use SpNeigh:

    To cite SpNeigh in publications, please use:

      Cheng J, Chow P, Liu N (2025). "SpNeigh: spatial neighborhood and
      differential expression analysis for high-resolution spatial
      transcriptomics." _bioRxiv_. doi:10.1101/2025.11.07.687304
      <https://doi.org/10.1101/2025.11.07.687304>.

    A BibTeX entry for LaTeX users is

      @Article{,
        title = {SpNeigh: spatial neighborhood and differential expression analysis for high-resolution spatial transcriptomics},
        author = {Jinming Cheng and Pierce Kah Hoe Chow and Nan Liu},
        journal = {bioRxiv},
        year = {2025},
        doi = {10.1101/2025.11.07.687304},
      }
