# BBKNN algorithm

## bbknn version 0.1.0

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/bbknn)](https://cran.r-project.org/package=bbknn)
[![Travis Build Status](https://travis-ci.org/TomKellyGenetics/bbknn.svg?branch=master)](https://travis-ci.org/TomKellyGenetics/bbknn)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/TomKellyGenetics/bbknn?branch=master&svg=true)](https://ci.appveyor.com/project/TomKellyGenetics/bbknn)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![codecov](https://codecov.io/gh/TomKellyGenetics/bbknn/branch/master/graph/badge.svg)](https://codecov.io/gh/TomKellyGenetics/bbknn)

## Batch effect correction with the BBKNN algorithm in R

This package allows calling the BBKNN algorthm for batch effect correction from R. See the Python implementation for more details: 

https://github.com/Teichlab/bbknn

## Install

This package requires the 'bbknn', 'anndata', and 'scanpy' modules for python (3) to be installed on your system. For example:

``pip3 install bbknn anndata scanpy.api``

If you do not have root access, you can use `pip3 install --user` or `pip3 install --prefix` to install these in your user directory (which you have write permissions for) and ensure that this directory is in your PATH so that Python can find it.

The 'devtools' package will be used to install 'bbknn' and the dependancies (igraph and reticulate):

```R
if (!requireNamespace("devtools"))
    install.packages("devtools")
devtools::install_github("TomKellyGenetics/bbknn")
```

## Usage

This package provides a function to perform clustering with the BBKNN algorithm:

```R
data_matrix <- matrix(rnorm(1000), 20, 50)
batches <- c(rep(1, 20), rep(2, 20), rep(3, 10))
corrected_matrix <- bbknn(data_matrix, batches)
```

This matrix can then be used for plotting with tSNE or UMAP and further analysis of clusters.

### Citation

Please cite this implementation R in if you use it:

```
To cite the bbknn package in publications use:

  S. Thomas Kelly (2018). bbknn: R implementation of the BBKNN algorithm. R
  package version 0.1.0 https://github.com/TomKellyGenetics/bbknn

A BibTeX entry for LaTeX users is

  @Manual{,
    title = {bbknn: R implementation of the BBKNN algorithm},
    author = {S. Thomas Kelly},
    year = {2018},
    note = {R package version 0.1.0},
    url = {https://github.com/TomKellyGenetics/bbknn},
  }
 ```

Please also cite the original publication of this algorithm.

```
Park, Jong-Eun and Polanski, Krzysztof and Meyer, Kerstin and Teichmann, Sarah A (2018) Fast Batch Alignment of Single Cell Transcriptomes Unifies Multiple Mouse Cell Atlases into an Integrated Landscape. `bioRxiv:397042 <https://doi.org/10.1101/397042>
```
