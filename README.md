# RBM
RBM: a R package for microarray and RNA-Seq data analysis

Use A Resampling-Based Empirical Bayes Approach to Assess Differential Expression in Two-Color Microarrays and RNA-Seq data sets.

Bioconductor version: Release (3.9)  
Author: Dongmei Li and Chin-Yuan Liang  
Maintainer: Dongmei Li <Dongmei_Li@urmc.rochester.edu>  

[RBM on Bioconductor](https://bioconductor.org/packages/release/bioc/html/RBM.html)

**Citation** (from within R, enter citation("RBM")):
>Li D, Liang C (2019). RBM: RBM: a R package for microarray and RNA-Seq data analysis. R package version 1.16.0.

## Installation
To install this package, start R (version "3.6") and enter:
```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    
BiocManager::install("RBM")
```

For older versions of R, please refer to the appropriate Bioconductor release.

## Documentation
To view documentation for the version of this package installed in your system, start R and enter:
```R
browseVignettes("RBM")
```
