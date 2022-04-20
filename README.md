# SUITOR: selecting the number of mutational signatures through cross-validation
<br/>

### Introduction
For the  _de novo_ mutational signature analysis, 
estimating the correct number of signatures is the crucial starting point, 
since it influences all the downstream steps, including extraction of signature profiles,
 estimation of signature activities and classification of tumors based on the andestimated activities.
 Here we present an **R** package `SUITOR`, an unsupervised cross-validation tool to select 
the optimal number of signatures. 
The two main functions are `suitor()` and `suitorExtractWH()`.

For more information please refer to the package user maunual and vignette.
<br/>

### Installation
To install from Bioconductor:
```r
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager") 
  BiocManager::install("SUITOR") 
```
