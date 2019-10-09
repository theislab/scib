# LISI
Methods to compute Local Inverse Simpson's Index (LISI)

Check out how we use LISI to measure single cell integration methods in the Harmony paper: [bioRxiv](https://www.biorxiv.org/content/early/2018/11/04/461954)

# System Requirements 

LISI has been tested on R versions >= 3.4. Please consult the DESCRIPTION file for more details on required R packages. LISI has been tested on Linux, OS X, and Windows platforms.

# Installation 

Install LISI with devtools: 

```
install.packages('devtools')
devtools::install_github("immunogenomics/lisi")
```

Installation may include compiling C++ code from source, so it can take a few minutes.

# Demo

With a matrix of PCA embeddings and meta data table, compute LISI with respect to one or more cell labels. For demonstration, we include a matrix of cell coordinates `X` and a matrix of cell meta data `meta_data`. 

```
library(lisi)
lisi_res <- lisi::compute_lisi(lisi::X, lisi::meta_data, c('label1', 'label2'))
head(lisi_res)
```

For more information about the `compute_lisi` function, use `?compute_lisi` in R. 

