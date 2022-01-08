---
jupyter:
  jupytext:
    formats: ipynb,R:light,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.13.5
  kernelspec:
    display_name: R - U20-v4.1
    language: R
    name: iru20v41
---

# Configure

<!-- #region tags=[] -->
## Load renv and init
<!-- #endregion -->

```R
sessionInfo()
```

```R
source("renv/activate.R")
```

```R
library(renv)
Sys.setenv(RENV_PATHS_CACHE = "/lab/solexa_weissman/cweng/Packages/R/x86_64-pc-linux-gnu-library/4.1-focal")
renv::paths$cache()
#renv::init()
```

## install my commonly used packages 

```R
.libPaths()
```

```R

```

<!-- #region tags=[] -->
### Basics and plotting
<!-- #endregion -->

```R tags=[]
#install.packages("gridExtra")
#remove.packages("devtools")
#install.packages("Matrix")
#install.packages("Matrix.utils")
#install.packages("devtools",dependencies=T)
#devtools::install_github("caleblareau/BuenColors")
#
#devtools::install_github("chenweng1991/EZsinglecell")
```

```R tags=[]
library(ggplot2)
library(gridExtra)
library(plyr)
library(dplyr)
library(Matrix)
library(Matrix.utils) 
library(BuenColors)
library(EZsinglecell)
```

<!-- #region tags=[] -->
### Single cell RNA and ATAC
<!-- #endregion -->

```R tags=[]
#renv::install("Signac")
#BiocManager::install("EnsDb.Hsapiens.v86")
```

```R
library(Seurat)
library(Signac)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
```

```R
#install.packages("Bioconductor/GenomeInfoDb@c460a6278710f3e2b31f6e67d99b17ce89609959")
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")

#BiocManager::install("EnsDb.Hsapiens.v86")
# BiocManager::install("BiocGenerics")
```



<!-- #region tags=[] -->
### Make snapshots
<!-- #endregion -->

```R
renv::dependencies()
```

```R
renv::snapshot() 
```

```R

```
