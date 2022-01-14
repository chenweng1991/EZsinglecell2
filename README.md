
# EZsinglecell2

<!-- badges: start -->
<!-- badges: end -->

The goal of EZsinglecell2 is to simply the rotine single cell analysis
The testing jupyter notebook is here:
lab/solexa_weissman/cweng/Packages/MyMultiome/Helpers/EZsinglecell2/EZsinglecell2Running.ipynb

## Installation

You can install the development version of EZsinglecell2 like so:

``` r
# Install development version from GitHub:
# install.packages("devtools")
devtools::install_github("chenweng1991/EZsinglecell2")
library(EZsinglecell2)
```

## Example


Wrap Seurat RNA clustering

``` r
bmmc.data <- Read10X(data.dir = "/lab/solexa_weissman/cweng/Projects/MitoTracing_Velocity/SecondaryAnalysis/Donor01_BMMC_1/CellRanger/Donor01_BMMC_1/outs/filtered_feature_bc_matrix")
bmmc.ob<-GEM_Wrapper(mtx=bmmc.data$`Gene Expression`,exp="DN1_BMMC1",res=0.5)
```


Wrap Seurat Multiomics clustering

``` r
Donor01_CD34_1_Multiome_wrapper<-Multi_Wrapper(path="/lab/solexa_weissman/cweng/Projects/MitoTracing_Velocity/SecondaryAnalysis/Donor01_CD34_1_Multiomekit/CellRanger/Donor01_CD34_1/outs")
```

Seurat Multiomics clustering with subsetting
``` r
Donor01_CD34_1_Multiome_wrapper<-readRDS("/lab/solexa_weissman/cweng/Projects/RDS/Donor01_CD34_1_Multiome_wrapper.RDS")
HSC.meta<-subset(Donor01_CD34_1_Multiome_wrapper$seurat@meta.data,CellType=="HSC")
Donor01_HSC_Multiome_wrapper<-Multi_Wrapper(path="/lab/solexa_weissman/cweng/Projects/MitoTracing_Velocity/SecondaryAnalysis/Donor01_CD34_1_Multiomekit/CellRanger/Donor01_CD34_1/outs",atacmin=1000,umimin =1000,CellID = row.names(HSC.meta) )
```


Wrap Seurat ATAC clustering

``` r
# Need generate peak ~ cell sparse Matrix, if starting with long sparse matrix
library(Matrix.utils)
setwd('/lab/solexa_weissman/cweng/Projects/MitoTracing_Velocity/SecondaryAnalysis/Donor01_BMMC_1')

# Read and generate ATAC_RNA meta table
ATAC_RNA.summary<-readRDS("/lab/solexa_weissman/cweng/Projects/RDS/Donor1_bmmc01_ATAC_RNA.summary.RDS")
ATAC_RNA.summary.filtered<-subset(ATAC_RNA.summary,nCount_RNA>800 & nCount_RNA<25000 & UniqFragment>1500)

# build the dgCMatrix
FragmentsOnPeak<-read.table("COMBINE.3_Multiome_Donor01_BMMC_ATAC_Nova.uniqmapped.Peak_bc_sparse_mtx")
FragmentsOnPeak.filtered<-subset(FragmentsOnPeak,V2 %in% ATAC_RNA.summary.filtered$ATACBC)
names(FragmentsOnPeak.filtered)<-c("Peaks","Cells","Cts")
PeakVSCell.filtered.Mtx<-dMcast(FragmentsOnPeak.filtered,Peaks~Cells)

# Run ATAC clustering
bmmc.filtered.atac<-ATAC_Wrapper(PeakVSCell.filtered.Mtx)
```
