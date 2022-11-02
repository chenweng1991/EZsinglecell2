
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
## Data

```r 
data(ATACWhite) 
data(RNAWhite) 
data(Griffin_Signatures)
```

## Hem dateset convinience

Standard cell types
```r
STD.CellType.lv<-c("CD4","CD8","NK","B","Plasma","ProB","CLP","LMPP","pDC","MPP","HSC","CMP","MKP","MEP","EryP","cDC","MDP","CDP","GMP","Mono")
```

Plot for known progenitor markers  and BMMC markers
```r
## A function to plot progenitor marker genes
PlotProgMarker<-function(ob){
options(repr.plot.width=20, repr.plot.height=8)
DefaultAssay(ob)<-"SCT"
p1<-FeaturePlot(ob,"CD34",reduction = 'wnn.umap') ## HSPC
p2<-FeaturePlot(ob,"HLF",reduction = 'wnn.umap')  ## HSC
p3<-FeaturePlot(ob,"CRHBP",reduction = 'wnn.umap')  ## HSC
p4<-FeaturePlot(ob,"MPO",reduction = 'wnn.umap')    ## GMP
p5<-FeaturePlot(ob,"PLEK",reduction = 'wnn.umap')   ## MkP
p6<-FeaturePlot(ob,"ARPP21",reduction = 'wnn.umap')  ## Lym
p7<-FeaturePlot(ob,"GATA1",reduction = 'wnn.umap')   ## MEP
p8<-FeaturePlot(ob,"SPI1",reduction = 'wnn.umap')   ## GMP
grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,ncol=5)
}
## example to plot BMMC markers
FeaturePlot(Donor01_BMMC_Multiome_wrapper.filtered, reduction = 'wnn.umap',features = c("CCR7","CD8A","GNLY","MS4A1","LYZ","SLC4A1"),ncol=6)
```

Add and plot the signatures
```r
Donor01_BMMC_Multiome_wrapper.filtered<-AddHemSignature(Donor01_BMMC_Multiome_wrapper.filtered)
FeaturePlot(object = Donor01_BMMC_Multiome_wrapper.filtered, features = c("Sig.cDC1","Sig.Prog1","Sig.EarlyE1","Sig.LateE1","Sig.ProMono1","Sig.Mono1","Sig.ncMono1","Sig.cDC1","Sig.pDC1","Sig.ProB1","Sig.PreB1","Sig.B1","Sig.Plasma1","Sig.T1","Sig.CTL1","Sig.NK1"),reduction = 'wnn.umap',pt.size =1 )

## To add category
AddSTD_Cat<-function(ob=Donor01_BMMC_Multiome_wrapper.filtered){
ob@meta.data$STD_Cat<-recode(ob@meta.data$STD.CellType,
HSC="Stem",
MPP="EarlyP",CMP="EarlyP",MKP="EarlyP",
MEP="Mye_P",GMP="Mye_P",MDP="Mye_P",CDP="Mye_P", EryP="Mye_P",
LMPP="Lym_P",CLP="Lym_P",ProB="Lym_P",Plasma="Lym_P",
Mono="Mye",Ery="Mye",mDC="Mye",
CD4="Lym",CD8="Lym" ,NK="Lym"  ,B="Lym" ,pDC="Lym") 
ob@meta.data$STD_Cat<-factor(ob@meta.data$STD_Cat,levels=c("Stem","EarlyP","Mye_P","Lym_P","Mye","Lym"))
return(ob)
}
                             
AddSTD_Cat2<-function(ob){
ob@meta.data$STD_Cat2<-recode(ob@meta.data$STD.CellType,
HSC="Stem",
MPP="EarlyP",CMP="EarlyP",MKP="EarlyP",
MEP="LateP",GMP="LateP",MDP="LateP",CDP="LateP",
LMPP="LateP",CLP="LateP",ProB="LateP",EryP="LateP",
Mono="mature",Ery="mature",mDC="mature",
CD4="mature",CD8="mature" ,NK="mature",B="mature" ,pDC="mature",Plasma="mature") 
ob@meta.data$STD_Cat2<-factor(ob@meta.data$STD_Cat2,levels=c("Stem","EarlyP","LateP","mature"))
return(ob)
}    
```

Subset and reclustering
```r
Cluster4.ob<- subset(x = Donor01_BMMC_Multiome_wrapper.filtered,subset= seurat_clusters==4)
Cluster4.ob<-Reclustering(Cluster4.ob)
```

Find and show markers
```r
DefaultAssay(Donor01_BMMC_Multiome_wrapper.filtered) <- "SCT"
DN1_BMMC.markers <- FindAllMarkers(Donor01_BMMC_Multiome_wrapper.filtered, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
DN1_BMMC.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)-> top5
subset(top5,cluster==10)
```

Change Ident
```r
Donor01_BMMC_Multiome_wrapper.filtered<- SetIdent(Donor01_BMMC_Multiome_wrapper.filtered, value = "seurat_clusters")
DimPlot(Donor01_BMMC_Multiome_wrapper.filtered, reduction = "wnn.umap", label = T, label.size = 5, repel = TRUE) + ggtitle("WNN")
```
