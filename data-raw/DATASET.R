## code to prepare `DATASET` dataset goes here
ATACWhite<-read.table("/lab/solexa_weissman/cweng/Genomes/10X/WhiteList_10X_Multiome.ATAC")
RNAWhite<-read.table("/lab/solexa_weissman/cweng/Genomes/10X/WhiteList_10X_Multiome.RNA")



usethis::use_data(ATACWhite, RNAWhite, overwrite = TRUE)
