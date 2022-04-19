## code to prepare `DATASET` dataset goes here
ATACWhite<-read.table("/lab/solexa_weissman/cweng/Genomes/10X/WhiteList_10X_Multiome.ATAC")
RNAWhite<-read.table("/lab/solexa_weissman/cweng/Genomes/10X/WhiteList_10X_Multiome.RNA")

## A gift list from Peter Van Galen
Griffin_Signatures<-read.csv("Griffin_Signatures.csv")

usethis::use_data(ATACWhite, RNAWhite, overwrite = TRUE)
usethis::use_data(Griffin_Signatures, overwrite = TRUE)