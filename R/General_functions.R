#' Tomerge_v2
#'
#' This function is to quickly merge two dataframe by rownames, but can choose to leave A or B all information
#' @param A  dataframe A
#' @param B  dataframe B
#' @return  return a data frame with merged information
#' @export
#' @examples
#' Tomerge_v2(A,B)

Tomerge_v2<-function(A,B,leavex=T,leavey=F){
	mergeAB<-merge(A,B,by="row.names",all.x=leavex,all.y=leavey)
	row.names(mergeAB)<-mergeAB[,1]
	mergeAB<-mergeAB[,-1]
	return(mergeAB)
}



#' Function to translate the RNA barcode into ATAC barcode and add a column
#'
#' This function allows you to input the metadata with row name as cell barcode
#' @param meta  a dataframe with the row names as the RNA cell barcode usually with the post -1
#' @param RNAclusterPost  Usually it is -1(default), but can be changed accordingly
#' @return p from ggplot2
#' @examples
#' Translate_RNA2ATAC(meta)
#' @export
Translate_RNA2ATAC<-function(meta=bmmc.filtered@meta.data,RNAclusterPost="-1"){
data(ATACWhite)
data(RNAWhite)
# ATACWhite<-read.table("/lab/solexa_weissman/cweng/Genomes/10X/WhiteList_10X_Multiome.ATAC")
# RNAWhite<-read.table("/lab/solexa_weissman/cweng/Genomes/10X/WhiteList_10X_Multiome.RNA")
Dic2<-ATACWhite$V1
names(Dic2)<-as.character(RNAWhite$V1)
meta$ATACName<-Dic2[gsub(RNAclusterPost,"",row.names(meta))]
return(meta)
}


#' Function to Merge sparse Matrix
#'
#' This function allows you to input a list of sparse matrix and merge by rownames, return a new sparse matrix
#' @param mtx.list  A list of sparse matrix to be merged
#' @param postfix  a vector of postfix (Usually are numbers that added at the end of cell names). Better be consistent with a merged MitoTracing object orders
#' @return new sparse matrix
#' @examples
#' Donor4_HSC_HPC_BMMC.Mtx<-MergeMtx(list(Donor04_BMMC_Multiome_wrapper$seurat@assays$RNA@counts,Donor04_HPC_Multiome_wrapper$seurat@assays$RNA@counts,Donor04_HSC_Multiome_wrapper$seurat@assays$RNA@counts),c(3,2,1))
#' Donor4_HSC_HPC_BMMC.RNA.seurat<-GEM_Wrapper(Donor4_HSC_HPC_BMMC.Mtx)
#' @export
MergeMtx<-function(mtx.list,postfix){
colnames(mtx.list[[1]])<-strsplit(colnames(mtx.list[[1]]),"-") %>% lapply(.,function(x){x[1]}) %>% unlist %>% paste(.,postfix[1],sep="-")
Merged.mtx<-as.matrix(mtx.list[[1]])
for(i in 2:length(mtx.list)){
    colnames(mtx.list[[i]])<-strsplit(colnames(mtx.list[[i]]),"-") %>% lapply(.,function(x){x[1]}) %>% unlist %>% paste(.,postfix[i],sep="-")
    Merged.mtx<-Tomerge_v2(Merged.mtx,as.matrix(mtx.list[[i]]),leavex = T, leavey = T)
}
Merged.mtx[is.na(Merged.mtx)]<-0
Merged.mtx<-Matrix(as.matrix(Merged.mtx))
return(Merged.mtx)
}
