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
#' @param bclength The cell barcode length, default is 16
#' @param from A vector of the postfix,  usually is c(1,2,3,...), it depends on how many samples are aggregated in Cellranger RNA part
#' @param to A vector of the postfix, those cooresponds to the postfix added in scMitoTracing, in general, if it matches, then simply c(1,2,3,...), 
#' but in case not match, here provides a way to transform into scMitoTracing order
#' @return meta a dataframe
#' @examples
#' Translate_RNA2ATAC(meta)
#' @export
Translate_RNA2ATAC<-function(meta=bmmc.filtered@meta.data,PostFix=T,bclength=16,from=c(1,2,3),to=c(1,2,3)){
data(ATACWhite)
data(RNAWhite)
L<-nchar(row.names(meta))
post<-substr(row.names(meta),bclength+2,L)
post<-post %>% plyr::mapvalues(.,from=from,to=to)
NakeName<-substr(row.names(meta),1,bclength)
Dic2<-ATACWhite$V1
names(Dic2)<-as.character(RNAWhite$V1)
if(PostFix){
	meta$ATACName<-paste(Dic2[NakeName],post,sep="_")
}else{
	meta$ATACName<-Dic2[NakeName]
}
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
