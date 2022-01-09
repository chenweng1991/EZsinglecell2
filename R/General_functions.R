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
