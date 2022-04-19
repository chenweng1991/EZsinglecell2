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

#' Function to add hematopoietic signatures from Griffin_Signatures
#'
#' This function allows you to input a seurat object, add the signatures and return an seurat object
#' @param object a seurat object
#' @return a seurat object
#' @export
AddHemSignature<-function(object=Donor01_BMMC_Multiome_wrapper.filtered){
require(Seurat)
data(Griffin_Signatures)
Sig.HSC<-list(Griffin_Signatures$Griffin_BPDCN_HSC %>% .[!is.na(.)])
Sig.Prog<-list(Griffin_Signatures$Griffin_BPDCN_Prog %>% .[!is.na(.)])
Sig.EarlyE<-list(Griffin_Signatures$Griffin_BPDCN_EarlyE %>% .[!is.na(.)])
Sig.LateE<-list(Griffin_Signatures$Griffin_BPDCN_LateE	 %>% .[!is.na(.)])
Sig.ProMono<-list(Griffin_Signatures$Griffin_BPDCN_ProMono %>% .[!is.na(.)])
Sig.Mono<-list(Griffin_Signatures$Griffin_BPDCN_Mono %>% .[!is.na(.)])
Sig.ncMono<-list(Griffin_Signatures$Griffin_BPDCN_ncMono %>% .[!is.na(.)])
Sig.cDC<-list(Griffin_Signatures$Griffin_BPDCN_cDC %>% .[!is.na(.)])
Sig.pDC<-list(Griffin_Signatures$Griffin_BPDCN_pDC %>% .[!is.na(.)])
Sig.ProB<-list(Griffin_Signatures$Griffin_BPDCN_ProB %>% .[!is.na(.)])
Sig.PreB<-list(Griffin_Signatures$Griffin_BPDCN_PreB %>% .[!is.na(.)])
Sig.B<-list(Griffin_Signatures$Griffin_BPDCN_B %>% .[!is.na(.)])
Sig.Plasma<-list(Griffin_Signatures$Griffin_BPDCN_Plasma %>% .[!is.na(.)])
Sig.T<-list(Griffin_Signatures$Griffin_BPDCN_T %>% .[!is.na(.)])
Sig.CTL<-list(Griffin_Signatures$Griffin_BPDCN_CTL %>% .[!is.na(.)])
Sig.NK<-list(Griffin_Signatures$Griffin_BPDCN_NK %>% .[!is.na(.)])
DefaultAssay(object)<-"SCT"
object<-AddModuleScore(object = object, features = Sig.HSC, name = "Sig.HSC")
object<-AddModuleScore(object = object, features = Sig.Prog, name = "Sig.Prog")
object<-AddModuleScore(object = object, features = Sig.EarlyE, name = "Sig.EarlyE")
object<-AddModuleScore(object = object, features = Sig.LateE, name = "Sig.LateE")
object<-AddModuleScore(object = object, features = Sig.ProMono, name = "Sig.ProMono")
object<-AddModuleScore(object = object, features = Sig.Mono, name = "Sig.Mono")
object<-AddModuleScore(object = object, features = Sig.ncMono, name = "Sig.ncMono")
object<-AddModuleScore(object = object, features = Sig.cDC, name = "Sig.cDC")
object<-AddModuleScore(object = object, features = Sig.pDC, name = "Sig.pDC")
object<-AddModuleScore(object = object, features = Sig.ProB, name = "Sig.ProB")
object<-AddModuleScore(object = object, features = Sig.PreB, name = "Sig.PreB")
object<-AddModuleScore(object = object, features = Sig.B, name = "Sig.B")
object<-AddModuleScore(object = object, features = Sig.Plasma, name = "Sig.Plasma")
object<-AddModuleScore(object = object, features = Sig.T, name = "Sig.T")
object<-AddModuleScore(object = object, features = Sig.CTL, name = "Sig.CTL")
object<-AddModuleScore(object = object, features = Sig.NK, name = "Sig.NK")
return(object)    
}


#' Function to reclustering a seurat object
#'
#' This function allows you to input a seurat object(multiome), redo clustering. Usually this is after subset
#' @param ob a seurat object
#' @return a seurat object
#' @export
Reclustering<-function(ob){
DefaultAssay(ob) <- "RNA"
ob <- SCTransform(ob, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
DefaultAssay(ob) <- "ATAC"
ob <- RunTFIDF(ob)
ob <- FindTopFeatures(ob, min.cutoff = 'q0')
ob <- RunSVD(ob)
ob <- RunUMAP(ob, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
ob <- FindMultiModalNeighbors(ob, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
ob <- RunUMAP(ob, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
ob <- FindClusters(ob, graph.name = "wsnn", algorithm = 3, verbose = FALSE)
return(ob)
}