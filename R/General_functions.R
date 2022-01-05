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
