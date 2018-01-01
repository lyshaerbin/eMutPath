#' @title Get dysregulated edges in cancer
#' @description Identify the dysregulated edges for each sample
#' @param DysCN The dysregulated network in cancer vs normal
#' @param DysC The dysregulated network for specific cancer sample
#' @usage #Identifying the dysregulated edges in both two conditions
#' @aliases
#' Input.exp=GetExampleData(exampleData="Exp.input")
#' #obtain the sample label
#' Cancer_s=GetExampleData(exampleData="Cancer_s")
#' Normal_s=GetExampleData(exampleData="Normal_s")
#' #obtain the protein interaction network
#' Network=GetExampleData(exampleData="network")
#' #identify the dysregulated edges
#' DysCN=EdgeticDys_CN(Input.exp,Network,thr=0.01)
#' Input.exp=Input.exp[,Cancer_s]
#' DysC=EdgeticDys_CN(Input.exp,Network,thr=0.01)
#' Dys.net=EdgeticDys_both(DysCN,DysC)
#'@export
EdgeticDys_both <- function(DysCN,DysC){
  n1=dim(DysCN)[2]
  n2=dim(DysC)[2]
  Cname1=colnames(DysCN)[1:3]
  Cname2=colnames(DysC)[1:3]
  if(n1==n2){
    Both.dys=merge(DysCN,DysC,by.x=Cname1,by.y=Cname2)
    return(Both.dys)
  } else {
    print("error: The columns of two inputs should be consistent")
  }
}
