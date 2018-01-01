#' @title Identify the dysregulated network edges
#' @description Identify the dysregulated edges for each sample
#' @param Input.exp The gene expression matrix
#' @param Network The protein interaction network
#' @param thr The significant level for grubbs's test
#' @return All.sam.dys The dyregulated edges for each sample
#' @usage #Identifying the dysregulated edges between cancer and normal samples
#' @aliases EdgeticDys_CN
#' Input.exp=GetExampleData(exampleData="Exp.input")
#' #obtain the sample label
#' Cancer_s=GetExampleData(exampleData="Cancer_s")
#' Normal_s=GetExampleData(exampleData="Normal_s")
#' #obtain the protein interaction network
#' Network=GetExampleData(exampleData="network")
#' #identify the dysregulated edges
#' DysCN=EdgeticDys_CN(Input.exp,Network,thr=0.01)
#'@export
EdgeticDys_CN <- function(Input.exp,Network,thr){
  library(outliers)
  print("Remove the duplicated edges in network\n")
  Gene.input=rownames(Input.exp)
  Net.input=data.frame(N1=Network[,1],N2=Network[,2])
  Net.input=Net.input[!duplicated(Net.input),]
  Net.input=igraph::graph_from_data_frame(Net.input, directed = F, vertices = NULL)
  Net.input=igraph::simplify(Net.input,remove.loops = T)
  Network.input= as.data.frame(igraph::get.edgelist(Net.input))
  #---------------remove edges without expression--------
  print("Now, remove the edges without gene expression")
  Input.net=c()
  for(i in 1:dim(Network.input)[1]){
    if(i%%1000==0) {print(i)}
    x1=which(Gene.input==Network.input$V1[i])
    x2=which(Gene.input==Network.input$V2[i])
    if(length(x1)==1&length(x2)==1){
      Input.net=rbind(Input.net,Network.input[i,])
    }
  }
  ###########################
  grubbs.flag <- function(x,thr) {
    outliers <- NULL
    test <- x
    grubbs.result <- outliers::grubbs.test(test)
    pv <- grubbs.result$p.value
    # throw an error if there are too few values for the Grubb's test
    if (length(test) < 3 ) stop("Grubb's test requires > 2 input values")
    while(pv < thr) {
      outliers <- c(outliers,as.numeric(strsplit(grubbs.result$alternative," ")[[1]][3]))
      test <- x[!x %in% outliers]
      # stop if all but two values are flagged as outliers
      if (length(test) < 3 ) {
        warning("All but two values flagged as outliers")
        break
      }
      grubbs.result <- outliers::grubbs.test(test)
      pv <- grubbs.result$p.value
    }
    return(data.frame(X=x,Outlier=(x %in% outliers)))
  }
  All.dys.net=c() ##### row: net edges, col: cancer+normal sample
  Ma.dis=c()
  M=dim(Input.net)[1]
  print("begin to evaluate the edges:")
  for(i in 1:M){
    #print(i)
    if(i%%1000==0) {print(i)}
    x1=which(Gene.input==Input.net$V1[i])
    x2=which(Gene.input==Input.net$V2[i])
    Edg.exp=rbind(Input.exp[x1,],Input.exp[x2,])
    Edg.exp=t(Edg.exp)
    m_dist <- mahalanobis(Edg.exp[, 1:2], colMeans(Edg.exp[, 1:2]), cov(Edg.exp[, 1:2]))
    DD <- round(m_dist, 2)
    XXA=cbind(colnames(Input.exp),DD)
    xxb=grubbs.flag(as.numeric(as.character(XXA[,2])),thr)
    All.dys.net=rbind(All.dys.net,t(xxb$Outlier))
    Ma.dis=rbind(Ma.dis,t(xxb$X))
  }
  ##############output: sample dysregulated edges(A-B) distance
  All.sam.dys=c()
  All.sam.names=colnames(Input.exp)
  print("begin to prepare the output:")
  for(k1 in 1:length(All.sam.names)){
    if(k1%%100==0) {print(k1)}
    xa1=which(All.dys.net[,k1]==TRUE)
    if(length(xa1)>0){
      SS.dy=cbind(All.sam.names[k1],Input.net[xa1,],Ma.dis[xa1,k1])
      colnames(SS.dy)=c("Sam.ID","GeneA","GeneB","Ma.dis")
      All.sam.dys=rbind(All.sam.dys,SS.dy)
    }
  }
  return(All.sam.dys)
}
