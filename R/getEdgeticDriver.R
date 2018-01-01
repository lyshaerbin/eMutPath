#'@title Identify the edgetic driver mutations
#'@description Identify the driver mutations that perturbe the protein interaction network
#'@param Dysnet Dysregulated network for samples
#'@param Mut The mutation files for samples
#'@param alpha The significance level
#'@param n.sim The simulation times defined by users
#'@return MutPPI the mutation perturbed PPIs in cancer samples
#'@usage #Identifying the mutation mediated-dysregulated edges in both two conditions
#'@aliases
#' #obtain the data for gene expression.
#' Input.exp=GetExampleData(exampleData="Exp.input")
#' #obtain the sample label
#' Cancer_s=GetExampleData(exampleData="Cancer_s")
#' Normal_s=GetExampleData(exampleData="Normal_s")
#' #obtain the protein interaction network
#' Network=GetExampleData(exampleData="network")
#' #obtain the mutation data
#' Mut=GetExampleData(exampleData="mut")
#' DysCN=EdgeticDys_CN(Input.exp,Network,thr=0.01)
#' Input.exp2=Input.exp[,Cancer_s]
#' DysC=EdgeticDys_CN(Input.exp2,Network,thr=0.01)
#' Dys.net=EdgeticDys_both(DysCN,DysC)
#' Driver.mut=getEdgeticDriver(Dys.net,Mut,alpha=0.05,n.sim=1000)
#'@export
getEdgeticDriver <- function(Dysnet,Mut,alpha,n.sim){
  U.edge=unique(Dysnet[,2:3])
  #######################keep the interactions with at least one gene mutated
  Mutdys=c()
  for(i in 1:dim(U.edge)[1]){
    if(i%%100==0) {print(i)}
    xa=which(Mut$Hugo_Symbol==U.edge$GeneA[i])
    xb=which(Mut$Hugo_Symbol==U.edge$GeneB[i])
    if(length(xa)>0|length(xb)>0){
      Mutdys=rbind(Mutdys,U.edge[i,])
    }
  }
  #########################Mutation, PPI, #sample-1, #sample-2, #overlap, p
  Three.sam=function(fusion_sample){
    fusionthree=c()
    for (i in 1:length(fusion_sample)) {
      aa=strsplit(as.character(fusion_sample[i]),"-")
      bb=paste(aa[[1]][1],aa[[1]][2],aa[[1]][3],sep="-")
      fusionthree=rbind(fusionthree,bb)
    }
    return(fusionthree)
  }
  N11=length(unique(Dysnet$Sam.ID))#######the number of samples with dysregulated edges
  N22=length(unique(Mut$Tumor_Sample_Barcode))###########the number of samples with mutation
  FF.net=c()
  for(ii in 1:dim(Mutdys)[1]){
    if(ii%%100==0) {print(ii)}
    xx=which(Dysnet$GeneA==Mutdys$GeneA[ii]&Dysnet$GeneB==Mutdys$GeneB[ii])
    Dys.sam=unique(Dysnet$Sam.ID[xx])
    A=Three.sam(Dys.sam)
    Xmut1=which(Mut$Hugo_Symbol==Mutdys$GeneA[ii])
    Xmut2=which(Mut$Hugo_Symbol==Mutdys$GeneB[ii])
    Xmut=union(Xmut1,Xmut2)
    Gene.mut=Mut[Xmut,]
    Gene.mut=unique(Gene.mut[,1:7])
    for(j in 1:dim(Gene.mut)[1]){
      y1=which(Mut$Hugo_Symbol==Gene.mut$Hugo_Symbol[j]&
                 Mut$Chromosome==Gene.mut$Chromosome[j]&
                 Mut$Start_Position==Gene.mut$Start_Position[j]&
                 Mut$End_Position==Gene.mut$End_Position[j]&
                 Mut$Variant_Classification==Gene.mut$Variant_Classification[j]&
                 Mut$Tumor_Seq_Allele1==Gene.mut$Tumor_Seq_Allele1[j]&
                 Mut$Tumor_Seq_Allele2==Gene.mut$Tumor_Seq_Allele2[j])
      Mut.sam=unique(Mut$Tumor_Sample_Barcode[y1])
      B=Three.sam(Mut.sam)
      #Simulation p-value
      sim=unlist(lapply(1:n.sim,
                        function(i){AA=sample(1:N11,length(A));BB=sample(1:N22,length(B));return(sum(AA %in% BB))}))
      P=length(which(sim>=length(intersect(A,B))))/n.sim

      WKK=cbind(Gene.mut[j,],Mutdys[ii,],length(B),length(A),length(intersect(A,B)),P)
      FF.net=rbind(FF.net,WKK,stringsAsFactors = FALSE)
    }
  }
  Sig=which(FF.net$P<alpha)
  MutPPI=FF.net[Sig,]
  return(MutPPI)
}
