#' @title  Load example datasets
#' @description The default example files
#' @param exampleData Indicate which type of default datasets to import
#' @return  return the default datasets required by the users
#' @usage #obtain the data for gene expression.
#' @aliases
#' library(eMutPath)
#' Input.exp=GetExampleData(exampleData="Exp.input")
#' #obtain the sample label
#' Cancer_s=GetExampleData(exampleData="Cancer_s")
#' Normal_s=GetExampleData(exampleData="Normal_s")
#' #obtain the protein interaction network
#' Network=GetExampleData(exampleData="network")
#' #obtain the mutation data
#' Mut=GetExampleData(exampleData="mut")
#'@export
GetExampleData <- function(exampleData){
  if(!exists("envData")) envData<-initialize_data()
  ##############get the expression file
  if (exampleData=="Exp.input")
  {
    Exp.input<-get("Exp.input",envir=envData)

    return(Exp.input)
  }
  #############get the network file
  if (exampleData=="network")
  {
    network<-get("network",envir=envData)

    return(network)
  }
  ##############get the mutation file
  if (exampleData=="mut")
  {
    mut<-get("mut",envir=envData)

    return(mut)
  }
  #############get the columns of cancer samples
  if (exampleData=="Cancer_s")
  {
    Cancer_s<-get("Cancer_s",envir=envData)

    return(Cancer_s)
  }
  #############get the columns of normal samples
  if (exampleData=="Normal_s")
  {
    Normal_s<-get("Normal_s",envir=envData)

    return(Normal_s)
  }
}
