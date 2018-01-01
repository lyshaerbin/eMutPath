# eMutPath
Package: eMutPath

Type: R Package

Title: eMutPath (edge-based mutation prioritization based on network perturbation paths)

Version: 0.1.0

Author: Yongsheng Li

Maintainer: Yongsheng Li <YLi42@mdanderson.org>

Description: Identification of edgetic driver mutations by integrating gene
        expression and protein interaction network or functional pathways.
        
License: GPL (>= 2)

Encoding: UTF-8

LazyData: true

Imports: outliers,igraph,stats,utils

RoxygenNote: 6.0.1

NeedsCompilation: no

Packaged: 2017-12-18 14:19:49 UTC; YLi42

Usages:

Input.exp=GetExampleData(exampleData="Exp.input")

#' #obtain the sample label

Cancer_s=GetExampleData(exampleData="Cancer_s")

Normal_s=GetExampleData(exampleData="Normal_s")

#' #obtain the protein interaction network

Network=GetExampleData(exampleData="network")

#' #obtain the mutation data

Mut=GetExampleData(exampleData="mut")

DysCN=EdgeticDys_CN(Input.exp,Network,thr=0.01)

Input.exp2=Input.exp[,Cancer_s]

DysC=EdgeticDys_CN(Input.exp2,Network,thr=0.01)

Dys.net=EdgeticDys_both(DysCN,DysC)

Driver.mut=getEdgeticDriver(Dys.net,Mut,alpha=0.05,n.sim=1000)
