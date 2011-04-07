##------------------------------
## Load NCIgraph
##------------------------------

library(NCIgraph)

##------------------------------
## Get some raw networks
##------------------------------

data("NCIgraphVignette", package="NCIgraph")

##------------------------------
## Parse them
##------------------------------

grList <- getNCIPathways(cyList=NCI.demo.cyList, parseNetworks=TRUE, entrezOnly=TRUE, verbose=TRUE)$pList

##----------------------------------------------------------------
##
## Get the subtype of the second network. Some activation and some
## inhibition edges.
##
##----------------------------------------------------------------

getSubtype.NCIgraph(grList[[2]])

