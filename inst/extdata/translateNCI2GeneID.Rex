##------------------------------
## Load NCIgraph
##------------------------------

library(NCIgraph)

## Get some raw networks

data("NCIgraphVignette", package="NCIgraph")
 
## Parse them

grList <- getNCIPathways(cyList=NCI.demo.cyList, parseNetworks=TRUE, entrezOnly=TRUE, verbose=TRUE)$pList

## Get the gene ids for the first of them

gids <- translateNCI2GeneID(grList[[1]])