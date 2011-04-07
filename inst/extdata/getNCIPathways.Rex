##------------------------------
## Load NCIgraph
##------------------------------

library(NCIgraph)

##------------------------------
## Example 1: with Cytoscape
##------------------------------

## Must have Cytoscape running with some networks open and Cytoscape RPC plugin started.

## In this case, getNCIPathways will both read the raw networks from Cytoscape and parse them.

\dontrun{
grList <- getNCIPathways(cyList=NULL, parseNetworks=TRUE, entrezOnly=TRUE, verbose=TRUE)$pList
}

##------------------------------
## Example 2: without Cytoscape
##------------------------------

## Get some raw networks

data("NCIgraphVignette", package="NCIgraph")
 
## When passed a non null cyList argument (a list of networks),
## getNCIPathways will simply parse the list of networks

grList <- getNCIPathways(cyList=NCI.demo.cyList, parseNetworks=TRUE, entrezOnly=TRUE, verbose=TRUE)$pList
