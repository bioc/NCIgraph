## Load NCIgraph

library(NCIgraph)

## Get some raw networks

data("NCIgraphVignette", package="NCIgraph")

## Parse the first of them

parsedNetwork <- parseNCInetwork(NCI.demo.cyList[[1]],propagateReg=TRUE,separateEntrez=TRUE,mergeEntrezCopies=TRUE,entrezOnly=TRUE)
