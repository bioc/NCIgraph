%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  ./getSubtype.NCIgraph.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getSubtype.NCIgraph}
\alias{getSubtype.NCIgraph}


\title{Returns a list of @KEGGEdgeSubType objects describing each
edge of the NCI network}


\description{
 Returns a list of @KEGGEdgeSubType objects describing each
edge of the NCI network.
}

\usage{getSubtype.NCIgraph(object)}

\arguments{
  \item{object}{An \code{\link[=NCIgraph-class]{NCIgraph}} object.}
}

\value{A \code{\link[base]{list}} of KEGGEdgeSubType objects.}

\author{Laurent Jacob}

\examples{
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

}

