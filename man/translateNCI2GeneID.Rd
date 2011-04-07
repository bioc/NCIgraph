%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  ./translateNCI2GeneID.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{translateNCI2GeneID}
\alias{translateNCI2GeneID}


\title{Gives the entrezID corresponding to the nodes of a graph}


\description{
 Gives the entrezID corresponding to the nodes of a graph.
}

\usage{translateNCI2GeneID(g)}

\arguments{
  \item{g}{A \code{\link[=graph-class]{graph}} object.}
}

\value{A vector of \code{\link[base]{character}} giving the entrez ID of the nodes of g.}

\author{Laurent Jacob}

\seealso{
  \code{\link{parseNCInetwork}}()
}

\examples{
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
}

