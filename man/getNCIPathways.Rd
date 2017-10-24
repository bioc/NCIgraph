%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  ./getNCIPathways.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getNCIPathways}
\alias{getNCIPathways}


\title{Loads networks from Cytoscape and parses them}


\description{
 Loads networks from Cytoscape and parses them.
}

\usage{getNCIPathways(cyList=NULL, parseNetworks=TRUE, entrezOnly=TRUE, verbose=FALSE)}

\arguments{
  \item{cyList}{a \code{\link[base]{list}} providing the networks loaded from Cytoscape. If \code{\link[base]{NULL}}, the function will try to build the \code{\link[base]{list}} from Cytoscape.}
  \item{verbose}{If \code{\link[base:logical]{TRUE}}, extra information is output.}
  \item{parseNetworks}{A \code{\link[base]{logical}}. If \code{\link[base:logical]{FALSE}}, the raw NCI networks are
  returned as graphNEL objects. If \code{\link[base:logical]{TRUE}}, some additional parsing
  is performed by the parseNCInetwork function.}
  \item{entrezOnly}{A \code{\link[base]{logical}}. If \code{\link[base:logical]{TRUE}}, only keep nodes with an entrezID property.}
}

\value{A \code{\link[base]{list}} of two elements: pList, a \code{\link[base]{list}} of graphNEL objects,
and failedW a \code{\link[base]{list}} containing the names of the networks that R
failed to read from cytoscape.}

\author{Laurent Jacob}

\seealso{
  \code{\link{parseNCInetwork}}()
}

\examples{
##------------------------------
## Load NCIgraph
##------------------------------

library(NCIgraph)

##------------------------------
## Example 1: with Cytoscape
##------------------------------

## Must have Cytoscape running with some networks open and CyREST plugin started.

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
}
