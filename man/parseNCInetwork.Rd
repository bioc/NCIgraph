%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  ./parseNCInetwork.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{parseNCInetwork}
\alias{parseNCInetwork}


\title{Takes a NCI network and transforms it into a simpler graph only
representing inhibition/activation relationships between genes}


\description{
 Takes a NCI network and transforms it into a simpler graph only
representing inhibition/activation relationships between genes.
}

\usage{parseNCInetwork(g, propagateReg=TRUE, separateEntrez=TRUE, mergeEntrezCopies=TRUE, entrezOnly=TRUE)}

\arguments{
  \item{g}{A \code{\link[=graph-class]{graph}} object.}
  \item{propagateReg}{A \code{\link[base]{logical}}. If \code{\link[base:logical]{TRUE}}, use propagateRegulation to transform the network before parsing it.}
  \item{separateEntrez}{A \code{\link[base]{logical}}. If \code{\link[base:logical]{TRUE}}, don't merge two nodes with entrezID.}
  \item{mergeEntrezCopies}{A \code{\link[base]{logical}}. If \code{\link[base:logical]{TRUE}}, merge resulting nodes that have the same entrezID.}
  \item{entrezOnly}{A \code{\link[base]{logical}}. If \code{\link[base:logical]{TRUE}}, only keep nodes with an entrezID property.}
}

\value{The new \code{\link[=graph-class]{graph}} object.}

\author{Laurent Jacob}

\examples{
## Load NCIgraph

library(NCIgraph)

## Get some raw networks

data("NCIgraphVignette", package="NCIgraph")

## Parse the first of them

parsedNetwork <- parseNCInetwork(NCI.demo.cyList[[1]],propagateReg=TRUE,separateEntrez=TRUE,mergeEntrezCopies=TRUE,entrezOnly=TRUE)
}

