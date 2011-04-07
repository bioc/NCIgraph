%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  ./mergeNodes.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{mergeNodes}
\alias{mergeNodes}


\title{Merges a given list of nodes in a graph}


\description{
 Merges a given list of nodes in a graph.
}

\usage{mergeNodes(g, mEdges, separateEntrez=TRUE, entrezOnly=TRUE)}

\arguments{
  \item{g}{A \code{\link[=graph-class]{graph}} object.}
  \item{mEdges}{A \code{\link[base]{list}} of nodes to be merged.}
  \item{separateEntrez}{A \code{\link[base]{logical}}. If \code{\link[base:logical]{TRUE}}, don't merge two nodes with entrezID.}
  \item{entrezOnly}{A \code{\link[base]{logical}}. If \code{\link[base:logical]{TRUE}}, only keep nodes with an entrezID property.}
}

\value{The updated \code{\link[=graph-class]{graph}} object}

\author{Laurent Jacob}

\seealso{
  \code{\link{parseNCInetwork}}()
}

