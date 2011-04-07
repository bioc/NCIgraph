%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  ./directedBFS.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{directedBFS}
\alias{directedBFS}


\title{Uses a breadth first search on a directed graph to identify
which genes are regulated by a particular node in the graph}


\description{
 Uses a breadth first search on a directed graph to identify
which genes are regulated by a particular node in the graph.
}

\usage{directedBFS(g, node)}

\arguments{
  \item{g}{A \code{\link[=graph-class]{graph}} object.}
  \item{node}{A node of g.}
}

\value{A structured \code{\link[base]{list}} containing the regulated genes and the
type of interaction between node and each gene.}

\author{Laurent Jacob}

\seealso{
  \code{\link{propagateRegulation}}()
}

