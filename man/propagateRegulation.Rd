%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  ./propagateRegulation.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{propagateRegulation}
\alias{propagateRegulation}


\title{Transforms the network in a way that each Biochemical
Reaction node pointing to a Complex points to what is regulated by
the complex and updates the interaction types accordingly}


\description{
 Transforms the network in a way that each Biochemical
Reaction node pointing to a Complex points to what is regulated by
the complex and updates the interaction types accordingly.
}

\usage{propagateRegulation(g)}

\arguments{
  \item{g}{A \code{\link[=graph-class]{graph}} object.}
}

\value{The updated \code{\link[=graph-class]{graph}} object}

\author{Laurent Jacob}

\seealso{
  \code{\link{parseNCInetwork}}()
}

