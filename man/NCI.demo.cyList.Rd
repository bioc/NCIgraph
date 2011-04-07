\name{NCI.demo.cyList}
\docType{data}
\alias{NCI.demo.cyList}
\title{10 raw NCI networks from Nature curated pathways and BioCarta
  imported as graphNEL objects, for demonstration purpose}
\description{These are the ten first elements of the full list of raw
  networks that can be downloaded using the downloadCyLists.R script.}
\usage{NCI.demo.cyList}
\format{A list of 10 graphNEL objects.}
\examples{
data("NCIgraphVignette")
length(NCI.demo.cyList)

library(Rgraphviz)
plot(NCI.demo.cyList[[1]])
}

\author{Laurent Jacob}
\keyword{datasets}
