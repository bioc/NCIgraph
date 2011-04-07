## Copyright 2011 Laurent Jacob

## This file is part of NCIgraph.

## NCIgraph is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.

## NCIgraph is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with NCIgraph.  If not, see <http://www.gnu.org/licenses/>.

#########################################################################/**
# @RdocFunction parseNCInetwork
##
## @title "Takes a NCI network and transforms it into a simpler graph only
## representing inhibition/activation relationships between genes"
##
##
## \description{
##  @get "title".
## }
##
## @synopsis
##
## \arguments{
##   \item{g}{A \code{\link[=graph-class]{graph}} object.}
##   \item{propagateReg}{A @logical. If @TRUE, use propagateRegulation to transform the network before parsing it.}
##   \item{separateEntrez}{A @logical. If @TRUE, don't merge two nodes with entrezID.}
##   \item{mergeEntrezCopies}{A @logical. If @TRUE, merge resulting nodes that have the same entrezID.}
##   \item{entrezOnly}{A @logical. If @TRUE, only keep nodes with an entrezID property.}
## }
##
## \value{The new \code{\link[=graph-class]{graph}} object.}
##
## @author
##
## @examples "../inst/extdata/parseNCInetwork.Rex"
##
##*/########################################################################

##source('edgesToMerge.R')
##source('mergeNodes.R')
##source('propagateRegulation.R')

parseNCInetwork <- function(g,propagateReg=TRUE,separateEntrez=TRUE,mergeEntrezCopies=TRUE,entrezOnly=TRUE)
  {

    ## If there is no gene in the network, don't bother parsing it and return an empty graph
    entrezNodes <- nodes(g)[unlist(lapply(g@nodeData@data,FUN=function(e) 'biopax.xref.ENTREZGENE' %in% names(e)))]
    if(entrezOnly && length(entrezNodes) == 0)
      return(subGraph(entrezNodes,g))

    if(propagateReg)
      g <- propagateRegulation(g)
    
    ng <- nodes(g)
    p <- length(ng)

    ## Merge edges
    mEdges <- edgesToMerge(g)
    if(length(mEdges) > 0)
      pg <- mergeNodes(g,mEdges,separateEntrez=separateEntrez,entrezOnly=entrezOnly)
    else
      pg <- g

    if(mergeEntrezCopies)
      {
        ## Collapse nodes corresponding to the same entrez id
        mEdges <- list()
        nodeNames <- names(pg@nodeData@data)
        for(n in seq(along=pg@nodeData@data))
          {
            nodeE <- pg@nodeData@data[[n]]$biopax.xref.ENTREZGENE
            for(nn in seq(along=pg@nodeData@data)[-n])
              {
                nnodeE <- pg@nodeData@data[[nn]]$biopax.xref.ENTREZGENE
                names(nodeE) <- names(nnodeE) <-  NULL
                if(isTRUE(nodeE == nnodeE))
                  mEdges <- c(mEdges,list(c(nodeNames[n],nodeNames[nn])))
              }          
          }
        if(length(mEdges) > 0)
          pg <- mergeNodes(pg,mEdges,separateEntrez=FALSE,entrezOnly=entrezOnly) # By definition, we don't want to separate entrezID nodes in this step
      }

    return(pg)
  }


############################################################################
# HISTORY:
# 2011-02-14
# o Created.
############################################################################
