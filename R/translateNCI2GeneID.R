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
# @RdocFunction translateNCI2GeneID
##
## @title "Gives the entrezID corresponding to the nodes of a graph"
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
## }
##
## \value{A vector of @character giving the entrez ID of the nodes of g.}
##
## @author
##
## \seealso{
##   @see "parseNCInetwork"
## }
##
## @examples "../inst/extdata/translateNCI2GeneID.Rex"
##
##*/########################################################################

translateNCI2GeneID <- function(g)
  {
    lambdaf <- function(node)
      {
        eid <- node$biopax.xref.ENTREZGENE
        if(length(unique(eid))>1)
          warning(sprintf('More than one entrez ID associated with node %s\n',node$nodeName))
        if(!is.null(eid))
          eid[1]
        else
          NA
      }    
    eids <- unlist(lapply(g@nodeData@data,FUN=lambdaf))
    names(eids) <- nodes(g)
    return(eids)
  }
