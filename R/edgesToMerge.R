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
# @RdocFunction edgesToMerge
##
## @title " Identifies edges that should be merged to parse a NCI network"
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
## \value{A @list of edges to be merged}
##
## @author
##
## \seealso{
##   @see "parseNCInetwork"
## }
##
##*/########################################################################

edgesToMerge <- function(g)
  {
    ng <- nodes(g)
    p <- length(ng)
    eData <- g@edgeData
    rTypes <- unlist(lapply(eData@data, FUN=function(e) e$edgeType))
    ##edgeList <- matrix(unlist(lapply(names(eData),FUN=function(s) strsplit(s,"\\|"))),2,length(eData@data))
    edgeList <- lapply(names(eData),FUN=function(s) unlist(strsplit(s,"\\|")))

    ## Step 1: Merge controller and contains associations
    mEdges <- edgeList[rTypes %in% c("CONTROLLER","CONTAINS")]

    ## For LEFT-RIGHT, it can be more subtle: don't include those
    ## corresponding to complex assembly/disassembly

    lrEdges <- edgeList[rTypes %in% c("LEFT","RIGHT")]
    
    isAComplexAssembly <- function(e)
      "Complex Assembly" %in%
    c(g@nodeData@data[[e[1]]]$biopax.entity_type,g@nodeData@data[[e[2]]]$biopax.entity_type)

    if(length(lrEdges)>0)
      mEdges <- c(mEdges,lrEdges[!unlist(lapply(lrEdges, FUN=isAComplexAssembly))])
  }
    
############################################################################
# HISTORY:
# 2011-21-14
# o Created.
############################################################################
