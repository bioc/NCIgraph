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
# @RdocFunction getSubtype.NCIgraph
##
## @title "Returns a list of @KEGGEdgeSubType objects describing each
## edge of the NCI network"
##
##
## \description{
##  @get "title".
## }
##
## @synopsis
##
## \arguments{
##   \item{object}{An \code{\link[=NCIgraph-class]{NCIgraph}} object.}
## }
##
## \value{A @list of KEGGEdgeSubType objects.}
##
## @author
##
## @examples "../inst/extdata/getSubtype.NCIgraph.Rex"
##
##*/########################################################################

getSubtype.NCIgraph <- function(object)
{

  ## - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - -
  ## Argument 'object'
  if (!validGraph(object)) {
    throw("Argument 'object' is not a valid graph object")
  }

  stFromNodeData <- function(e)
    {
      stName <- unique(tolower(e$edgeType))
      if(length(stName) > 1)
        stName <- 'mixt'
      else
        stName <- as.character(stName)
      arr <- '???'
      if(stName == 'activation')
        arr <- '<--'
      if(stName == 'inhibition')
        arr <- '|--'
      
      st <- new('KEGGEdgeSubType',name=stName,value= arr)
      list(subtype=st) # For consistency with KEGGgraph
    }
  
  lapply(object@edgeData@data, FUN=stFromNodeData)  
}

setMethod('getSubtype', signature('NCIgraph'), getSubtype.NCIgraph)


############################################################################
## HISTORY:
## 2011-03-03
## o Created
############################################################################
