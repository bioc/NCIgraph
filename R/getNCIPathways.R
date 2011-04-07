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
# @RdocFunction getNCIPathways
##
## @title "Loads networks from Cytoscape and parses them"
##
##
## \description{
##  @get "title".
## }
##
## @synopsis
##
## \arguments{
##   \item{cyList}{a @list providing the networks loaded from Cytoscape. If @NULL, the function will try to build the @list from Cytoscape.}
##   \item{verbose}{If @TRUE, extra information is output.}
##   \item{parseNetworks}{A @logical. If @FALSE, the raw NCI networks are
##   returned as graphNEL objects. If @TRUE, some additional parsing
##   is performed by the parseNCInetwork function.}
##   \item{entrezOnly}{A @logical. If @TRUE, only keep nodes with an entrezID property.}
## }
##
## \value{A @list of two elements: pList, a @list of graphNEL objects,
## and failedW a @list containing the names of the networks that R
## failed to read from cytoscape.}
##
## @author
##
## \seealso{
##   @see "parseNCInetwork"
## }
##
## @examples "../inst/extdata/getNCIPathways.Rex"
##
##*/########################################################################

getNCIPathways <- function(cyList=NULL, parseNetworks=TRUE, entrezOnly=TRUE, verbose=FALSE)
{

  failedW <- c()
  if(is.null(cyList))
    {      
      ## Get a list of all cytoscape windows
      tt <- try({cy <- CytoscapeConnection(); wList <- getWindowList(cy)})
      
      if (class(tt)=="try-error")
        stop("Failed to read graphs from Cytoscape. Check that Cytoscape is open and the Cytoscape RPC plugin is loaded.");
      
      ## Convert each cytoscape graph to a graphNEL
      cyList <- list()
      
      for(w in wList)
        {
          tt <- try(cw.cellCycle <- existing.CytoscapeWindow(w, copy.graph=TRUE))          
          if (class(tt)=="try-error")
            {
              warning(sprintf("Failed to load network %s to R",w))
              failedW <- c(failedW,w)
            }
          else
            cyList[[w]] <- cw.cellCycle@graph
        }
    }

  ## List of edge types present in the database
  ## unique(unlist(lapply(cyList, FUN= function(g) unique(unlist(lapply(g@edgeData@data, FUN=function(e) e$BIOPAX_EDGE_TYPE))))))


  pList <- list()
  
  if(parseNetworks) # Parse each graph to a more usable format
    for(w in names(cyList))
      {
        if(verbose)
          cat("Loading network",w,"\n")
        pList[[w]] <- parseNCInetwork(cyList[[w]],entrezOnly=entrezOnly)
        class(pList[[w]]) <- 'NCIgraph'
      }

  return(list(pList=pList,cyList=cyList,failedW=failedW))
}

############################################################################
# HISTORY:
# 2011-02-02
# o Created.
############################################################################

## To generate the Cytoscape script (to be used through cytoscape command tool plugin):
##for a in networkData/pid.nci.nih.gov/browse_pathways.shtml/*.owl; do echo "network import file=\"$a\"" >> loadNCI.cy; done
