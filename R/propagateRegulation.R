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
# @RdocFunction propagateRegulation
##
## @title "Transforms the network in a way that each Biochemical
## Reaction node pointing to a Complex points to what is regulated by
## the complex and updates the interaction types accordingly"
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
## \value{The updated \code{\link[=graph-class]{graph}} object}
##
## @author
##
## \seealso{
##   @see "parseNCInetwork"
## }
##
##
##*/########################################################################


propagateRegulation <- function(g)
  {
    nd <- g@nodeData@data

    goRegURL <- "<LI>- <A class=\"link\" HREF=\"http://www.godatabase.org/cgi-bin/amigo/go.cgi?open_1=0006350\">GO: 0006350</A></LI>"
    
    toAdd <- c()
    toRemove <- c()
    edgeSign <- c()
    
    ## List Biochemical Reaction nodes which are not transcription regulations
    findNonRegulations <- function(e)
      {
        isAReaction <- isTRUE(e$biopax.entity_type == 'Biochemical Reaction')
        isNotARegulation <- !(isTRUE(e$biopax.xref.GO == '0006350')
                             || isTRUE(e$biopax.relationship_references == goRegURL))
        return(isAReaction && isNotARegulation)
      }
    nd <- nd[unlist(lapply(nd,FUN=findNonRegulations))]
    for(n in seq(along=nd))
      {
        cReac <- names(nd)[n]
        
        ## Get genes which are regulated by the target of the current reaction
        ntgtStr <- directedBFS(g,cReac)
        newTgts <- ntgtStr$tgtNodes
        regSign <- ntgtStr$regSign

        ## Anything pointing to the reaction should now point to these targets, not the reaction anymore
        parents <- names(edges(g))[unlist(lapply(edges(g),FUN=function(e) cReac %in% e))]
        for(pp in parents)
          {
            oldType <- edgeData(g,pp,cReac,'edgeType')
            if(!oldType %in% c('INHIBITION','ACTIVATION'))
              next
            for(nn in newTgts)
              {
                toAdd <- cbind(toAdd,c(pp,nn))
                if(oldType == c('INHIBITION','ACTIVATION')[1 + (regSign[[nn]] + 1)/2])
                  edgeSign <- c(edgeSign,'ACTIVATION')
                else
                  edgeSign <- c(edgeSign,'INHIBITION')
              }
            toRemove <- cbind(toRemove,c(pp,cReac))
          }
      }

    ## Add and remove edges
    for(e in seq(along=edgeSign))
      {
        g <- addEdge(toAdd[1,e],toAdd[2,e],g)
        edgeData(g,toAdd[1,e],toAdd[2,e],'edgeType') <- edgeSign[e]
      }
    for(e in seq(along=toRemove[1,]))
      g <- removeEdge(toRemove[1,e],toRemove[2,e],g)
    
    return(g)
  }
