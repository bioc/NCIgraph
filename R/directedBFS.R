## Copyright 2010 Laurent Jacob

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
# @RdocFunction directedBFS
##
## @title "Uses a breadth first search on a directed graph to identify
## which genes are regulated by a particular node in the graph"
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
##   \item{node}{A node of g.}
## }
##
## \value{A structured @list containing the regulated genes and the
## type of interaction between node and each gene.}
##
## @author
##
## \seealso{
##   @see "propagateRegulation"
## }
##
##
##*/########################################################################

directedBFS <- function(g,node)
  {
    el <- edgeL(g)
    nnames <- nodes(g)
    ndata <- g@nodeData@data
    queued <- c(node)
    visited <- c(node)
    tgtNodes <- c()
    regSign <- list()
    regSign[[node]] <- 1

    goRegURL <- "<LI>- <A class=\"link\" HREF=\"http://www.godatabase.org/cgi-bin/amigo/go.cgi?open_1=0006350\">GO: 0006350</A></LI>"
    
    while(length(queued != 0))
      {
        v <- queued[1] # Unqueue one node v
        ##if(node == "-305921021789615134-pid_i_103676")
        ##  print(v)
        queued <- queued[-1]
        ##if(v == "-305921021789615161-pid_i_103682")
        ##  print(nnames[el[[v]]$edges] %in% visited)
        for(vv in nnames[el[[v]]$edges]) # Loop on neighbors of v
          if(! vv %in% visited)
            {
              visited <- c(visited,vv)

              ## Sign of regulation from root
              if(isTRUE(g@edgeData@data[[paste(v,'|',vv,sep='')]]$edgeType == 'INHIBITION'))
                regSign[[vv]] <- -regSign[[v]]
              else
                regSign[[vv]] <- regSign[[v]]
              
              isRegulated <- (isTRUE(ndata[[v]]$biopax.xref.GO == '0006350')
              || isTRUE(ndata[[v]]$biopax.relationship_references == goRegURL))

              isGene <- isTRUE(ndata[[vv]]$biopax.entity_type == 'Protein')
              
              if(isRegulated && isGene)
                tgtNodes <- c(tgtNodes,vv)
              else
                queued <- c(queued,vv)
            }
      }
    
    return(list(tgtNodes=tgtNodes,regSign=regSign[tgtNodes]))
  }
