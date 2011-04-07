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
# @RdocFunction mergeNodes
##
## @title "Merges a given list of nodes in a graph"
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
##   \item{mEdges}{A @list of nodes to be merged.}
##   \item{separateEntrez}{A @logical. If @TRUE, don't merge two nodes with entrezID.}
##   \item{entrezOnly}{A @logical. If @TRUE, only keep nodes with an entrezID property.}
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

mergeNodes <- function(g,mEdges,separateEntrez=TRUE,entrezOnly=TRUE)
  {
    ng <- nodes(g)
    ed <- g@edgeData@data
    p <- length(ng)

    ##---------------------------------------------------    
    ## Build a graph with the same nodes as g and edges
    ## between any two nodes that should be merged.
    ##---------------------------------------------------
    A <- matrix(0,p,p)
    rownames(A) <- colnames(A) <- ng
    for(c in seq(along=mEdges))
      A[mEdges[[c]][1],mEdges[[c]][2]] <- 1
    A <- A  | t(A)
    lg <- new('graphAM', A)
    cc <- connectedComp(as(lg,'graphNEL')) # The connected components of this graph
                                        # will be the nodes of the new one.

    if(separateEntrez)
      {
        ##--------------------------------------------------------------
        ## Filter: if several entities in the cc have an entrez ID, they
        ## must be assigned a different node each.
        ##--------------------------------------------------------------

        newCC <- list()
        oldCC <- c() # Keep track of original cc to avoid creating spurious edges
        
        ii <- 0
        for(c in cc)
          {
            ii <- ii+1
            eiList <- lapply(g@nodeData@data[c],FUN=function(node) node$biopax.xref.ENTREZGENE)
            eids <- names(unlist(eiList))

            if((length(c) == 1) || (length(eids) <= 1)) # length(eids) = 0 if no node in cc has an entregene
              {
                oldCC <- c(oldCC,ii)
                newCC <- c(newCC,list(c))
                next
              }
            
            ## For each entrez id, remove nodes that are in a different
            ## connected component when removing the other entrez
            ## ids. Duplicate the others.
            for(n in seq(along=eids))
              {
                ## Connected comp without all entrez nodes except n
                sg <- subGraph(c[-match(eids[-n],c)],g)
                ccsg <- connectedComp(sg)
                ccc <-  ccsg[which(unlist(lapply(ccsg,FUN=function(e) eids[n] %in% e)))] # conn comp containing current entrez
                oldCC <- c(oldCC,ii)
                newCC <- c(newCC,ccc) 
              }
          }
        
        cc <- newCC
      }
    else
      oldCC <- seq(along=cc)
    
    ##-------------------------------------------------------
    ## Define the new nodes (build list of nodes and nodeData)
    ##-------------------------------------------------------
    ndata <- list()
    for(c in cc)
      {

        ## Is there an entrez node in the cc?
        eidx <- which(unlist(sapply(c, FUN=function(s) 'biopax.xref.ENTREZGENE' %in% names(g@nodeData@data[[s]]))))
        if(length(eidx)>1)
          {
            if(separateEntrez)
              warning('One of the formed nodes contains more than one entrez id')
            eidx <- eidx[1]
          }
        ## If so, make sure the corresponding node is the first element of the list (used for naming)
        if(length(eidx) != 0)
          c = c(c[eidx],c[-eidx])        

        ## What properties are present in the nodes of the connected comp
        nppts <- unique(unlist(sapply(c, FUN=function(s) names(g@nodeData@data[[s]]))))
        d <- list()
        ## Merge the nodes of the cc for each property
        for(ppt in nppts)
          d[[ppt]] <- unlist(sapply(c, FUN=function(s) g@nodeData@data[[s]][[ppt]]))

        d$nodeNames <- c

        cname <- d$canonicalName[1]
        
        ## Make sure the name is unique
        ii <- 1
        while(cname %in% names(ndata))
          {
            cname <- paste(d$canonicalName[1],ii,sep='')
            ii <- ii+1
          }
        ndata[[cname]] <- d
      }

    ##-----------------------------------------------------------------
    ## Deal with edges between the new nodes (build edgeL and edgeData)
    ##-----------------------------------------------------------------    
    eL <- as.list(rep(list(edges=NULL),length(cc))) ##vector("list",length(cc))
    names(eL) <- names(ndata)
    edata <- list()
    for(i in seq(along=cc))
      {
        for(j in seq(along=cc)[-i])
          {
            ccjTargets <- nodes(g)[unlist(edgeL(g)[cc[[j]]])]
            if(any(ccjTargets %in% cc[[i]]) && (oldCC[i] != oldCC[j])) ## Any node of cc[[j]] points to any node of cc[[i]]
              {
                ## eList part
                eL[[j]]$edges <- c(eL[[j]]$edges,i)#names(ndata)[i])
                ## edgeData part
                ## Edge ids
                allEdgesL <- sapply(names(ed),FUN=function(s) strsplit(s,'\\|'))
                jiEdges <- unlist(lapply(allEdgesL,FUN=function(e) return((e[1] %in% cc[[j]]) && (e[2] %in% cc[[i]]))))
                ## What properties are present among the original edges
                eppts <- unique(unlist(lapply(ed[jiEdges], FUN=function(e) names(e))))
                d <- list()
                ## Merge the edges for each property
                for(ppt in eppts)
                    d[[ppt]] <- unlist(sapply(ed[jiEdges], FUN=function(e) e[[ppt]]))
                edata[[paste(names(ndata)[j],'|',names(ndata)[i],sep='')]] <- d
              }
          }
      }
    eL <- lapply(eL,FUN=function(e) if(is.null(e$edges)) list(edges=numeric(0)) else e)

    ##--------------------------------------------------
    ## Build the graph
    ##--------------------------------------------------
    pg <- new('graphNEL',names(eL),eL,'directed')

    ##--------------------------------------------------
    ## Set up its edgeData
    ##--------------------------------------------------
    eppts <- unique(unlist(lapply(edata,FUN=function(e) names(e)))) # Existing edgeData properties (must be initialized)
    for(ppt in eppts)
      edgeDataDefaults(pg,ppt) <- ''
    for(e in seq(along=edata))
      {
        
        fromTo <- strsplit(names(edata)[e],'\\|')[[1]]
        for(p in seq(along=edata[[e]]))
          edgeData(pg,fromTo[1],fromTo[2],names(edata[[e]])[p]) <- list(edata[[e]][[p]])
      }

    ##--------------------------------------------------
    ## Set up its nodeData
    ##--------------------------------------------------
    nppts <- unique(unlist(lapply(ndata,FUN=function(e) names(e)))) # Existing nodeData properties (must be initialized)
    for(ppt in nppts)
      nodeDataDefaults(pg,ppt) <- ''
    for(e in seq(along=ndata))
      {        
        for(p in seq(along=ndata[[e]]))
          nodeData(pg,names(ndata)[e],names(ndata[[e]])[p]) <- list(ndata[[e]][[p]])
      }

    if(entrezOnly) # Remove nodes with no entrez id field
      {
        entrezNodes <- nodes(pg)[unlist(lapply(pg@nodeData@data,FUN=function(e) 'biopax.xref.ENTREZGENE' %in% names(e)))]
        pg <- subGraph(entrezNodes,pg)        
      }
    
    return(pg)
  }


############################################################################
# HISTORY:
# 2011-02-15
# o Created.
############################################################################
