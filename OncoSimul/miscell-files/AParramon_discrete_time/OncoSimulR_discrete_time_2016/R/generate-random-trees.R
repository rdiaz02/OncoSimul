## Copyright 2013, 2014, 2015 Ramon Diaz-Uriarte

## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.



simOGraph <- function(n, h = 4, conjunction = TRUE, nparents = 3,
                      multilevelParent = TRUE,
                      removeDirectIndirect = TRUE,
                      rootName = "Root") {
    ## Returns an adjacency matrix
    if(h > n)
        stop("h > n")

    if( (conjunction == FALSE) & (nparents > 1) )
        nparents <- 1
    if( nparents == 1) {
        ## indirect issues do not affect if nparents == 1
        removeDirectIndirect <- FALSE
    }
    if(h == 1) { ## obviously, since all just descend from root
        multilevelParent <- FALSE
        removeDirectIndirect <- FALSE
    }

    adjMat <- matrix(0L, ncol = n + 1, nrow = n + 1)
    ## split into h groups

    ## each element in a level can be connected to one or more (if
    ## conjunction) elements in a previous level

    if(h > 1) {
        grs <- mysplitb(n, h)
        ## Most papers do not show both direct and indirect dependencies (but
        ## Farahani does), and most connect nodes with the immediately
        ## superior level, not to several levels. Note this thing of direct
        ## and indirect not an issue if no conjunctions allowed. Semantics if
        ## both do differ for Monotone and semimonotone, of course, if both
        ## direct and indirect.

        ## All, except first group --- those directly connected to root
        for(i in (2:h)) {
            if(multilevelParent) {
                parents <- unlist(grs[1:(i - 1)])
            } else {
                parents <- grs[[i - 1]]
            }
            adjMat[ , grs[[i]] + 1 ] <- connect(parents, grs[[i]], conjunction,
                                                 nparents, n) 
        }
        ## Those in the first group, the ones connected to root
        adjMat[1, grs[[1]] + 1 ] <- 1L
    } else { ## h = 1, so all connected to root
        adjMat[1, 2:(n+1) ] <- 1L
    }
    
    colnames(adjMat) <- rownames(adjMat) <- c(rootName, 1:n)
    

    ## Prune to remove indirect connections
    if(multilevelParent & removeDirectIndirect)
        adjMat <- removeIndirectConnections(adjMat)
    return(adjMat)
}

## simOG <- simOGraph

mysplitb <- function(n, gr) {
    ## Split nodes (n) in gr sets.

    if(gr == 1)
        return(list(seq.int(n)))
    if(gr == n)
        return(as.list(seq.int(n)))
    else {
        breaks <- sort(sample(seq.int(n - 1), gr - 1))
        dd <- c(breaks[1], diff(breaks), n - breaks[gr - 1])
        split(seq.int(n), factor(rep(seq.int(gr), dd)))
    }
}


connect <- function(parents, descendants, conjunction, nparents, n) {
    ## Returns the columns of the adjacency matrix
    
    ## How many parents will each descendant have?
    if(conjunction) {
        np <- replicate(length(descendants),
                        sample(seq.int(nparents), 1))
        np <- vapply(np, function(x) min(x, length(parents)), numeric(1))
    }
    else
        np <- rep(1, length(descendants))
 
    ## ml <- Map(function(i, np) connectIndiv(indiv = i, parents, np = np, n),
    ##           descendants, np)
    ## a loop is a lot simpler
    mret <- matrix(0L, nrow = n + 1, ncol = length(descendants))
    for(i in seq_along(descendants)) {
        mret[, i] <- connectIndiv(parents, np[i], n)
    }
    mret
}



connectIndiv <- function(parents, nparents, n) {
    if(length(parents) > 1) 
        thep <- sample(parents, nparents, replace = FALSE)
    else
        thep <- parents
    v <- rep(0L, n)
    v[thep] <- 1L
    return(c(0L,v)) ## added root
}

## Not used
## findSuperParents <- function(x, adjMat) {
##     parents <- which(adjMat[, x + 1]  == 1) - 1
##     allP <- findAllParents(x, adjMat)
##     return(setdiff(allP, parents))
## }

findAllParents <- function(x, adjMat) {
    if(x == 0)
        return(NULL)
    else{
        p <- which(adjMat[, x + 1] == 1) - 1
        p1 <- unlist(lapply(p, function(x) findAllParents(x, adjMat)))
        return(c(p, p1))
    }
}

repeatedParents <- function(x, adjMat) {
    ap <- findAllParents(x, adjMat)
    dups <- duplicated(ap)
    dupP <- setdiff(ap[dups], 0)
    dupP
}


removeIndirectConnections <- function(adjMat) {
    ## This is a bad name: we remove the direct connections. How? We
    ## search, for each node, for the set of all
    ## parents/grandparents/grandgranparents. If any of those ancestors is
    ## repeated, it means you go from that ancestor to the node in
    ## question through at least two different routes. Thus, ensure the
    ## direct is 0 (it might already be, no problem). Once you do that,
    ## you know there are not both indirect AND direct connections.
    for(i in ncol(adjMat):2) {
        dp <- repeatedParents( i - 1, adjMat)
        if(length(dp))
            adjMat[cbind(dp + 1, i)] <- 0L
    }
    return(adjMat)
}
