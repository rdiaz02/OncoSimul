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


## Posets, restriction tables, etc.


## This is all too complicated. Have here just the minimal set we need,
## because we will soon be changing formats.



## When we go poset -> rT or
##  adjMat -> rT

## we have there all there is to be, or the poset follows the usual
## (messy) rules. In particular, CBN does understand that OK. And a linear
## poset, for example, places all nodes there. And CBN returns all nodes.


## Problem is generating a poset from an adjmat that is the output of OT,
## since that adjmat need not contain all nodes, as OT sometimes drops
## some of them. So we need to preserve names AND when we write the poset,
## we need to specify, as we always do, the correct number of columns.
## Thus, we use the special intAdjMatToPosetPreserveNames.
## Change name! OTtoPoset


## In newer, future, versions of OncoSimulR all this will be cleared out
##    by having posets that are strict. This means:

##  - if we pass things to CBN, we will clean that to the messy format in
##  write.poset.

## And we will relax adjMats, to allow for non-integer names? Or pass that
## as additional labels, so we can interchange with posets without ambiguity?

## adjMatToPoset is not really used for anything at all in the code now.

## Remember from adjmat we go to rT directly.

## Now, posetToGraph also includes checking strict adjMat, except if
## called from plotPoset, where arbitrary labels can be used.





checkProperMinimalAdjMat <- function(x,
                                     rootNames = c("0", "root", "Root"),
                                     root = TRUE,
                                     orderedNames = TRUE) {

    ## We go from poset to rt through adjmat, but the adjmat has no root,
    ## as the rT has no root. But this must have ordered names, as any
    ## poset uses only integers from 1 to whatever.

    ## Any adjmat should fulfill what follows, with the option of not
    ## checking for root and for orderedNames
    
    if(is.null(colnames(x)) && is.null(rownames(x)))
        stop("column and/or row names are null")
    if(!identical(colnames(x), rownames(x)))
        stop("colnames and rownames not identical")
    if(!is.integer(x))
        stop("Non-integer values")
    if( !all(x %in% c(0L, 1L) ))
        stop("Values not in {0L, 1L}")
    posRoot <- which(colnames(x) %in% rootNames)
    cs <- colSums(x)
    if(root) {
        ## only applies to rooted adj mat, as a matrix without root where
        ## all descend from root will have not a single 1
        if(any( (cs + rowSums(x)) < 1))
            stop("You have nodes that have no incoming and no outgoing connections")
    }
    if(length(posRoot) > 0) {
        if(length(posRoot) > 1)
            stop("Ambiguous location of root")
        if( posRoot != 1)
            stop("Root must be in first row and column")
        if(cs[1] > 0)
            stop("Nothing can have an edge to Root")
        cs <- cs[-1]
        if(any(cs < 1))
            stop(paste("In trees with Root,",
                       " there must be at least one incoming edge to every node (except Root)"))
        if(!root)
            stop("You said there is no root, but there is one")
    } else {
        if(root)
            stop("No Root, and you said one should be present")
    }

    if(orderedNames) {
        if(length(posRoot) == 1) {
            ncx <- ncol(x) - 1
            scn <- sort(as.integer(colnames(x)[-1]))
        }
        else {
            ncx <- ncol(x)
            scn <- sort(as.integer(colnames(x)))
        }
        if(!identical(scn, seq.int(ncx)))
            stop(paste("Either non-integer column names,",
                       " or non-successive integers, ",
                       " or not sorted integer names") )
    }
}


checkProperFullAdjMat <- function(x,
                              rootNames = c("0", "root", "Root")) {
    checkProperMinimalAdjMat(x, rootNames = rootNames,
                             root = TRUE, orderedNames = TRUE)
}


checkProperOTAdjMat <- function(x) {
    ## oncotree ALWAYS returns things that have a rootName that says "Root"
    checkProperMinimalAdjMat(x, rootNames = "Root", root = TRUE,
                             orderedNames = FALSE)
}



poset.to.restrictTable <- function(x) {
    x1 <- posetToGraph(x, names = 1:max(x), addroot = FALSE, type = "adjmat")
    adjmat.to.restrictTable(x1, root = FALSE,
                            rootNames = c("0", "root", "Root"))
}


adjmat.to.restrictTable <- function(x, root = FALSE,
                                    rootNames = c("0", "root", "Root")) {

    checkProperMinimalAdjMat(x, rootNames = rootNames,
                             root = root,
                             orderedNames = TRUE)
    if(root)
        x <- x[-1, -1, drop = FALSE]
    ## ## we have the zero
    ## if( any(colnames(x) %in% c("0", "root", "Root")) & !root)
    ##     warning("Looks like the matrix has a root but you specified root = FALSE")

    ## if(!identical(colnames(x), rownames(x)))
    ##     stop("colnames and rownames not identical")
    ## if(root) {
    ##     posRoot <- which(colnames(x) %in% rootNames)
    ##     if(!length(posRoot))
    ##         stop("No column with the root name")
    ##     if(length(posRoot) > 1)
    ##         stop("Ambiguous location of root")
    ##     x <- x[-posRoot, -posRoot]
    ## }

    ## if(typeof(x) != "integer")
    ##     warning("This is not an _integer_ adjacency matrix")
    ## if( !all(x %in% c(0, 1) ))
    ##     stop("Values not in [0, 1]")

    ## if(!is.null(colnames(x))) {
    ##     ## FIXME: this makes sense with numeric labels for columns, but
    ##     ## not ow.
    ##     oi <- order(as.numeric(colnames(x)))
    ##     if(any(oi != (1:ncol(x)))) {
    ##         warning("Reordering adjacency matrix")
    ##         x <- x[oi, oi]
    ##     }
    ## }
    
    num.deps <- colSums(x)
    max.n.deps <- max(num.deps)
    rt <- matrix(-9L, nrow = nrow(x),
                 ncol = max.n.deps + 2)
    for(i in 1:ncol(x)) {
        if( num.deps[ i ])
            rt[i , 1:(2 + num.deps[ i ])] <- c(i, num.deps[i ],
                                               which(x[, i ] != 0))
        else
            rt[i , 1:2] <- c(i , 0L)
    }
    class(rt) <- "restrictionTable"
    return(rt)
}




## this preserves the names in the adj mat. For OT -> poset.
## intAdjMatToPosetPreserveNames <- function(x, dropRoot = TRUE) {

OTtoPoset <- function(x) {
    checkProperOTAdjMat(x)

    ## ## root must always be dropped, as we are creating a poset
    ## dropRoot <- TRUE
    ## if(!dropRoot)
    ##     warning("Are you sure you do not want dropRoot?")
    ## if(dropRoot) {
    ##     ncx <- ncol(x)
    ##     x <- x[-1, -1]
    ##     ## y <- (which(x == 1L, arr.ind = TRUE) )
    ##     ## if(nrow(y) == 0) ## all nodes descend from 0
    ##     ##     y <- cbind(0L, ncx - 1L)
    ##     namesInts <- type.convert(colnames(x), as.is = TRUE)
        
    ## } else {
    ##     ## y <- (which(x == 1L, arr.ind = TRUE) )
    ##     namesInts <- c(0L, type.convert(colnames(x)[-1], as.is = TRUE))
    ## }

    ncx <- ncol(x)
    x <- x[-1, -1, drop = FALSE]
    ## y <- (which(x == 1L, arr.ind = TRUE) )
    ## if(nrow(y) == 0) ## all nodes descend from 0
    ##     y <- cbind(0L, ncx - 1L)
    namesInts <- type.convert(colnames(x), as.is = TRUE)
    
    if(!is.integer(namesInts))
        stop("cannot convert to poset adj mat with non-int colnames")
    y <- (which(x == 1L, arr.ind = TRUE) )
    if(nrow(y) == 0) ## all nodes descend from 0
        y <- cbind(0L, ncx - 1L)

    p2 <- cbind(namesInts[ y[, 1] ], namesInts[ y[, 2] ])
    storage.mode(p2) <- "integer"
    if(ncol(p2) == 1) ## the hack for when all from root
        p2 <- matrix(nrow = 0, ncol = 2)
    return(p2)
}    

## the next is just a convenience
## sortAdjMat <- function(am) {
##     cn <- colnames(am)
##     rootpos <- grep("^Root$", cn) 
##     if(length(rootpos) != 1)
##         stop("No root in adj mat, or multiple Roots")
##     cn <- c("Root", sort(colnames(am)[-rootpos]))
##     return(am[cn, cn])
## }


## No longer used
## sortAdjMat <- function(am) {
##     ## If column names, except Root, are integers, sort as integers. O.w.,
##     ## general lexicog. sort.
##     cn <- colnames(am)
##     rootpos <- grep("^Root$", cn) 
##     if(length(rootpos) != 1)
##         stop("No root in adj mat, or multiple Roots")
##     cn0 <- colnames(am)[-rootpos]
##     namesInts <- type.convert(cn0, as.is = TRUE)
##     if(is.integer(namesInts)) {
##         cn <- c("Root", sort(namesInts))
##     } else {
##         cn <- c("Root", sort(cn0))
##     }
##     return(am[cn, cn])
## }




convertRestrictTable <- function(x) {
    ## Only used to put the restrictTable in the format that the C++ code
    ## uses.
    t.restrictTable <- matrix(as.integer(x),
                              ncol = nrow(x), byrow = TRUE)

    t.restrictTable[-2, ] <- t.restrictTable[-2, ] - 1
    return(t.restrictTable)
}


posetToGraph <- function(x, names,
                         addroot = FALSE,
                         type = "graphNEL",
                         strictAdjMat = TRUE) {

    ## Using addroot = FALSE returns adjacency matrices without root,
    ## which is not neat.
    
    ## But this is a complicated mess that can lead to confusion. And
    ## cannot easily handle posets with names of nodes. So I fix all that.
    
    ## Intermediate nodes, if no ancestor or descendant, need not
    ## be in the poset.
    ## Any node with index largest than any node with ancestor or descendant
    ## needs to be in the file.

    ## E.g., rbind(c(1,2), c(4, 5)) will get three as a child-less and
    ## parent-less node.  But to place 6 you need c(6, NA)

    
    ## We could use something like in run.oncotree,
    ## as a poset also is a set of parent-child,
    ## and we would then need to add the root connections
    ## as is done below in no.ancestor.

    ## But we do not for now. Note we show lonely nodes, which oncotrees
    ## do not.  wait: when using root, we do not have "lonely nodes"
    ## anymore.  But that is irrelevant for metrics based on transitive
    ## closure. Not for Diff, etc.

    ## In fact, this is all OK, but is confussing, because I can
    ## have two kinds of posets: ones that are full, with NAs, etc, if
    ## needed. Others that are not, but are fixed by the code here. But if
    ## using the later, the user needs to make sure that the last node is
    ## in the poset. This can be used as a shortcut trick, but in the docs
    ## I do not do it, as it is bad practice.

    if(!is.integer(x))
        stop("poset contains non-integers")
    
    m <- length(names) 
    m1 <- matrix(0L, nrow = m, ncol = m)
    colnames(m1) <- names
    rownames(m1) <- names
    if(is.null(dim(x)) ) {
        if(length(x) != 1)
            stop("If poset is not a matrix, it must be a vector of length 1")
    } else { ## a matrix so a non-null poset
        ## a debugging check. 
        if(nrow(x) > 0) {
            if(addroot) {
                if(max(x, na.rm = TRUE) != (m - 1))
                    stop("\n in poset, max(x) != (m - 1)")
            } else { ## this was missing!
                if(max(x, na.rm = TRUE) != m)
                    stop("\n in poset, max(x) != m")
            }
        }
        if(nrow(x) > 0) {
            if(addroot)
                m1[x + 1] <- 1L
            else
                m1[x] <- 1L ## this will remove all entries with a 0
                            ## index. So posets where explicit the dep. on
                            ## 0.
        }
        if((length(names) > 1) & addroot) {
            no.ancestor <- which(apply(m1, 2, function(x) all(x == 0)))
            no.ancestor <- no.ancestor[-1]
            m1[cbind(1, no.ancestor)] <- 1L
        } ## o.w. do nothing
    }
    if(addroot)
        m1[1, 1] <- 0L
    ## yes, we must skip this check when plotting a poset (because we can
    ## have arbitrary non-integer names )
    if(strictAdjMat)
        checkProperMinimalAdjMat(m1, root = addroot)
    
    if(type == "adjmat") return(m1)
    else if (type == "graphNEL") return(as(m1, "graphNEL"))
    ## does not show the labels
    ## else if (type == "igraph") return(graph.adjacency(m1))
}

