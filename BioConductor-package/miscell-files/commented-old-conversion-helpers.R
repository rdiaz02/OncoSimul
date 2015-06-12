## Posets, restriction tables, etc.

## ## When we go poset -> rT or
## ##  adjMat -> rT

## ## we have there all there is to be, or the poset follows the usual
## ## (messy) rules. In particular, CBN does understand that OK. And a linear
## ## poset, for example, places all nodes there. And CBN returns all nodes.


## ## Problem is generating a poset from an adjmat that is the output of OT,
## ## since that adjmat need not contain all nodes, as OT sometimes drops
## ## some of them. So we need to preserve names AND when we write the poset,
## ## we need to specify, as we always do, the correct number of columns.


## ## In newer, future, versions of OncoSimulR all this will be cleared out
## ##    by having posets that are strict. This means:

## ##  - if we pass things to CBN, we will clean that to the messy format in
## ##  write.poset.

## ##  - when using OTs, ??, maybe force all nodes to be there, connected
## ##  from root? Think this later. Or use the current the mess?? Or never
## ##  use the simple poset for nothing, except calling CBN.



## graphToPoset <- function(g) {
##     m <- as(g, "matrix") 
##     mi <- m
##     storage.mode(mi) <- "integer"
##     stopifnot(all.equal(m, mi))
##     return(adjMatToPoset(mi))
## }

## poset.to.restrictTable <- function(x) {
##     ## x1 <- posetToAdjmat(x) ## If I had used the new ones I wrote
##     x1 <- posetToGraph(x, names = 1:max(x), addroot = FALSE, type = "adjmat")
##     adjmat.to.restrictTable(x1)
## }


## checkProperAdjMat <- function(x,
##                               rootNames = c("0", "root", "Root")) {
##     if(is.null(colnames(x)) && is.null(rownames(x)))
##         stop("column and/or row names are null")
##     if(!identical(colnames(x), rownames(x)))
##         stop("colnames and rownames not identical")
##     posRoot <- which(colnames(x) %in% rootNames)
##     if(!length(posRoot)) ## adjmat.to.restrictTable depends on this being true
##         stop("No column with the root name")
##     if(length(posRoot) > 1)
##         stop("Ambiguous location of root")
##     if(!is.integer(x))
##         stop("Non-integer values")
##     if( !all(x %in% c(0, 1) ))
##         stop("Values not in [0, 1]")
##     if( posRoot != 1)
##         stop("Root must be in first row and column")
##     ## The following is not reasonable
##     ## scn <- sort(as.integer(colnames(x)[-1]))
##     ## if(!identical(scn, seq.int(ncol(x) - 1)))
##     ##     stop("Either non-integer column names, or non-successive integers")
## }


## adjmat.to.restrictTable <- function(x, root = FALSE,
##                                     rootNames = c("0", "root", "Root")) {
##     ## we have the zero
##     if( any(colnames(x) %in% c("0", "root", "Root")) & !root)
##         warning("Looks like the matrix has a root but you specified root = FALSE")

##     if(!identical(colnames(x), rownames(x)))
##         stop("colnames and rownames not identical")
##     if(root) {
##         posRoot <- which(colnames(x) %in% rootNames)
##         if(!length(posRoot))
##             stop("No column with the root name")
##         if(length(posRoot) > 1)
##             stop("Ambiguous location of root")
##         x <- x[-posRoot, -posRoot]
##     }

##     if(typeof(x) != "integer")
##         warning("This is not an _integer_ adjacency matrix")
##     if( !all(x %in% c(0, 1) ))
##         stop("Values not in [0, 1]")
##     if(!is.null(colnames(x))) {
##         ## FIXME: this makes sense with numeric labels for columns, but
##         ## not ow.
##         oi <- order(as.numeric(colnames(x)))
##         if(any(oi != (1:ncol(x)))) {
##             warning("Reordering adjacency matrix")
##             x <- x[oi, oi]
##         }
##     }
    
##     num.deps <- colSums(x)
##     max.n.deps <- max(num.deps)
##     rt <- matrix(-9L, nrow = nrow(x),
##                  ncol = max.n.deps + 2)
##     for(i in 1:ncol(x)) {
##         if( num.deps[ i ])
##             rt[i , 1:(2 + num.deps[ i ])] <- c(i, num.deps[i ],
##                                                which(x[, i ] != 0))
##         else
##             rt[i , 1:2] <- c(i , 0L)
##     }
##     class(rt) <- "restrictionTable"
##     return(rt)
## }

## ## Attempt at new version, but does not work, as we need to be able to
## ## deal with matrices without root. Because we often go from a poset.  So
## ## when we do poset -> restrictTable, from the poset we generate an
## ## adjmat wthout root

## ## adjmat.to.restrictTable <- function(x,
## ##                                     rootNames = c("0")) {
    
## ##     null <- checkProperAdjMat(x, rootNames = rootNames)

## ##     ## For now, we only accept matrices with integer column names
## ##     if(!is.integer(type.convert(colnames(x), as.is = TRUE)))
## ##         stop("For now, we only convert to restrictTable if integer row/col. names")
    
## ##     ## restrictTable has no root, ever, so remove it from adj matrix
## ##     ## as adjacency matrix has root, always
## ##     x <- x[-1, -1]
    
## ##     oi <- order(as.numeric(colnames(x)))
## ##     if(any(oi != (1:ncol(x)))) {
## ##         warning("Reordering adjacency matrix")
## ##         x <- x[oi, oi]
## ##     }
    
## ##     num.deps <- colSums(x)
## ##     max.n.deps <- max(num.deps)
## ##     rt <- matrix(-9L, nrow = nrow(x),
## ##                  ncol = max.n.deps + 2)
## ##     for(i in 1:ncol(x)) {
## ##         if( num.deps[ i ])
## ##             rt[i , 1:(2 + num.deps[ i ])] <- c(i, num.deps[i ],
## ##                                                which(x[, i ] != 0))
## ##         else
## ##             rt[i , 1:2] <- c(i , 0L)
## ##     }
## ##     class(rt) <- "restrictionTable"
## ##     return(rt)
## ## }


## adjMatToPoset <- function(x, rootNames = c("0", "root", "Root"),
##                           dropRoot = FALSE) {
##     null <- checkProperAdjMat(x, rootNames)
##     ## rootNames is not used for anything else. We simply check it is
##     ## there, so that the poset will have the root.

    
##     p1 <- intAdjMatToPoset(x, dropRoot = dropRoot)

##     return(p1)
    
##     ## ## Well, nice, but posets now only accept integer labels
##     ##  And if back, we will need to think the dropRoot argument here.
##     ## namNodes <- colnames(x)
##     ## ## Map back to non-integer labels if any used in the adjacency matrix
##     ## namesInts <- is.integer(type.convert(namNodes, as.is = TRUE))
    
##     ## Actually, this is wrong, as you can have nameInts but non
##     ## consecutive indices in the adjMat, and you will get a poset that
##     ## starts at 1, always, etc.

##     ## if(namesInts) {
##     ##     return(p1)
##     ## } else {
##     ##     p2 <- cbind(namNodes[p1[, 1] + 1], namNodes[p1[, 2] + 1])
##     ##     return(p2)
##     ## }
## }

## intAdjMatToPoset <- function(x, dropRoot = FALSE) {
##     if(dropRoot) {
##         ncx <- ncol(x)
##         x <- x[-1, -1]
##         y <- (which(x == 1L, arr.ind = TRUE) )
##         if(nrow(y) == 0) ## all nodes descend from 0
##             y <- cbind(0L, ncx - 1L)

##     } else {
##         y <- (which(x == 1L, arr.ind = TRUE) - 1L)
##     }
##     rownames(y) <- colnames(y) <- NULL
##     storage.mode(y) <- "integer"
##     return(y)
## }



## intAdjMatToPosetPreserveNames <- function(x, dropRoot = FALSE) {
##     if(dropRoot) {
##         ncx <- ncol(x)
##         x <- x[-1, -1]
##         ## y <- (which(x == 1L, arr.ind = TRUE) )
##         ## if(nrow(y) == 0) ## all nodes descend from 0
##         ##     y <- cbind(0L, ncx - 1L)
##         namesInts <- type.convert(colnames(x), as.is = TRUE)
        
##     } else {
##         ## y <- (which(x == 1L, arr.ind = TRUE) )
##         namesInts <- c(0L, type.convert(colnames(x)[-1], as.is = TRUE))
##     }
##     if(!is.integer(namesInts))
##         stop("cannot convert to poset adj mat with non-int colnames")
##     y <- (which(x == 1L, arr.ind = TRUE) )
##     if(nrow(y) == 0) ## all nodes descend from 0
##         y <- cbind(0L, ncx - 1L)

##     p2 <- cbind(namesInts[ y[, 1] ], namesInts[ y[, 2] ])
##     storage.mode(p2) <- "integer"
##     return(p2)
## }    



## convertRestrictTable <- function(x) {
##     ## Only used to put the restrictTable in the format that the C++ code
##     ## uses.
##     t.restrictTable <- matrix(as.integer(x),
##                               ncol = nrow(x), byrow = TRUE)

##     t.restrictTable[-2, ] <- t.restrictTable[-2, ] - 1
##     return(t.restrictTable)
## }


## posetToGraph <- function(x, names,
##                          addroot = FALSE,
##                          type = "graphNEL") {

##     ## Using addroot = FALSE returns adjacency matrices without root,
##     ## which is not neat.
    
##     ## But this is a complicated mess that can lead to confusion. And
##     ## cannot easily handle posets with names of nodes. So I fix all that.
    
##     ## Intermediate nodes, if no ancestor or descendant, need not
##     ## be in the poset.
##     ## Any node with index largest than any node with ancestor or descendant
##     ## needs to be in the file.

##     ## E.g., rbind(c(1,2), c(4, 5)) will get three as a child-less and
##     ## parent-less node.  But to place 6 you need c(6, NA)

    
##     ## We could use something like in run.oncotree,
##     ## as a poset also is a set of parent-child,
##     ## and we would then need to add the root connections
##     ## as is done below in no.ancestor.

##     ## But we do not for now. Note we show lonely nodes, which oncotrees
##     ## do not.  wait: when using root, we do not have "lonely nodes"
##     ## anymore.  But that is irrelevant for metrics based on transitive
##     ## closure. Not for Diff, etc.

##     ## In fact, this is all OK, but is confussing, because I can
##     ## have two kinds of posets: ones that are full, with NAs, etc, if
##     ## needed. Others that are not, but are fixed by the code here. But if
##     ## using the later, the user needs to make sure that the last node is
##     ## in the poset. This can be used as a shortcut trick, but in the docs
##     ## I do not do it, as it is bad practice.

##     if(!is.integer(x))
##         stop("poset contains non-integers")
    
##     m <- length(names) 
##     m1 <- matrix(0L, nrow = m, ncol = m)
##     colnames(m1) <- names
##     rownames(m1) <- names
##     if(is.null(dim(x)) ) {
##         if(length(x) != 1)
##             stop("If poset is not a matrix, it must be a vector of length 1")
##     } else { ## a matrix so a non-null poset
##         ## a debugging check. 
##         if(nrow(x) > 0) {
##             if(addroot) {
##                 if(max(x, na.rm = TRUE) != (m - 1))
##                     stop("\n in poset, max(x) != (m - 1)")
##             } else { ## this was missing!
##                 if(max(x, na.rm = TRUE) != m)
##                     stop("\n in poset, max(x) != m")
##             }
##         }
##         if(nrow(x) > 0) {
##             if(addroot)
##                 m1[x + 1] <- 1L
##             else
##                 m1[x] <- 1L ## this will remove all entries with a 0
##                             ## index. So posets where explicit the dep. on
##                             ## 0.
##         }
##         if((length(names) > 1) & addroot) {
##             no.ancestor <- which(apply(m1, 2, function(x) all(x == 0)))
##             no.ancestor <- no.ancestor[-1]
##             m1[cbind(1, no.ancestor)] <- 1L
##         } ## o.w. do nothing
##     }
##     if(addroot)
##         m1[1, 1] <- 0L
    
##     if(type == "adjmat") return(m1)
##     else if (type == "graphNEL") return(as(m1, "graphNEL"))
##     ## does not show the labels
##     ## else if (type == "igraph") return(graph.adjacency(m1))
## }


## ## adjmat.to.restrictTableOld <- function(x, root = FALSE,
## ##                                     rootNames = c("0", "root", "Root")) {
## ##     ## we have the zero
## ##     if( any(colnames(x) %in% c("0", "root", "Root")) & !root)
## ##         warning("Looks like the matrix has a root but you specified root = FALSE")

## ##     if(!identical(colnames(x), rownames(x)))
## ##         stop("colnames and rownames not identical")
## ##     if(root) {
## ##         posRoot <- which(colnames(x) %in% rootNames)
## ##         if(!length(posRoot))
## ##             stop("No column with the root name")
## ##         if(length(posRoot) > 1)
## ##             stop("Ambiguous location of root")
## ##         x <- x[-posRoot, -posRoot]
## ##     }

## ##     if(typeof(x) != "integer")
## ##         warning("This is not an _integer_ adjacency matrix")
## ##     if( !all(x %in% c(0, 1) ))
## ##         stop("Values not in [0, 1]")
## ##     if(!is.null(colnames(x))) {
## ##         ## FIXME: this makes sense with numeric labels for columns, but
## ##         ## not ow.
## ##         oi <- order(as.numeric(colnames(x)))
## ##         if(any(oi != (1:ncol(x)))) {
## ##             warning("Reordering adjacency matrix")
## ##             x <- x[oi, oi]
## ##         }
## ##     }
    
## ##     num.deps <- colSums(x)
## ##     max.n.deps <- max(num.deps)
## ##     rt <- matrix(-9L, nrow = nrow(x),
## ##                  ncol = max.n.deps + 2)
## ##     for(i in 1:ncol(x)) {
## ##         if( num.deps[ i ])
## ##             rt[i , 1:(2 + num.deps[ i ])] <- c(i, num.deps[i ],
## ##                                                which(x[, i ] != 0))
## ##         else
## ##             rt[i , 1:2] <- c(i , 0L)
## ##     }
## ##     class(rt) <- "restrictionTable"
## ##     return(rt)
## ## }







## ## ## names are preserved in igraph, graphNEL, interconversions, etc
## ## uu <- matrix(c(1, 0, 0, 1), ncol = 2)
## ## colnames(uu) <- rownames(uu) <- c("root", "A")
## ## get.adjacency(graph.adjacency(uu), sparse = FALSE)


## ## get.adjacency(igraph.from.graphNEL(as(uu, "graphNEL")), sparse = FALSE)
## ## as(igraph.to.graphNEL(graph.adjacency(uu)), "matrix")



## ## ## These look nice, etc, but I am not really using them that much, and I
## ## ## am not exposing them to the user. Thus, return to former things.
## ## ### z1: these two do not use keeproot. It is always true.
## ## posetToAdjmatNew <- function(x, keeproot = TRUE, rootNames = c("0",
## ## "root", "Root")) { internalPosetToGraph(x, type = "adjmat", keeproot =
## ## keeproot, rootNames = rootNames) }

## ## posetToGraphNew <- function(x,
## ##                          keeproot = TRUE,
## ##                          rootNames = c("0", "root", "Root")) {
## ##     internalPosetToGraph(x, type = "graphNEL", keeproot = keeproot,
## ##                  rootNames = rootNames)
## ## }

## ## checkProperPoset <- function(x, rootNames = c("0", "root", "Root")) {
## ##     ## checks if proper poset. If it is, returns poset and sorted unique
## ##     ## names

## ##     ## I return things to avoid doing the same work twice.
## ##     if(ncol(x) != 2)
## ##         stop("Posets must have two columns")
## ##     nn <- unique(as.vector(x))
## ##     matchRoot <- which(nn %in% rootNames)
## ##     if(length(matchRoot) > 1)
## ##         stop("More than one possible root")
## ##     if(length(matchRoot) == 0)
## ##         stop("No root in the poset")
## ##     n1 <- nn[matchRoot]
## ##     nn <- c(n1, nn[-matchRoot])
## ##     return(list(x, uniqueNames = nn))
## ## }

## ## internalPosetToGraph <- function(x,
## ##                                  type,
## ##                                  keeproot = TRUE,
## ##                                  rootNames = c("0", "root", "Root"),
## ##                                  addroot = FALSE,
## ##                                  names = NULL
## ##                                  ) {

## ##     ## So there is a simple, straightforward version, but we need to
## ##     ## account for posets that include no root. That one is called when we
## ##     ## issue addroot = TRUE, and we must then pass a valid names vector.
    
    
    
## ##     pp <- checkProperPoset(x, rootNames = rootNames)
    
## ##     m1 <- matrix(0L, nrow = length(pp$uniqueNames),
## ##                  ncol = length(pp$uniqueNames))
## ##     colnames(m1) <- pp$uniqueNames
## ##     rownames(m1) <- pp$uniqueNames

## ##     m1[pp$x] <- 1L

## ##     if(!keeproot)
## ##         m1 <- m1[-1, -1]
## ##     if(type == "adjmat") return(m1)
## ##     else if (type == "graphNEL") return(as(m1, "graphNEL"))
## ##     ## else if (type == "igraph") return(graph.adjacency(m1))
## ##     ## I am not using igraph now
## ## }







## ## this is blatantly wrong
## ## graph.to.poset <- function(x) {
## ##   ## FIXME: this are characters, not numeric
## ##   return(matrix(as.numeric(unlist(edgeL(x))), ncol = 2,
## ##                 byrow = TRUE))
## ## }





## ## why was this working? when converting adjmat to poset to adjmat?
## m1 <- structure(c(0L, 0L, 0L, 1L, 0L, 0L, 0L, 1L, 0L), .Dim = c(3L, 
## 3L), .Dimnames = list(c("Root", "2", "4"), c("Root", "2", "4"
##                                              )))

## m1b <- structure(c(0L, 0L, 0L, 0L, 0L, 1L, 1L, 0L, 0L), .Dim = c(3L, 
## 3L), .Dimnames = list(c("Root", "2", "4"), c("Root", "2", "4"
## )))
