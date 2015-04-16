## this is blatantly wrong
## graph.to.poset <- function(x) {
##   ## FIXME: this are characters, not numeric
##   return(matrix(as.numeric(unlist(edgeL(x))), ncol = 2,
##                 byrow = TRUE))
## }

graphToPoset <- function(g) {
    return(intAdjMatToPoset(as(g, "matrix")))
}

poset.to.restrictTable <- function(x) {
    x1 <- posetToGraph(x, names = 1:max(x), addroot = FALSE, type = "adjmat")
    adjmat.to.restrictTable(x1)
}


## have a graph without a null??

checkProperPoset <- function(x, rootNames = c("0", "root", "Root")) {
    ## checks if proper poset. If it is, returns poset and sorted unique
    ## names
    if(ncol(x) != 2)
        stop("Posets must have two columns")
    nn <- unique(as.vector(x))
    matchRoot <- which(nn %in% rootNames)
    if(length(matchRoot) > 1)
        stop("More than one possible root")
    if(length(matchRoot) == 0)
        stop("No root in the poset")
    n1 <- nn[matchRoot]
    nn <- c(n1, nn[-matchRoot])
    return(list(x, uniqueNames = nn))
}


checkProperAdjMat <- function(x) {

}

posetToGraph <- function(x,
                         type = "graphNEL",
                         keeproot = TRUE,
                         rootNames = c("0", "root", "Root")) {

    pp <- checkProperPoset(x)
    
    m1 <- matrix(0L, nrow = length(pp$uniqueNames),
                 ncol = length(pp$uniqueNames))
    colnames(m1) <- pp$uniqueNames
    rownames(m1) <- pp$uniqueNames

    m1[x] <- 1L

    if(!keeproot)
        m1 <- m1[-1, -1]
    if(type == "adjmat") return(m1)
    else if (type == "graphNEL") return(as(m1, "graphNEL"))
    ## does not show the labels
    ## else if (type == "igraph") return(graph.adjacency(m1))
}


posetToGraphOld <- function(x, names,
                            addroot = FALSE,
                            type = "graphNEL") {

    ## This used to be the old code. But this is a complicated mess that
    ## can lead to confusion. And cannot easily handle posets with names
    ## of nodes. So I fix all that.
    
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
    
    if(type == "adjmat") return(m1)
    else if (type == "graphNEL") return(as(m1, "graphNEL"))
    ## does not show the labels
    ## else if (type == "igraph") return(graph.adjacency(m1))
}

adjmat.to.restrictTable <- function(x, root = FALSE,
                                    rootNames = c("0", "root", "Root")) {
    ## we have the zero
    if( any(colnames(x) %in% c("0", "root", "Root")) & !root)
        warning("Looks like the matrix has a root but you specified root = FALSE")

    if(!identical(colnames(x), rownames(x)))
        stop("colnames and rownames not identical")
    if(root) {
        posRoot <- which(colnames(x) %in% rootNames)
        if(!length(posRoot))
            stop("No column with the root name")
        if(length(posRoot) > 1)
            stop("Ambiguous location of root")
        x <- x[-posRoot, -posRoot]
    }

    if(typeof(x) != "integer")
        warning("This is not an _integer_ adjacency matrix")
    if( !all(x %in% c(0, 1) ))
        stop("Values not in [0, 1]")
    if(!is.null(colnames(x))) {
        ## FIXME: this makes sense with numeric labels for columns, but
        ## not ow.
        oi <- order(as.numeric(colnames(x)))
        if(any(oi != (1:ncol(x)))) {
            warning("Reordering adjacency matrix")
            x <- x[oi, oi]
        }
    }
    
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


adjMatToPoset <- function(x, rootNames = c("0", "root", "Root")) {
    p1 <- intAdjMatToPoset(x, rootNames)
    namNodes <- colnames(x)
    ## Map back to non-integer labels if any used in the adjacency matrix
    namesInts <- is.integer(type.convert(namNodes, as.is = TRUE))
    if(!namesInts) {    
        p <- cbind(namNodes[p1[, 1] + 1], namNodes[p1[, 2] + 1])
        return(p)
    } else {
        return(p1)
    }
}

intAdjMatToPoset <- function(x, 
                             rootNames = c("0", "root", "Root")) {

     if(!identical(colnames(x), rownames(x)))
         stop("colnames and rownames not identical")
     posRoot <- which(colnames(x) %in% rootNames)
     if(!length(posRoot))
         stop("No column with the root name")
     if(length(posRoot) > 1)
         stop("Ambiguous location of root")
     y <- (which(x == 1, arr.ind = TRUE) - 1L)
     rownames(y) <- colnames(y) <- NULL
     return(y)
}    



convertRestrictTable <- function(x) {
    ## Only used to put the restrictTable in the format that the C++ code
    ## uses.
    t.restrictTable <- matrix(as.integer(x),
                              ncol = nrow(x), byrow = TRUE)

    t.restrictTable[-2, ] <- t.restrictTable[-2, ] - 1
    return(t.restrictTable)
}

