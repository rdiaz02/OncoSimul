




## just for me; I assume everything starts with rt* and ends in RData
## and the object is called tmp
plotdriversdir0 <- function(...){
  op <- par(ask = TRUE)
  files <- dir(pattern = "^rt.*RData$")
  for(fi in files) {
    load(fi)
    plotDrivers0(tmp, main = fi, ...)
  }
  par(op)
}
## using rds
plotdriversdir <- function(...){.R
  op <- par(ask = TRUE)
  files <- dir(pattern = "^rt.*rds$")
  for(fi in files) {
    tmp <- readRDS(fi)
    plotDrivers0(tmp, main = fi, ...)
  }
  par(op)
}


plotDrivers <- function(z, na.subs = TRUE, log = "y", type = "l",
                        lty = 1:9, col = c(8, "orange", 6:1),
                        lwd = 2, ...) {
  ## we pass only the driver data frame
  y <- z[, 2:ncol(z)]
  if(na.subs){
    y[y == 0] <- NA
  }

  matplot(x = z[, 1],
          y = y,
          type = type, log = log, lty = lty, col = col, lwd = lwd,
          ...)
  ## will need to add a legend
  legend(x = "topleft",
         title = "Number of drivers",
         lty = lty, col = col, lwd = lwd,
         legend = (1:ncol(y)) - 1)
}


sampleZZ <- function(zz, seed = "auto"){
  if(seed == "auto") {
    ## paste as numeric(hostname()) y el segundo from time
  }

  fname <- paste(fileroot, randomstgring)
  save(...., compress = FALSE)

  
}


## FIXME: rename as "selectJointFreq"
selectSB <- function(x, threshold = 5, maxGenes = 12, weighted = TRUE,
                     plot = FALSE, returnData = TRUE) {
    ## x: the genotypes, with columns as genotypes, rows as genes
    ## threshold is percentage here.
    ## use tcrossprod? 
  xx <- x %*% t(x)
  ## ceiling used to ensure threshold at least 1
  threshold <- ceiling((threshold/100) * ncol(x))
  xx[xx < threshold] <- 0
  if(!weighted)
    xx[xx >= threshold] <- 1
  ## diag(xx) <- 0
  if(weighted)
    g1 <- graph.adjacency(xx,
                          diag = FALSE,
                          weighted = TRUE,
                          mode = "undirected")
  else
    g1 <- graph.adjacency(xx,
                          diag = FALSE,
                          weighted = NULL,
                          mode = "undirected")
  if(plot) plot(g1)
  ## print(xx)
  lc <- largest.cliques(g1)
  if(length(lc) > 1) {
    cat("\n WARNING: more than one largest clique\n")
    if(!weighted)
        lcx <- lc[[1]]
    else {
        ## return the clique with largest number of connections
        ## this is arbritrary, of course
        sums.w <- sapply(lc,
                         function(z) sum(xx[z, z]))
        the.lc <- which.max(sums.w)
        lcx <- lc[[the.lc]]
        message("   returning largest clique number ", the.lc)
    }

  } else
      lcx <- lc[[1]]
  
  if(length(lcx) > maxGenes)
    lcx <- lcx[order(diag(xx)[lcx], decreasing = TRUE)][1:maxGenes]
  if(returnData)
      return(x[lcx, , drop = FALSE])
  else
      return(lcx)
}



## rename as "selectMarginalFreq"
selectPercent <- function(x, threshold = 5, maxGenes = 12, returnData = TRUE) {
      ## x: the genotypes, with columns as genotypes, rows as genes
    ## threshold is percentage here.
  xx <- rowSums(x)
  threshold <- ceiling((threshold/100) * ncol(x))
  lcx <- which(xx >= threshold)
  if(length(lcx) > maxGenes)
    lcx <- lcx[order(xx[lcx], decreasing = TRUE)][1:maxGenes]
  if(returnData)
      return(x[lcx, , drop = FALSE])
  else
      return(lcx)
}

selectGenes <- function(x, method, threshold = 5, maxGenes = 12,
                        returnData = TRUE, plot = FALSE) {
  if(method == "SB")
    return(selectSB(x, threshold, maxGenes = maxGenes, returnData, plot))
  if(method == "Percent")
    return(selectPercent(x, threshold, maxGenes = maxGenes, returnData))
}



veltcs <- function(x, remove.root = TRUE) {
  ## return a vector from edge list of transitive closure
  ## and yes, edges used to be in igraph, but now also in graph, and the
  ## graph authors used the same name. Great!

    ## we really want to remove root, since otherwise we would be inflating
    ## the number of connections.
  if(length(x) > 1) {
    if(class(x[[1]]) == "graphNEL")
      tc0 <- lapply(x, transitive.closure)
    else ## adjacency matrix
      tc0 <- lapply(x, transClos)
    el <- unlist(lapply(tc0, graph::edges), recursive = FALSE)

  } else{
    if(class(x) == "graphNEL")
      tc0 <- transitive.closure(x)
    else ## adjacency matrix
      tc0 <- transClos(x)
    el <- graph::edges(tc0)
  }
  if(remove.root) {
      el <- el[-which(names(el) == "Root")]
      if(!length(el))
          warning("edge list of length zero. Probably using a tree without Root")
  }
  ll <- lapply(el, length)
  el <- el[ll > 0]

  ell <- unlist(lapply(seq_along(el),
                       function(i) paste(names(el)[[i]], el[[i]], sep = "_")))

  ## ell <- unlist(sapply(names(sapply(el2, names)),
  ##                      function(x) {paste(x, el2[[x]], sep = "_")}   ))
  
  return(unique(ell)) ## important with mixtures
}



metrics1 <- function(x, y) {
  ## My PFD and PND
  ## x is model
  ## y is true

    ## removing the root seems coherent with what Gerstung et al do
    ## we are trying to discover dependency relationships
  ex <- veltcs(x, remove.root = TRUE)
  ey <- veltcs(y, remove.root = TRUE)
  pfd <- length(setdiff(ex, ey))/length(ex)
  pnd <- length(setdiff(ey, ex))/length(ey)
  return(c("PFD" = pfd, "PND " = pnd))  
}




linear.poset <- function(x) {
  ## a direct translation of linear_poset and write_poset
  ## in cbn.py by
  ## Niko Beerenwinkel and Moritz Gerstung

  nr <- nrow(x)
  nc <- ncol(x)
  sorted <- order(colMeans(x), decreasing = TRUE)
  poset <- matrix(0, ncol = nc, nrow = nc)
  s <- sorted[1]
  for (t in sorted[2:nc]) {
    poset[s, t] <- 1
    s <- t
  }
 
  ## now, translate write_poset.
  ## posetw <- matrix(0, ncol = 2, nrow = nr)
  ## for (i in 1:nc) {
  ##   for(j in 1:nc) {
  ##     if(poset[i, j])
  ##       posetw[i, ] <- c(i, j)
  ##   }
  ## }
  
  ## do the R way
  posetw <- which(poset == 1, arr.ind = TRUE)
  posetw <- posetw[order(posetw[, 1]), ]
  return(posetw)
}

write.linear.poset <- function(x, filename) {
  poset <- linear.poset(x)
  fn <- paste(filename, ".poset", sep = "")
  write(ncol(x), file = fn)
  write(t(poset), file = fn, append = TRUE,
        sep = " ", ncolumns = 2)
  write(0, file = fn, append = TRUE)
 
}

graph.to.poset <- function(x) {
  ## FIXME: this are characters, not numeric
  return(matrix(as.numeric(unlist(edgeL(x))), ncol = 2,
                byrow = TRUE))
}

## FIXME: a plot restriction table
## which works by converting restriction to poset


## to plot a null poset (no restrictions)
## just pass the number of drivers as a scalar
plot.poset <- function(x, names = NULL, addroot = FALSE,
                       box = FALSE, ...) {
  if(is.null(names)) {
    if(addroot) names <- c("Root", 1:max(x))
    else names <- 1:max(x)
  }
  plot(poset.to.graph(x, names, addroot), ...)
  if(box)
    box()
}



plot2.oncotree <- function(x, weights = "estimated",
                           edgeFontSize = 8, roundEdgeWeight = 2,
                           addToLabel = 15, ## how many tabs to add, to
                                           ## separate the labels from edge
                           ## angleLabel = NULL,
                           box = TRUE,
                           nodeFontSize = 10,
                         ...) {
    ## Produce nicer plots than available with default of oncotree
    ## but less cumbersome than going through the pst-tree route

    ## sometimes no est.weight, if first fit failed
    if( (weights == "estimated") && (is.null(x$parent$est.weight))) {
        warning("Setting weights in the plot to observed, as no estimated component")
        weights <- "observed"
    }
        
    wg <- switch(weights,
                 "estimated" = x$parent[["est.weight"]][-1],
                 "observed" = x$parent[["obs.weight"]][-1],
                  stop("unknown value for weights")
                 )
    
    parents <- x$parent$parent[-1]
    children <- x$parent$child[-1]
    edgeNamesForgraphNEL <- paste0(parents, "~", children)
    weights <- round(wg, roundEdgeWeight)
    ## tabs will not work when creating pdfs or eps
    addEmpty <- paste(rep(" ", addToLabel), collapse="")
    weights <- paste0(as.character(weights), addEmpty)
    names(weights) <- edgeNamesForgraphNEL
    eattrs <- list(label = weights)
    ## nope, this is for node labels
    ## if(!is.null(angleLabel)) {
    ##     labelangle <- rep(angleLabel, length(weights))
    ##     names(labelangle) <- edgeNamesForgraphNEL
    ##     eattrs <- list(label = weights, labelangle = labelangle)
    ## }
    ## gdf is igraph object
    gdf <- graph.data.frame(data.frame(parents = x$parent$parent[-1],
                                       children = x$parent$child[-1]),
                            ## if weight is added, the edge curves
                            ## which is very ugly. Add via edgeAttrs
                            ## weight = wg),
                            directed = TRUE,
                            vertices = NULL)
    gfn <- igraph.to.graphNEL(gdf)
    plot(gfn, edgeAttrs = eattrs,
         attrs = list(edge = list(fontsize = edgeFontSize),
             node = list(fontsize = nodeFontSize)),
         ...)
    if(box)
        box()
}

read.diprog <- function(dirname, nodenames){
    ## similar to what happens with posets, there is no explicit zero.
    f <- dir(path = dirname, pattern = "edges", full.names = TRUE)
    if(length(f) != 1)
        stop("The number of edges files is different from 1")
    if(file.info(f)$size > 0) {
        m1 <- as.matrix(read.csv(f,
                                 header = FALSE, stringsAsFactors = FALSE))
        ## closely based on read.poset
        colnames(m1) <- rownames(m1) <- NULL
        missing.nodes <- setdiff(nodenames, unique(as.vector(m1)))
        if(length(missing.nodes) == length(nodenames))
            warning("Did you use the same names? length(missing.nodes) == length(nodenames)")
    } else {
        missing.nodes <- nodenames
        m1 <- NULL
    }
    if(length(missing.nodes)) {
        mnm <- matrix(c(missing.nodes, rep(NA, length(missing.nodes))),
                      ncol = 2)
        m1 <- rbind(m1, mnm)
    }
    return(m1) ## this is a like the output from read.poset
}

read.bic.diprog <- function(dirname) {
    f <- dir(path = dirname, pattern = "summary.csv", full.names = TRUE)
    return(as.numeric(scan(f, skip = 1, sep = ",", what = "", quiet = TRUE)[4]))
}

example.create.data <- function(n = 100, p = 7) {
    x <- matrix( sample(c(1, 0), n * p, replace = TRUE), ncol = p )
    x[x[, 1] == 0, 2] <- 0
    x[x[, 1] == 1, 2] <- 1
    x[1:10, 2] <- sample(c(1, 0), 10, replace = TRUE)
    colnames(x) <- 1:p
    x
}

example.create.dep.data <- function(n = 1000, p = 5) {
## a simple one with dependency
    x1 <- sample(c(1, 0), n, replace = TRUE)
    x2 <- sample(c(1, 0), n, replace = TRUE)
    x3 <- rep(0, n)
    x3[(x1 == 1) & (x2 == 1)] <- 1
    x <- cbind(x1, x2, x3)
    x <- cbind(x,
               matrix( sample(c(1, 0), n * (p - 3), replace = TRUE),
                      ncol = (p - 3 )))
    
    colnames(x) <- 1:ncol(x)
    return(x)
}



run.one.diprog <- function(data, n = "MPN", e = 0.2, p = 3, t = 20, m = 1000,
                       callDiProg = "python ~/Sources/diprog/DiProg.py",
                       addname = NULL, ## I will definitely use this
                       dirname = NULL,
                       rmfile = TRUE,
                       silent = TRUE,
                       type.out = "graphNEL") {
    ddp <- data
    rownames(ddp) <- NULL
    ## colnames(ddp) <- 1:ncol(ddp) ## so we can use "poset.to.graph"
    ## we want it to be unique
    if(is.null(dirname)) {
        dirname <- tempfile()
        dirname0 <- NULL
        if(!is.null(addname)) {
            dirname0 <- dirname
            ## FIXME: do this in OS indep way
            dirname <- paste0(dirname, "/", addname, "p", p)
        }
        dir.create(dirname, recursive = TRUE)
        if(!silent)
            cat(paste("\n Created dir", dirname))
        
    }
    ## FIXME: do this in OS indep way
    fname <- paste0(dirname, "/input.csv")
    write.table(ddp, sep = ",", file = fname,
                quote = FALSE, row.names = FALSE, col.names = TRUE)

    thecall <- paste(callDiProg,
                 paste(c("-n", "-e", "-p", "-t", "-m"),
                       c(n, e, p, t, m), collapse = " "),
                 paste("-d", fname),
                 paste("-o", dirname)
                 )
    zzz <- system(thecall, ignore.stdout = silent)
    if(!silent) cat("\n\n")
    
    poset <- read.diprog(dirname, colnames(ddp))
    bic <- read.bic.diprog(dirname)
    if(!silent)
        cat("\n  BIC = ", bic, "\n")
    ## we want to use poset.to.graph, below, which expects integers But it
    ## is better, since DiProg deals with them, to leave full names, in
    ## case we want to check.
    nn <- 1:ncol(ddp)
    names(nn) <- colnames(ddp)
    poset.int <- matrix(nn[poset], ncol = ncol(poset))
    gr <- poset.to.graph(poset.int, names = c("Root", colnames(data)),
                         addroot = TRUE, type = type.out)
    if(rmfile) {
        files <- dir(dirname, full.names = TRUE)
        sapply(files, file.remove)
        file.remove(dirname)
        if(!is.null(dirname0))
            file.remove(dirname0)
    }
    return(list(graph = gr, BIC = bic, p = p))
}


run.diprog <- function(x, n = "MPN", e = 0.2, p.range = 1:4,
                       t = 200, m = 1000,
                       callDiProg = "python ~/Sources/diprog/DiProg.py",
                       addname = NULL, ## I will definitely use this
                       dirname = NULL,
                       rmfile = TRUE,
                       silent = FALSE,
                       type.out = "graphNEL",
                       cores = 1,
                       returnAll = TRUE) {

    outs <- list()
    outs <- mclapply(p.range,
                    function(z) {
                        run.one.diprog(data = x,
                                       n = n,
                                       e = e,
                                       p = z,
                                       t = t,
                                       m = m,
                                       callDiProg = callDiProg,
                                       addname = addname,
                                       dirname = dirname,
                                       rmfile = rmfile,
                                       silent = silent,
                                       type.out = type.out)
                    },
                    mc.cores = cores
                     )
    ## yes, it is the max, and we want the smallest network with the best
    best <- which.max(lapply(outs, function(x) x$BIC))
    if(!silent)
        cat("\n Best solution with BIC ", outs[[best]]$BIC,
            " at p = ", p.range[best], "\n")

    if(returnAll)
        return(c(BestSolution = outs[best], OtherSolutions = outs[-best])) ## first element is the best
    else
        return(BestSolution = outs[[best]])
}



call.external.cbn <- function(data, file = "testcbn", eparam = 0.05,
                              temp = 1, steps = 200, silent = FALSE,
                              init.poset = TRUE, cores = NULL) {
  ## I assume h-cbn and ct-cbn are available

  data2 <- cbind(1, data)
  write(c(nrow(data2), ncol(data2)),
        file = paste(file, ".pat", sep = ""),
        sep = " ")
  write(t(data2), file = paste(file, ".pat", sep = ""),
        ncolumns = ncol(data2),
        append = TRUE, sep = " ")
  write(c("null", colnames(data)),
        file = paste(file, ".prf", sep = ""),
        sep = " ")
  if(is.null(cores)) {
    OMPthreads <- detectCores()
  } else{
    OMPthreads <- cores
  }
  ompt <- paste("export OMP_NUM_THREADS=", OMPthreads, "; ", sep = "")
  
  dir.create(file)
  
  if(init.poset) {
    write.linear.poset(data, file)
  } else { ## Use ct-cbn to search and create starting poset;
    ## possibly eternal. NOT RECOMMENDED
      warning("Not using an intial poset can take VERY long")
      writeLines(as.character(c(ncol(data), 0)),
                 con = paste(file, ".poset", sep = ""))
      ## First create the lambda file
      zzz <- system(paste(ompt , paste("h-cbn -f",  file, "-w")),
                    ignore.stdout = silent)
      if(!silent) cat("\n\n")
      ## this call requires a lambda file
      zzz <- system(paste(ompt, paste("h-cbn -f",  file, "-e", eparam,
                                      "-w -m")), ignore.stdout = silent)
      if(!silent) cat("\n\n")
  }
  ## Remove option -m, the printing of most likely path as
  ##    - we do not use it now
  ##    - it can lead to strange problems getting millions of ceros printed out
  ## zzz <- system(paste(ompt, paste("h-cbn -f",  file, "-s", 
  ##                                 "-T", temp,  "-N", steps,
  ##                                 "-m -w")), ignore.stdout = silent)
  
  zzz <- system(paste(ompt, paste("h-cbn -f",  file, "-s", 
                                  "-T", temp,  "-N", steps,
                                  "-w")), ignore.stdout = silent)
  if(!silent) cat("\n\n")
  ## the final poset in file/00000.poset
}

read.poset <- function(dirname, maxn, verbose = FALSE) {
    ## Read a poset as generated by h-cbn
    tmp <- scan(paste(dirname, "/00000.poset", sep = ""),
                quiet = !verbose)
    tmp <- tmp[-(length(tmp))]
    nn <- tmp[1]
    tmp <- tmp[-1]
    tmp <- matrix(tmp, ncol = 2, byrow = TRUE)
    ## nodes with no ancestors or descendants
    missing.nodes <- setdiff(1:nn, unique(tmp))
    if(length(missing.nodes)) {
        if(verbose)
            message("Reading a poset with missing nodes") ## this is OK if a node not placed in the graph
        if(maxn != nn)
            warning("maxn != nn and missing nodes. Probably should not happen") ## FIXME: stop here?
        mnm <- matrix(c(missing.nodes, rep(NA, length(missing.nodes))),
                      ncol = 2) ## this is perfectly OK, and
                                ## poset.to.graph will do what it should
        tmp <- rbind(tmp, mnm)  ## ditto.
    } else {
        if(maxn != nn)  ## FIXME: I think I should stop in this case,
                        ## unless the user says o.w. But poset.to.graph will break anyway.
            warning("No missing nodes but maxn != nn. Probably should not happen")
    }
    return(tmp)
}



## Later, we will want much more info recovered, such as the probs, etc.
## But unclear where from. The .lambda file is written at the start, so it
## seems not to get updated at end, in contrast to the poset.

## I think I'd need to call external a second time (and a different
## external), a second time, once we are done estimating the poset. But
## this is all completely unclear.

run.cbn <- function(x, file,
                    temp = 1, steps = NULL,
                    silent = FALSE,
                    type.out = "graphNEL",
                    init.poset = TRUE,
                    eparam = 0.05,
                    rmfile = TRUE,
                    cores = NULL) {
    ## FIXME: allow for file to be NULL, as in run.diprog
    if(is.null(steps))
        steps <- ncol(x)^2 ## Their default

    ## their defaults are temperature = 1,
    ## and number of steps = number of genes ^ 2
    zzz <- call.external.cbn(x, file = file, eparam = eparam,
                             temp = temp, steps = steps,
                             silent = silent, init.poset = init.poset,
                             cores = cores)
    cnames <- colnames(x)
    poset <- read.poset(file, ncol(x))
    ## Actually, leave the root in there. But when getting the transitive
    ## closure, do not use.
    ## Why leave root? Because oncotree and treemix seem to use it.
    
    gr <- poset.to.graph(poset, names = c("Root", cnames),
                         addroot = TRUE, type = type.out)
    ## gr <- poset.to.graph(poset, names = cnames,
    ##                      addroot = FALSE, type = type.out)

    if(rmfile) {
        files.created <- paste(file, c(".pat", ".prf", ".log", ".poset",
                                       ".lambda"), sep = "")
        file.remove(files.created)
        file.remove(paste(file, "/00000.poset", sep = ""))
        file.remove(file)
    }
    return(gr)
}

## simple example
## x <- matrix(sample(c(0, 1), 500, replace = TRUE), ncol = 5)
## x[x[, 1] == 0, 2] <- 0
## x[x[, 1] == 1, 2] <- 1
## x[x[, 1] == 1, 4] <- 1
## x[x[, 1] == 0, 4] <- 0
## x[1:3, 4] <- c(1, 1, 1)
## x[1:3, 2] <- c(0, 1, 0)
## colnames(x) <- letters[4:8]

## oo <- run.cbn(x, file = "~/tmp/ff213", rnfile = FALSE)
## missing intermediate and final nodes
## look at this too: read.poset(dirname = "~/tmp/ff213", 5)

## Yes, CBN and DiProg, return all nodes, even if some has freq 0.
## z <- x
## z[, 3] <- 0

## run with oncotree, CBN, and DiProg.
## DiProg places them in the DAG file as nodes.

## CBN gives the total number as the number of nodes in poset, even if you
## pass events with no occurrences. And the .lambda file contains
## estimates of for all.

run.oncotree <- function(x, type.out = "graphNEL",
                         error.fun = function(x, y) { sum((x - y)^2)},
                         hack.all.occurrences = FALSE) {
    if(hack.all.occurrences) {
        cs <- colSums(x)
        nsubs <- nrow(x)
        all.occurr <- which(cs == nsubs)
        if(length(all.occurr)) {
            message(" Using the hack for all occurrences")
            rows.flip <- sample(seq.int(nsubs), length(all.occurr))
            mm <- cbind(rows.flip, all.occurr)
            x[mm] <- 0
        }
    }
    ## yes, ugly, but I do not want to specify it
    onco.fit <- oncotree.fit(x, error.fun = error.fun)
    
    gdf <- graph.data.frame(data.frame(parents = onco.fit$parent$parent[-1],
                                       children = onco.fit$parent$child[-1]
                                       ),
                            directed = TRUE,
                            vertices = NULL)
    if(type.out == "adjmat") return(get.adjacency(gdf, sparse = FALSE))
    else if (type.out == "graphNEL") return(igraph.to.graphNEL(gdf))
}

## require(Rtreemix)
## to get the oncogenetic tree call with K = 1 and noise = FALSE
run.rtreemix <- function(x, K = 3, noise = TRUE, only.graphnel = TRUE,
                         equal.edgeweights = TRUE) {

    rtm <- new("RtreemixData", Sample = x,
             Events = c("Root", colnames(x)))
    ot <- fit(data = rtm, K = K, noise = noise,
              equal.edgeweights = equal.edgeweights)

  if(only.graphnel) {
    if(K == 1) return(getTree(ot, k = 1))
    ltrees <- list()
    for(i in 1:K)
      ltrees[i] <- getTree(ot, k = i)
    return(ltrees)
  } else {
    return(ot)
  }
}

ssugar.true.graph <- function(poset) {
    ## return the true graph from a poset, with labeled nodes and a Root
    ## This is input to metrics1 This is just to avoid repeated calls.

    ## This function was first used and defined in directory-analysis-3.R
    ## and directory-analysis-2.R

    return(poset.to.graph(poset,
                          names = c("Root", paste("G.", 1:max(poset), sep = "")),
                          addroot = TRUE))
}


large.adj.mat <- function(x, names = c("Root", paste("G", 1:60, sep = ".")),
                          transclos = FALSE) {

    ## Return an adjacency matrix that contains up to all the nodes named
    ## in "names". If they are not there, rows and columns are added to
    ## the matrix, and filled with zeroes. Yes, zeroes. If a node does not
    ## exist originally we DO NOT add an edge from root to that node as we
    ## should not.
    
    ## x is the graphNEL object, and has root
    
    LA <- matrix(0L, nrow = length(names), ncol = length(names))
    rownames(LA) <- colnames(LA) <- names

    if(!inherits(x, "graphNEL"))
        return(LA)
    if(!transclos)
        A <- get.adjacency(igraph.from.graphNEL(x),
                           sparse = FALSE)
    else
        A <- get.adjacency(igraph.from.graphNEL(transitive.closure(x)),
                           sparse = FALSE)

    ## NO, we do NOT remove the root, as we want the adjacency matrix.
    ## Without the root, the adjacency matrix for p903 would
    ## be the same with and without node 4, as no one depends on it
    ## and it depends on none, except Root


    ## Keeping the root for adjacency matrices seems to be what Yan et
    ## al. 2006, SAMB, do
    
    ## proot <- which(rownames(A) == "Root")
    ## if(proot !=  which(colnames(A) == "Root"))
    ##     stop("Root should be in same col and pos")
    ## A <- A[-proot, -proot]
    
    w1 <- which(A == 1, arr.ind = TRUE)
    rr <- rownames(A)[w1[, 1]]
    cc <- colnames(A)[w1[, 2]]
    ## make sure names is what we think
    if(! all( rownames(A) %in% names ))
        stop("eh? not all the names in vector names")
    mi <- cbind(rr, cc)
    LA[mi] <- 1L
    return(LA)
}




DiffAdjMat <- function(t1, t2, transclos = FALSE) {
    ## computes the difference between adjacency matrices of t1 and t2,
    ## where t1 and t2 are two graphNEL objects. We expect one node to be
    ## called "Root"

    ## Lots of code below, but for now only use the simple
    ## difference. That is, by the way, the same as the edit distance of
    ## Hainke et al., 2012, when the two trees have the same number
    ## (should also be the same identity, but they don't check) of
    ## nodes. We do not require that. They could have different
    ## nodes. That is handled by large.adj.mat.

    ## large.adj.mat adds time, but ensures same nodes AND adjacency
    ## matrices in same order.

    ## the tc is for transitive closure
    ## nI <- norm(m1 - m2, type = "I") ## the statistic in Yin et al. 2006, SAGMB
    ## nS <- sum(abs(m1 - m2))
    ## nI.tc <- norm(m1.tc - m2.tc, type = "I") ## the statistic in Yin et al. 2006, SAGMB
    ## nS.tc <- sum(abs(m1.tc - m2.tc))

    ## return(data.frame(name, nI, nI/df.all$ndr[i], nS, nS/df.all$ndr[i],
    ##                   nI.tc, nI.tc/df.all$ndr[i], nS.tc, nS.tc/df.all$ndr[i],
    ##                   stringsAsFactors = FALSE))

    ## for aesthetic purposes, we want Root first, and then the rest, sorted.
    ## Yes, Root is already present, but this way I force it to be first.
    getnodes <- function(x) {
        if(inherits(x, "graphNEL"))
            return(nodes(x))
        else
            return(NULL)
    }
    nodes.t1 <- getnodes(t1)
    nodes.t2 <- getnodes(t2)
    allN <- union("Root", sort(union(nodes.t1, nodes.t2)))
    m1 <- large.adj.mat(t1, names = allN, transclos = transclos)
    m2 <- large.adj.mat(t2, names = allN, transclos = transclos)
    nS <- sum(abs(m1 - m2))
    return(nS)
}



FP.TP.FN.counts <- function(t1, t2) {
    ## FP, TP, FN. The TN are a different story.

    ## t1 and t2 are graphNEL objects. And we return the coutns using the
    ## transistive closure of the relations.

    if(!inherits(t1, "graphNEL")) {
        ex <- integer(0)
    } else {
        ex <- veltcs(t1, remove.root = TRUE)
    }

    if(!inherits(t2, "graphNEL")) {
        ey <- integer(0)
    } else {
        ey <- veltcs(t2, remove.root = TRUE)
    }

    FP <- length(setdiff(ex, ey))
    FN <- length(setdiff(ey, ex))
    TP <- length(intersect(ey, ex))
    return(c(FP = FP, FN = FN, TP = TP))
}


performance.stats <- function(t1, t2) {
    return(c(Diff = DiffAdjMat(t1, t2), FP.TP.FN.counts(t1, t2)))
}



## how to get depth of node from a poset
## To use metrics for order.
## See Youn and Simon for metrics; something like MSE?


order.muts <- function(x){
    ## This we use for oncog. trees, as nothing better.
    ## And we make the task particularly easy.

    nn <- nodes(x)
    rn <- which(nn == "Root")
    if(!length(rn))
        stop("There must be a node called Root")
    nn <- nn[-rn]
    ll <- sapply(get.shortest.paths(igraph.from.graphNEL(x),
                             from = "Root", to = nn,
                             output = "epath"),
                 length)
    names(ll) <- nn
    class(ll) <- "orderMuts"
    return(ll)
}
    
compare.order.muts <- function(x, reference) {
    ## FIXME:
    ##  - we do not take into account how many are not recovered
    ##    or how many are falsely recovered. But this is a filtering thing.
    ##  - what should the statistic be? MSE per true order? Or average MSE?
    nc <- intersect(names(x), names(reference))
    return(cbind( x = x[nc], reference = reference[nc]))
}
    

## d1 <- matrix(rbinom(1000 * 14, 1, 0.1), ncol = 14)
## rownames(d1) <- paste("ID", 1:1000, sep = "")
## colnames(d1) <- paste(1:14, sep = "")

## d1[, 1] <- rbinom(1000, 1, 0.8)
## d1[ d1[, 1] == 1, 2] <- rbinom(sum(d1[, 1] == 1) , 1, 0.99)
## d1[ d1[, 1] == 0, 2] <- rbinom(sum(d1[, 1] == 0) , 1, 0.01)

## d1[ d1[, 2] == 1, 3] <- rbinom(sum(d1[, 2] == 1) , 1, 0.99)
## d1[ d1[, 2] == 0, 3] <- rbinom(sum(d1[, 2] == 0) , 1, 0.01)

## tmp <- bcbn(d1, samples = 10000)



bcbn <- function(data, p.thresh = 0.6,
                 cores = 4, chains = 4, samples = 25000,
                 thin = 10, epsilon = 0.05,
                 poset.mode = TRUE) {

    ## based on the dosim function in the bcbn package from T. Sakoparnig
    
    ## My additions:

    ## p.thresh is the threshold so that we return the transitiveClosure
    ## when the posterior for a direct relationship is larger than the threshold.

    ## I return, or not (poset.mode) the poset mode
    
    registerDoMC(cores=cores)
    n_samples <- samples
    n_chains <- chains
    thin <- thin
    epsilon <- epsilon
    theta=0
    
    n<-dim(data)[2]
    n_cases <- dim(data)[1]

    mlist<-list()
    edgelist<-list()
    l=0
    converged = 0
    repeat {
        l=l+1
        rets <- foreach( i = 1:n_chains ) %dopar% {
            print(paste("chain:",i))
            if( theta == 0 ) {
                theta=as.double(runif(n))
            }
            
            if(length(mlist)!=0) {
                edges_in = c(t(edgelist[[i]][n_samples][[1]]))
                theta = as.double(mlist[[i]][n_samples,1:n])
                epsilon = mlist[[i]][n_samples,n+1]
            } else {
                edges_in = as.integer(rep(0,n*n))
            }
            
            ret<-.C("sample_full_cbn", theta, as.integer(n), as.double(epsilon), edges_in, as.integer(n_samples), as.integer(thin), as.integer(c(t(data))), as.integer(n_cases), theta_out=as.double(rep(0,n*n_samples)), epsilon_out=as.double(rep(0,n_samples)), edges_out=as.integer(rep(0,n_samples*n*n)), log_posterior_out=as.double(rep(0,n_samples)))
            
        } ## here the chains are done
        
        
        mlist<-list()
        mcmclist<-list()
        edgelist<-list()
        
        ## X11()
        ## par(mfrow=c(2,2))
        
        for( i in 1:n_chains) {
            theta_m<-matrix(rets[[i]]$theta_out,ncol=n,byrow=T)
            paramatrix<-cbind(theta_m,rets[[i]]$epsilon_out,rets[[i]]$log_posterior_out)
            mlist[[i]]<-paramatrix
            mcmclist[[i]]<-mcmc(paramatrix)
            sublist<-list()
            edgesum<-0
            for(k in 1:n_samples) {
                sublist[[k]]<-matrix(rets[[i]]$edges_out[((k-1)*n*n+1):(k*n*n)],ncol=n,byrow=T)
                edgesum<-edgesum+matrix(rets[[i]]$edges_out[((k-1)*n*n+1):(k*n*n)],ncol=n,byrow=T)
            }
            edgelist[[i]]<-sublist
                                        #print(diag(var(paramatrix)))
##        print(summary(paramatrix))
        
        ## image(edgesum)
        }
                                        #image(edgelist[[1]][[getmode(edgelist[[1]])]],main=paste("mode of chain 1 run ",l))
        mclist<-mcmc.list(mcmclist)
        
        terr<-try(gdiag<-gelman.diag(mclist))
        if (class(terr) != 'try-error') {
            if ( max( gdiag$psrf[,1] ) < 1.1 ) {
                if ( converged==1 ) {
                    break
                }
                converged=1
            }
            
            ## print(gdiag)
        }
        print(paste("finished run:",l))
    }
    ## this is expecting, I think, 4 chains. It is hard coded!!!
    
    elist<-c(edgelist[[1]],edgelist[[2]],edgelist[[3]],edgelist[[4]])
    ## X11()
    ## plotTop(elist,4) ## seems eternal
                                        #
    ## print(gdiag)
                                        #jpeg("400_0.1_dens.jpg",quality=100,height=600,width=600)
    
    ## par(mfrow=c(3,4))
    ## plot(mclist,trace=F,ask=T,auto.layout=F)
    
                                        #dev.off()
                                        #jpeg("400_0.1_structures.jpg",quality=100,height=200,width=1000)
                                        #plotTop(elist,4)
                                        #dev.off()
                                        #
    marginaledges<-marginal(elist)


    if(poset.mode)
        posetMode <- elist[getmode(elist)]
    else
        posetMode <- NA

    p.transitiveClosure <- transitiveClosure(marginaledges  > p.thresh)

    return(list(
        marginaledges = marginaledges,
        p.transitiveClosure = p.transitiveClosure,
        posetMode = posetMode))
}



## simulate from a poset. Like in non-timed OTs
## FIXME This is a 2 second piece of code. It is not correct!!

simposet <- function(poset, p) {

    ## poset should include the root, coded as 1, so all genes are
    ## displaced. And should be a complete poset, showing all deps,
    ## including from Root.
    
    num.genes <- max(poset) - 1 ## as root is not a gene
    genotype <- c(1, rep(NA, num.genes))

    poset <- data.frame(poset)
    poset$runif <- runif(nrow(poset))
    
    
    poset ## has several columns: first is parent, 2 is child, then we have
    ## the randon unifor number.
    
    ## p es la probabilidad que fijamos


    ## this.relation.prob.OK could be done outside, but having it inside
    ## the loop would allow to use different thresholds for different
    ## relationships
    
    for (i in (1:nrow(poset))) {
        child <- poset[i, 2]
        this.relation.prob.OK <- as.numeric(poset[i, "runif"] > p)
        the.parent <- genotype[ poset[i, 1] ]
        genotype[child] <- this.relation.prob.OK * the.parent
    }
    return(genotype)
}


