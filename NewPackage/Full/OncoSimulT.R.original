## Testing CBN programs available
.._OncoSimul_test.ctcbn <- system("ct-cbn -h", ignore.stdout = TRUE)
.._OncoSimul_test.hcbn <- system("h-cbn -h", ignore.stdout = TRUE)
if(.._OncoSimul_test.ctcbn || .._OncoSimul_test.hcbn) {
    warning(paste(
        "\n\n",
        "\n******************************************************",
        "\n******************************************************\n",
        "          OncoSimulT installation warning:\n",
        "\n",
        "The external programs h-cbn and ct-cbn were not found.",
        "You will not be able to use CBN.",
        "You can download them from http://www.bsse.ethz.ch/cbg/software/ct-cbn",
        "You might want to use the makefile Makefile-cbn-modified-RDU under",
        "the miscell-code directory of this package.\n",
        "\n******************************************************",
        "\n******************************************************",
        "\n\n"
        )
            )
}
rm(.._OncoSimul_test.ctcbn)
rm(.._OncoSimul_test.hcbn)
## Done with the testing

## FIXME add similar testing for DiProg



## FIXME: change typeFitness to "modelType"

## edgeList is no longer available. Now it is called with the
## (of course, much easier to decode) edgeL. And not even an alias
## to former name. Oh well.


## Store in object features of the call.
## Done from C++?


## wrapTi <- function(v1_, v2_, v3_, v4_, v5_,  v6_) {
##   .Call("wrap_ti", v1_, v2_, v3_, v4_, v5_, v6_, PACKAGE = "MatherR")
## }

## wrapAlgo2 <- function(v1_, v2_, v3_, v4_, v5_, v6_, vd_, v7_) {
##   .Call("wrap_Algo2", v1_, v2_, v3_, v4_, v5_, v6_,  vd_, v7_, PACKAGE = "MatherR")
## }

## wrapAlgo3 <- function(v1_, v2_, v3_, v4_, v5_,  v6_, vd_, v7_) {
##   .Call("wrap_Algo3", v1_, v2_, v3_, v4_, v5_, v6_, vd_, v7_, PACKAGE = "MatherR")
## }

## wrapFitnessLinearRcpp <- function(allGenotypes_, genNum_, birthRate_,
##                               s_, numDrivers_,  v6_) {
##   .Call("wrap_fitness_linear_Rcpp", allGenotypes_, genNum_, birthRate_,
##         s_, numDrivers_,  v6_, PACKAGE = "MatherR")
## }

## wrapFitnessCBNRcpp <- function(mutatedPos_, genotypes_,
##                            genNum_, restrictTable_,
##                            numDrivers_,
##                            birthRate_, s_,
##                            fitnessParent_, typeCBN_, retval_){
  
##   .Call("wrap_fitness_CBN_Rcpp",mutatedPos_, genotypes_,
##                            genNum_, restrictTable_,
##                            numDrivers_,
##                            birthRate_, s_,
##                            fitnessParent_, typeCBN_, retval_,
##         PACKAGE = "MatherR")
## }

## wrapFitnessCBNArma <- function(mutatedPos_, genotypes_,
##                            genNum_, restrictTable_,
##                            numDrivers_,
##                            birthRate_, s_,
##                            fitnessParent_, typeCBN_, retval_){
  
##   .Call("wrap_fitness_CBN_Arma",mutatedPos_, genotypes_,
##                            genNum_, restrictTable_,
##                            numDrivers_,
##                            birthRate_, s_,
##                            fitnessParent_, typeCBN_, retval_,
##         PACKAGE = "MatherR")
## }

## wrapFitnessCBNstd <- function(mutatedPos_, genotypes_,
##                               genNum_, restrictTable_,
##                               numDrivers_,
##                               birthRate_, s_,
##                               fitnessParent_, typeCBN_, retval_){
  
##   .Call("wrap_fitness_CBN_std",mutatedPos_, genotypes_,
##         genNum_, restrictTable_,
##         numDrivers_,
##         birthRate_, s_,
##         fitnessParent_, typeCBN_, retval_,
##         PACKAGE = "MatherR")
## }



## Algo5OP <- function(rt) {
## seed_gsl
## crt? or converted outside?
## numDrivers: ncol of crt
## }


## FIXME: how large should initSize_species be???
Algo5 <- function(restrict.table,
                  numGenes,
                  typeFitness,
                  typeCBN,
                  birth, 
                  s,
                  sh,
                  death,
                  mu,
                  initSize,
                  sampleEvery,
                  detectionSize,
                  mutatorGenotype,
                  finalTime,
                  initSize_species = 2000,
                  initSize_iter = 500,
                  seed_gsl = NULL,
                  verbosity = 1,
                  initMutant = -1,
                  ## next two set so no forced sampling by default
                  speciesFS = 40000,
                  ratioForce = 2,
                  max.memory = 20000,
                  max.wall.time = 3600,
                  keep.every = 20,
                  alpha = 0.0015,
                  K = 1000,
                  endTimeEvery = NULL,
                  finalDrivers = 1000,
                  silent = TRUE) {
  ## the value of 20000, in megabytes, for max.memory sets a limit of ~ 20 GB
  
  ## FIXME: check argument types for typeFitness and typeCBN

  if(initSize_species < 10) {
    warning("initSize_species too small?")
  }
  if(initSize_iter < 100) {
    warning("initSize_iter too small?")
  }
  if(keep.every < sampleEvery)
    warning("setting keep.every to sampleEvery")
  if(is.null(seed_gsl)) {## passing a null creates a random seed
    seed_gsl <- as.integer(round(runif(1, min = 0, max = 2^16)))
    if(!silent)
      cat(paste("\n Using ", seed_gsl, " as seed for GSL\n"))
  }

  numDrivers <- nrow(restrict.table)
  if(length(unique(restrict.table[, 1])) != numDrivers)
    stop("EH??!! length(unique(restrict.table[, 1])) != numDrivers)")
  ddr <- restrict.table[, 1]
  if(any(diff(ddr) != 1))
    stop(" any(diff(ddr) != 1")
  ## sanity checks
  if(max(restrict.table[, 1]) != numDrivers)
    stop("max(restrict.table[, 1]) != numDrivers")
  if(numDrivers > numGenes)
    stop("numDrivers > numGenes")
  
  ## this is biologically feasible, but silly and a mess for the C++ code.
  if(initMutant > numDrivers)
    stop("initMutant > numDrivers")

  non.dep.drivers <- restrict.table[which(restrict.table[, 2] == 0), 1]

  ## this makes no sense, and since we can make init size > 1, with a non-allowed
  ## mutant, it can lead to weird results if used carelessly
  if(initMutant > max(non.dep.drivers))
    warning("initMutant > max(non.dep.drivers)!!! This can lead to meaningless results.")
  
  if(initMutant == 0) {
    initMutant <- sample(non.dep.drivers, 1) - 1
  }


  if( (typeFitness == "bozic1") && (mutatorGenotype) )
    warning("Using fitness bozic1 with mutatorGenotype;",
            "this will have no effect.")

  if( (typeFitness == "exp") && (death != 1) )
    warning("Using fitness exp with death != 1")


  if( (is.null(endTimeEvery) || (endTimeEvery > 0)) &&
      (typeFitness %in% c("bozic1", "bozic2", "exp", "log", "linear") )) {
          warning(paste("endTimeEvery will take a positive value. ",
                        "This will make simulations not stop until the next ",
                        "endTimeEvery has been reached. Thus, in simulations ",
                        "with very fast growth, simulations can take a long ",
                        "time to finish, or can hit the wall time limit. "))
      }
  if(is.null(endTimeEvery))
    endTimeEvery <- keep.every
  if( (endTimeEvery > 0) && (endTimeEvery %% keep.every) )
    warning("!(endTimeEvery %% keep.every)")
  ## a sanity check in restricTable, so no neg. indices for the positive deps
  neg.deps <- function(x) {
    ## checks a row of restrict.table
    numdeps <- x[2]
    if(numdeps) {
      deps <- x[3:(3+numdeps - 1)]
      return(any(deps < 0))
    } else FALSE
  }
  if(any(apply(restrict.table, 1, neg.deps)))
    stop("Negative dependencies in restriction table")

  ## transpose the table
  rtC <- convertRestrictTable(restrict.table)

     
  ## return the matching call? call <- match.call()
  ## and then return(c(.Call(), call))
  call <- match.call()
  return(c(.Call("Algorithm5P",
                 rtC,
                 numDrivers,
                 numGenes,
                 typeCBN,
                 birth, 
                 s, 
                 death,
                 mu,
                 initSize,
                 sampleEvery,
                 detectionSize,
                 finalTime,
                 initSize_species,
                 initSize_iter,
                 seed_gsl,
                 verbosity,
                 speciesFS,
                 ratioForce,
                 typeFitness,
                 max.memory,
                 mutatorGenotype,
                 initMutant,
                 max.wall.time,
                 keep.every,
                 alpha,
                 sh,
                 K,
                 endTimeEvery,
                 finalDrivers,
               PACKAGE = "OncoSimulT"),
           call = call,
           NumDrivers = numDrivers,
           initMutant = initMutant))
}


colnames.to.pops.by.time <- function(pops.by.time) {
  if(prod(dim(pops.by.time)) > 1) {  
    ## colnames(pops.by.time) <- rep("", ncol(pops.by.time))
    colnames(pops.by.time) <- c("Time",
                                paste("Sp_", 1:tmp$NumSpecies, sep = ""))
  }
}


create.muts.by.time <- function(tmp) { ## tmp is the output from Algorithm5
  if(tmp$NumSpecies > 1) {
    NumMutations <- apply(tmp$Genotypes, 2, sum)
    muts.by.time <- cbind(tmp$pops.by.time[, c(1), drop = FALSE] ,
                          t(apply(tmp$pops.by.time[, -c(1), drop = FALSE], 1,
                                  function(x) tapply(x, NumMutations, sum))))
    colnames(muts.by.time)[c(1)] <- "Time"
  } else {
    muts.by.time <- tmp$pops.by.time
  }
  return(muts.by.time)
} 
  

create.drivers.by.time <- function(tmp, numDrivers) {
  CountNumDrivers <- apply(tmp$Genotypes[1:numDrivers, ,drop = FALSE], 2, sum)
  if(tmp$NumSpecies > 1) {
    if(length(unique(CountNumDrivers )) > 1) {
      drivers.by.time <- cbind(tmp$pops.by.time[, c(1), drop = FALSE] ,
                               t(apply(tmp$pops.by.time[, -c(1), drop = FALSE], 1,
                                       function(x) tapply(x, CountNumDrivers, sum)))) 
    } else {
      drivers.by.time <- cbind(tmp$pops.by.time[, c(1), drop = FALSE] ,
                               rowSums(tmp$pops.by.time[, -c(1), drop = FALSE]))
    }
    colnames(drivers.by.time) <- c("Time",
                                   paste("dr_", colnames(drivers.by.time)[-c(1)],
                                         sep = ""))
  } else {
    drivers.by.time <- NULL
  }
  return(drivers.by.time)
} 




sampleZZ <- function(zz, seed = "auto"){
  if(seed == "auto") {
    ## paste as numeric(hostname()) y el segundo from time
  }

  fname <- paste(fileroot, randomstgring)
  save(...., compress = FALSE)

  
}


## For plotting, this helps decrease huge file sizes, while still showing
## the start of each clone, if it was originally recorded.

thin.pops <- function(x, keep = 0.1, min.keep = 3) {
    norig <- nrow(x$pops.by.time)
    keep1 <- round(seq.int(from = 1, to = norig,
                           length.out = round(norig * keep)))
    keep2 <- apply(x$pops.by.time[, -1],
                   1, function(x) any((x[x > 0] < min.keep)))
    keep <- sort(union(keep1, keep2))
    x$pops.by.time <- x$pops.by.time[keep, ]
    return(x)
}

plotPopAndDrivers <- function(z, col = c(8, "orange", 6:1),
                              log = "y",
                              ltyPop = 2:6,
                              lwdPop = 0.2,
                              ltyDrivers = 1,
                              lwdDrivers = 3,
                              xlab = "Time units",
                              ylab = "Number of cells",
                              addtot = FALSE,
                              addtotlwd = 0.5,
                              yl = NULL,
                              thinPops = TRUE,
                              thinPops.keep = 0.1,
                              thinPops.min = 2,
                              ...
                              ) {

    ## uses both of plotPop and plotDrivers in a single plot.
    if(thinPops)
        z <- thin.pops(z, keep = thinPops.keep, min.keep = thinPops.min)
    gc()
    ndr <- apply(z$Genotypes[1:z$NumDrivers, , drop = FALSE], 2, sum)

    if(is.null(yl))
        yl <- c(1, max(apply(z$pops.by.time[, -1, drop = FALSE], 1, sum)))

    plotPop(z,
            ndr = ndr, 
            xlab = xlab,
            ylab = ylab,
            lty = ltyPop,
            col = col, 
            ylim = yl,
            lwd = lwdPop,
            axes = FALSE,
            log = log,
            ...)
    par(new = TRUE)
    plotDrivers0(z,
                 timescale = 1,
                 trim.no.drivers = FALSE,
                 xlab = "", ylab = "",
                 lwd = lwdDrivers,
                 lty = ltyDrivers,
                 col = col, 
                 addtot = addtot,
                 addtotlwd = addtotlwd,
                 log = log, ylim = yl,
                 ...)
}


plotPop <- function(z, ndr = NULL, na.subs = TRUE,
                    log = "y", type = "l",
                    lty = 1:8, col = 1:9, ...) {

    ## if given ndr, we order columns based on ndr, so clones with more
    ## drivers are plotted last

    y <- z$pops.by.time[, 2:ncol(z$pops.by.time)]
    
    if(na.subs){
        y[y == 0] <- NA
    }
  if(!is.null(ndr)) {
      ## could be done above, to avoid creating
      ## more copies
      oo <- order(ndr)
      y <- y[, oo]
      ndr <- ndr[oo]
      col <- col[ndr + 1]
  }
  matplot(x = z$pops.by.time[, 1],
          y = y,
          log = log, type = type,
          col = col, lty = lty,
          ...)
}


plotDrivers0 <- function(x,
                         timescale = 4,
                         trim.no.drivers = TRUE,
                         addtot = TRUE,
                         addtotlwd = 2,
                         na.subs = TRUE, log = "y", type = "l",
                         lty = 1:9, col = c(8, "orange", 6:1),
                         lwd = 2,
                         ...) {

  z <- create.drivers.by.time(x, x$NumDrivers)
  if(trim.no.drivers && x$MaxDriversLast) {
      fi <- which(apply(z[, -c(1, 2), drop = FALSE], 1, function(x) sum(x) > 0))[1]
      z <- z[fi:nrow(z), , drop = FALSE]
  }
  y <- z[, 2:ncol(z), drop = FALSE]
  if(na.subs){
      y[y == 0] <- NA
  }
  if(timescale != 1) {
      time <- timescale * z[, 1]
  } else {
      time <- z[, 1]
  }
  if(nrow(y) <= 2) type <- "b"
  matplot(x = time,
          y = y,
          type = type, log = log, lty = lty, col = col, lwd = lwd,
          ...)
  if(addtot) {
      tot <- rowSums(y, na.rm = TRUE)
      lines(time, tot, col = "black", lty = 1, lwd = addtotlwd)
  }
  ## will need to add a legend
  legend(x = "topleft",
         title = "Number of drivers",
         lty = lty, col = col, lwd = lwd,
         legend = (1:ncol(y)) - 1)
}


## I should create a single function that does both or a single one. Lots
## of repetition. Nope, not that much. Negligible. So forgotten for now
## plot2 <- function(z,
##                   col = c(8, "orange", 6:1),
##                   log = "y",
##                   ltyPop = 2:6,
##                   lwdPop = 0.2,
##                   ltyDrivers = 1,
##                   lwdDrivers = 3,
##                   xlab = "Time units",
##                   ylab = "Number of cells",
##                   addtot = FALSE,
##                   addtotlwd = 0.5,
##                   yl = NULL,
##                   usendr = TRUE,
##                   plotClones = TRUE,
##                   plotDrivers = TRUE,
##                   ...
##                   ) {
##     ndr <- apply(z$Genotypes[1:z$NumDrivers, , drop = FALSE], 2, sum)

##     if(is.null(yl))
##         yl <- c(1, max(apply(z$pops.by.time[, -1, drop = FALSE], 1, sum)))

##     x <- z$pops.by.time[, 1]
##     y <- z$pops.by.time[, -1, drop = FALSE]
##     if(na.subs){
##         y[y == 0] <- NA
##     }
##     plot(type = "n", x = c(min(x), max(x)), y = yl,
##          xlab = xlab, ylab = ylab, log = log)
##     if(plotClones) {
##         if(usendr) {
##             oo <- order(ndr)
##             y <- y[, oo]
##             ndr <- ndr[oo]
##             col <- col[ndr + 1]
##         }
##         par(new = TRUE)
##         matplot(x = x,
##                 y = y,
##                 log = log, type = type,
##                 col = col, lty = lty,
##                 axes = FALSE,
##                 ...)
##     }
##     if(plotDrivers)  {
##         par(new = TRUE)
##         z <- create.drivers.by.time(z, z$NumDrivers)
##         ## what was this for? Rarely used now
##         ## if(trim.no.drivers && x$MaxDriversLast) {
##         ##     fi <- which(apply(z[, -c(1, 2), drop = FALSE], 1, function(x) sum(x) > 0))[1]
##         ##     z <- z[fi:nrow(z), , drop = FALSE]
##         ## }
##         if(nrow(y) <= 2) type <- "b"
##         matplot(x = time,
##                 y = y,
##                 type = type, log = log, lty = lty, col = col, lwd = lwd,
##                 ...)


##     }
##     ## axis, etc
## }




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
plotdriversdir <- function(...){
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


rtNoDep <- function(numdrivers) {
  ## create a restriction table with no dependencies
  x <- matrix(nrow = numdrivers, ncol = 3)
  x[, 1] <- 1:numdrivers
  x[, 2] <- 0
  x[, 3] <- -9
  return(x)
}


## in the future: Algo5 takes either
## format. The C format is found
## if min value for first col is 0?
## or force to use either one of two args
## restrict.table or c.restrict.table

convertRestrictTable <- function(x) {
  ## to convert the table to the format for C
  ## as there the mutations are numbered from 0

  ## In R the format for a row is:
  ##  - the mutation,
  ##  - the number of mutations on which it depends
  ##  - the actual mutations on which it depends
  ##  - the rest are "-9"
  
  t.restrictTable <- matrix(as.integer(x),
                            ncol = nrow(x), byrow = TRUE)

  t.restrictTable[-2, ] <- t.restrictTable[-2, ] - 1
  return(t.restrictTable)
}

## rename as "selectJointFreq"
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


## veltcs0 <- function(x) {
##     ## now broken, because of changes to edgeList (now called edgeL),
##     ## and other changes in graphNEL, I think.
##   ## return a vector from edge list of transitive closure
##   if(length(x) > 1) {
##     if(class(x[[1]]) == "graphNEL")
##       tc0 <- lapply(x, transitive.closure)
##     else ## adjacency matrix
##       tc0 <- lapply(x, transClos)
##     el <- unlist(lapply(tc0, edgeL), recursive = FALSE)

##   } else{
##     if(class(x) == "graphNEL")
##       tc0 <- transitive.closure(x)
##     else ## adjacency matrix
##       tc0 <- transClos(x)
##     el <- edgeL(tc0)
##   }
##   ell <- sapply(el, function(x) paste(x[1], x[2], sep = "_"))
##   return(unique(ell)) ## important with mixtures
## }


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


adjmat.to.restrictTable <- function(x) {
  ## we have the zero
  ## x <- x[-1, -1]
  num.deps <- colSums(x)
  max.n.deps <- max(num.deps)
  rt <- matrix(-9, nrow = nrow(x),
               ncol = max.n.deps + 2)
  for(i in 1:ncol(x)) {
    if( num.deps[ i ])
      rt[i , 1:(2 + num.deps[ i ])] <- c(i, num.deps[i ], which(x[, i ] != 0))
    else
      rt[i , 1:2] <- c(i , 0)
  }
  return(rt)
}

poset.to.restrictTable <- function(x) {
  x1 <- poset.to.graph(x, names = 1:max(x), addroot = FALSE, type = "adjmat")
  adjmat.to.restrictTable(x1)
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

## have a graph without a null??


poset.to.graph <- function(x, names,
                           addroot = FALSE,
                           type = "graphNEL") {
    ## Intermediate nodes, if no ancestor or descendant, need not
    ## be in the poset.
    ## Any node with index largest than any node with ancestor or descendant
    ## needs to be in the file.
    ## E.g., rbind(c(1,2), c(4, 5)) will get three as a child-less and parent-less node.
    ## But to place 6 you need c(6, NA)

    
    ## We could use something like in run.oncotree,
    ## as a poset also is a set of parent-child,
    ## and we would then need to add the root connections
    ## as is done below in no.ancestor.

    ## But we do not for now. Note we show lonely nodes, which oncotrees
    ## do not.  wait: when using root, we do not have "lonely nodes"
    ## anymore.  But that is irrelevant for metrics based on transitive
    ## closure. Note for Diff, etc.

    ## FIXME: in fact, this is all OK, but is confussing, because I can
    ## have two kinds of posets: ones that are full, with NAs, etc, if
    ## needed. Others that are not, but are fixed by the code here. But if
    ## using the later, the user needs to make sure that the last node is
    ## in the poset.
    
  m <- length(names) 
  m1 <- matrix(0, nrow = m, ncol = m)
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
        m1[x + 1] <- 1
      else
        m1[x] <- 1
    }
    if((length(names) > 1) & addroot) {
      no.ancestor <- which(apply(m1, 2, function(x) all(x == 0)))
      no.ancestor <- no.ancestor[-1]
      m1[cbind(1, no.ancestor)] <- 1
    } ## o.w. do nothing
  }
  if(type == "adjmat") return(m1)
  else if (type == "graphNEL") return(as(m1, "graphNEL"))
  ## does not show the labels
  ## else if (type == "igraph") return(graph.adjacency(m1))
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



## t is in seconds.
## remmber to run with various p, and use BIC score to choose the best.
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


## FIXME: silen should be silent, not redirecto to intern
## a long time ago, here, for playing, temp=10, but I always used
## temp = 1 (see run.cbn). Set to 1 here for consistency.
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



## run.oncotree <- function(x, type.out = "graphNEL",
##                          lonely = FALSE) {
##     ## way too complicated!! simplify!!
##   m <- ncol(x)
##   onco.fit <- oncotree.fit(x)
##   cnames <- colnames(x)
##   ## root gets number 1
##   ## convert to adjacency matrix
##   ## but remember 1 is actually the root, so we add a 1, and later rename
##   parents <- onco.fit$parent$parent.num[-1]
##   children <- match(onco.fit$parent$child[-1], cnames) + 1
##   p.to.child <- cbind( parents , children )
##   m1 <- matrix(0, nrow = m + 1, ncol = m + 1)
##   m1[p.to.child] <- 1
##   colnames(m1) <- c("Root", cnames)
##   rownames(m1) <- c("Root", cnames)
##   if(!lonely) {
##       ## reproduce behavior of Oncotree, not showing nodes without ancestor
##       ## except for root, of course
##       rows.cols.to.keep <- c(1,
##                              unique(match(onco.fit$parent$child[-1], cnames) + 1))
##       m1 <- m1[rows.cols.to.keep, rows.cols.to.keep]
##   }
##   if(type.out == "adjmat") return(m1)
##   else if (type.out == "graphNEL") return(as(m1, "graphNEL"))
## }



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
## run.oncotree.Rtreemix <- function(x) {
##     ## obtain the oncogenetic tree using Rtreemix.
##     ## requires Rtreemix, of course.
##     dd <- new("RtreemixData", Sample = x, Events = c("Root", colnames(x)))
##     rtm <- fit(data = dd, K = 1, noise = FALSE)
##     return(Trees(rtm)[[1]])
## }


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


run.YounSimon <- function(sample.gene.mutation, N, parallel, MAX) {
    ## This is the function order_estimate, from the original code by Youn
    ## and Simon, but I've added MAX to the arguments
    
    a <- .YS_generatea(MAX,max(rowSums(sample.gene.mutation)))
    aa <- .YS_generatea(MAX-1,max(rowSums(sample.gene.mutation)))

## Run optimization with different inital values N times using parallel computing if the variable "parallel" is TRUE :
    if(parallel)
        tmp <- foreach (kk = 1:N) %dopar%  .YS_main.function(sample.gene.mutation,a,aa, MAX)
    else ## Otherwise, run for loop
    {
        tmp <- vector("list", N)
        for(i in 1:N) tmp[[i]] <- .YS_main.function(sample.gene.mutation,a,aa, MAX)
    }

    minusloglik=rep(Inf,N)
    for(l in 1:N)
        if(is.list(tmp[[l]])) minusloglik[l]=tmp[[l]][[2]]

    result <- tmp[[which(minusloglik==min(minusloglik))[1]]][[1]] #find the one giving the maximum likelihood

    return(result)
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



## d1 <- matrix(rbinom(1000 * 14, 1, 0.1), ncol = 14)
## rownames(d1) <- paste("ID", 1:1000, sep = "")
## colnames(d1) <- paste(1:14, sep = "")

## d1[, 1] <- rbinom(1000, 1, 0.8)
## d1[ d1[, 1] == 1, 2] <- rbinom(sum(d1[, 1] == 1) , 1, 0.99)
## d1[ d1[, 1] == 0, 2] <- rbinom(sum(d1[, 1] == 0) , 1, 0.01)

## d1[ d1[, 2] == 1, 3] <- rbinom(sum(d1[, 2] == 1) , 1, 0.99)
## d1[ d1[, 2] == 0, 3] <- rbinom(sum(d1[, 2] == 0) , 1, 0.01)

## tmp <- bcbn(d1, samples = 10000)










##### Clonal ordering

run.kind.of.clonalordering <- function(x, p.v.thresh = 0.05,
                               p.adjust.method = "none", type.out =
                               "graphNEL") {
  ## This just establishes an order, since higher freq. ought to happen
  ## earlier, but the actual counts alllow for reversion

    ## This ain't the clonal ordering in the usual sense
    ## Actually, this is an example of how not to try to
    ## establish an order
    omuts <- sort(colSums(x), decreasing = TRUE)
  ngenes <- length(omuts)
  namesg <- names(omuts)
  lps <- matrix(0, nrow = choose(ngenes, 2), ncol = 3)
  lpn <- matrix("no_name", nrow = choose(ngenes, 2), ncol = 2)
  k <- 1
 
  for(i in 1:(ngenes - 1)) {
    for(j in (i + 1):ngenes) {
      tmpx <- table(x[, i], x[, j])[c(2, 3)]
      tmp.p <- binom.test(tmpx)$p.value
      if(tmpx[1] > tmpx[2]) {
        lpn[k, ] <- c(namesg[i], namesg[j])
        lps[k, ] <- c(i, j, tmp.p)
      } else {
        cat("NOTE: reversion w.r.t. to total counts\n")
        lpn[k, ] <- c(namesg[j], namesg[i])
        lps[k, ] <- c(j, i, tmp.p)
      } 
      k <- k + 1
    }
  } 
  if(p.adjust.method != "none")
    lps[, 3] <- p.adjust(lps[, 3], method = p.adjust.method)
  i.select <- which(lps[, 3] < p.v.thresh)
  lps <- lps[i.select, ]
  lpn <- lpn[i.select, ]
    ## browser()
  ## lps contains the posets
    poset.to.graph(lps[, c(1, 2)], names = c("Root", namesg),
                 addroot = TRUE,
                 type = type.out)
}

## Trying to implement logic in Barrett et al., 1999.
##  - for each subject, find order of mutations (mut.order)
##  - with data for all subjects, do a test.
##  - these data are not timed. I.e., we do not take some samples
##    before other samples.

## But we could do that: have mut.order.timed and see which comes first.


mut.order <- function(x, y) {
  ## return codes
  ## 4: conflicting evidence or no evidence. Like an NA
  ## 3: same time
  ## 1: x before y
  ## 2: y before x
  t1 <- table(x, y)

  if(min(dim(t1)) == 1) {
      ## At least one of them is constant
      if(max(dim(t1)) == 1) {
          ## Both constant
          if(max(x) > max(y))
              return(1)
          else if(max(y) > max(x))
              return(2)
          else if(max(x) == 0)
              return(4)
          else
              return(3)
      } else {
          if(dim(t1)[1] == 2) {
              if(all(y == 1))
                  return(2)
              else
                  return(1)
          } else{
              if(all(x == 1))
                  return(1)
              else
                  return(2)
          }
      }
  }
  
  if( (t1[2] > 0) && (t1[3] > 0) ) {
    return(4)
  } else if( (t1[2] > 0) && (t1[4] > 0) ) {
    return(1)
  } else if( (t1[3] > 0) && (t1[4] > 0) ) {
    return(2)
  } else if( (t1[4] > 0)) { ## and t1[3] == 0 and t1[2] == 0
    return(3)
  } else return(4)
}

indiv.clonalordering.notime <- function(x) {
  ## combs <- t(combn(ncol(x), 2))
  combs <- t(combn(sort(colnames(x)), 2))
  oo <- apply(combs, 1, function(z) mut.order(x[, z[1]], x[ , z[2]]))
  return(data.frame(combs, oo))
}


clonal.ordering <- function(x) {
    ## x is a list of observations (rows) by genes (columns)
    ## Each list is for a different subject
    ## But each list should have the same genes.

    ## Output could be a data frame, with first two cols
    ## the genes (so the gene pair) and third column and successive
    ## the evidence codes, one for each sample

    ## Still need to figure out how to sample (single cells, groups of cells)
    ## and how to filter genes
}



## wrapFitnessLinearVerbose <- function(allGenotypes_, genNum_, birthRate_,
##                               s_, numDrivers_,  v6_) {
##   .Call("wrap_fitness_linear_verbose", allGenotypes_, genNum_, birthRate_,
##         s_, numDrivers_,  v6_, PACKAGE = "MatherR")
## }


## wrapf1 <- function(inmat, dummy) {
##   .Call("f1", inmat, dummy)
##   print(inmat)
## }



## Example:
## m1 <- matrix(1:15, ncol = 3)
## wrapf1(m1, 1)




get.mut.vector.whole <- function(filename, timeSample, threshold = 0.5,
                                 remove.offending = TRUE) {
    ## FIXME: make remove.offending = FALSE. It is extremely dangerous!!!
    ## Obtain, from a file with results from a simulation run, the vector
    ## of 0/1 corresponding to each gene.
    
    ## threshold is the min. proportion for a mutation to be detected
    ## We are doing whole tumor sampling here, as in Sprouffske

    ## timeSample: do we sample at end, or at a time point, chosen
    ## randomly, from all those with at least one driver?
    
    ## We can be using rds and RData

    ## This is fragile: if things fail below, but there is a tmp object in
    ## global env, we can get garbage.

    ## Can be made more robust by assigning to another name on load,
    ## and in the future I will only use rds
    ## or if RData, do as in ADaCGH2 with named RDatas
    rt <- try({
        if(length(grep("RData$", filename))) {
            load(filename) ## we expect the object to be called tmp
        } else {
            tmp <- readRDS(filename)
        }
    })
    if(inherits(rt, "try-error")) {
        if(remove.offending) {
            warning("Removing offending filename ", filename)
            file.remove(filename)
            return(NA)
        }
    } else {
        
        if(timeSample == "last") {
            return(as.numeric((tcrossprod(tmp$pops.by.time[nrow(tmp$pops.by.time), -1],
                                          tmp$Genotypes)/tmp$TotalPopSize) > threshold))
        } else if (timeSample == "uniform") {
            the.time <- sample(which(tmp$PerSampleStats[, 4] > 0), 1)
            pop <- tmp$pops.by.time[the.time, -1]
            popSize <- tmp$PerSampleStats[the.time, 1]
            return( as.numeric((tcrossprod(pop, tmp$Genotypes)/popSize) > threshold) )
        }
  }
}





get.mut.vector.singlecell <- function(filename, timeSample,
                                      remove.offending = TRUE) {
    ## FIXME: make remove.offending = FALSE. It is extremely dangerous!!!
    ## Obtain, from a file with results from a simulation run, the vector
    ## of 0/1 corresponding to each gene.
    
    ## No threshold, as single cell.

    ## timeSample: do we sample at end, or at a time point, chosen
    ## randomly, from all those with at least one driver?
    
    ## We can be using rds and RData

    ## This is fragile: if things fail below, but there is a tmp object in
    ## global env, we can get garbage.

    ## Can be made more robust by assigning to another name on load,
    ## and in the future I will only use rds
    ## or if RData, do as in ADaCGH2 with named RDatas
    rt <- try({
        if(length(grep("RData$", filename))) {
            load(filename) ## we expect the object to be called tmp
        } else {
            tmp <- readRDS(filename)
        }
    })
    if(inherits(rt, "try-error")) {
        if(remove.offending) {
            warning("Removing offending filename ", filename)
            file.remove(filename)
            return(NA)
        }
    } else {
        if(timeSample == "last") {
            the.time <- nrow(tmp$pops.by.time)
        } else if (timeSample == "uniform") {
            the.time <- sample(which(tmp$PerSampleStats[, 4] > 0), 1)
        }
        pop <- tmp$pops.by.time[the.time, -1]
        ##       popSize <- tmp$PerSampleStats[the.time, 1]
        ## genot <- sample(seq_along(pop), 1, prob = pop)
        return(tmp$Genotypes[, sample(seq_along(pop), 1, prob = pop)])
  }
}


## simulate from a poset. Like in non-timed OTs
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
