wrapTi <- function(v1_, v2_, v3_, v4_, v5_,  v6_) {
  .Call("wrap_ti", v1_, v2_, v3_, v4_, v5_, v6_, PACKAGE = "MatherR")
}

wrapAlgo2 <- function(v1_, v2_, v3_, v4_, v5_, v6_, vd_, v7_) {
  .Call("wrap_Algo2", v1_, v2_, v3_, v4_, v5_, v6_,  vd_, v7_, PACKAGE = "MatherR")
}

wrapAlgo3 <- function(v1_, v2_, v3_, v4_, v5_,  v6_, vd_, v7_) {
  .Call("wrap_Algo3", v1_, v2_, v3_, v4_, v5_, v6_, vd_, v7_, PACKAGE = "MatherR")
}

wrapFitnessLinearRcpp <- function(allGenotypes_, genNum_, birthRate_,
                              s_, numDrivers_,  v6_) {
  .Call("wrap_fitness_linear_Rcpp", allGenotypes_, genNum_, birthRate_,
        s_, numDrivers_,  v6_, PACKAGE = "MatherR")
}

wrapFitnessCBNRcpp <- function(mutatedPos_, genotypes_,
                           genNum_, restrictTable_,
                           numDrivers_,
                           birthRate_, s_,
                           fitnessParent_, typeCBN_, retval_){
  
  .Call("wrap_fitness_CBN_Rcpp",mutatedPos_, genotypes_,
                           genNum_, restrictTable_,
                           numDrivers_,
                           birthRate_, s_,
                           fitnessParent_, typeCBN_, retval_,
        PACKAGE = "MatherR")
}

wrapFitnessCBNArma <- function(mutatedPos_, genotypes_,
                           genNum_, restrictTable_,
                           numDrivers_,
                           birthRate_, s_,
                           fitnessParent_, typeCBN_, retval_){
  
  .Call("wrap_fitness_CBN_Arma",mutatedPos_, genotypes_,
                           genNum_, restrictTable_,
                           numDrivers_,
                           birthRate_, s_,
                           fitnessParent_, typeCBN_, retval_,
        PACKAGE = "MatherR")
}

wrapFitnessCBNstd <- function(mutatedPos_, genotypes_,
                              genNum_, restrictTable_,
                              numDrivers_,
                              birthRate_, s_,
                              fitnessParent_, typeCBN_, retval_){
  
  .Call("wrap_fitness_CBN_std",mutatedPos_, genotypes_,
        genNum_, restrictTable_,
        numDrivers_,
        birthRate_, s_,
        fitnessParent_, typeCBN_, retval_,
        PACKAGE = "MatherR")
}



## FIXME: how large should initSize_species be???
Algo5 <- function(restrict.table,
                  numGenes,
                  typeCBN,
                  birth, 
                  s, 
                  death,
                  mu,
                  initSize,
                  sampleEvery,
                  detectionSize,
                  finalTime = 1e90,
                  initSize_species = 2000,
                  initSize_iter = 500,
                  seed_gsl = NULL,
                  verbosity = 1,
                  func = "D",
                  typeFitness = "bozic",
                  ## next two set so no forced sampling by default
                  speciesFS = 40000,
                  ratioForce = 2) {

  ## FIXME: check argument types for typeFitness and typeCBN
  
  if(initSize_species < 10) {
    warning("initSize_species too small?")
  }
  if(initSize_iter < 100) {
    warning("initSize_iter too small?")
  }
  
  if(is.null(seed_gsl)) {## passing a null creates a random seed
    seed_gsl <- as.integer(round(runif(1, min = 0, max = 2^16)))
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

  tt <- system.time({
    if(func == "A") {
      tmp <- .Call("Algorithm5",
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
                   initSize_species,
                   initSize_iter,
                   seed_gsl,
                   verbosity,
                   PACKAGE = "MatherR")
    } else if(func == "B") {
      tmp <- .Call("Algorithm5B",
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
                   initSize_species,
                   initSize_iter,
                   seed_gsl,
                   verbosity,
                   PACKAGE = "MatherR")
    } else if(func == "C") {
      tmp <- .Call("Algorithm5C",
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
                   PACKAGE = "MatherR")
                                        # tmp$Extinction <- "NA in this method"
    } else if(func == "D") {
      tmp <- .Call("Algorithm5D",
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
                   PACKAGE = "MatherR")
    } else if(func == "E") {
      tmp <- .Call("Algorithm5E",
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
                   PACKAGE = "MatherR")
    } else if(func == "F") {
      tmp <- .Call("Algorithm5F",
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
                   PACKAGE = "MatherR")
    } else if(func == "G") {
      tmp <- .Call("Algorithm5G",
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
                   PACKAGE = "MatherR")
    } else if(func == "H") {
      tmp <- .Call("Algorithm5H",
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
                   PACKAGE = "MatherR")
    } else if(func == "I") {
      tmp <- .Call("Algorithm5I",
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
                   PACKAGE = "MatherR")
    } else if(func == "J") {
      tmp <- .Call("Algorithm5J",
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
                   PACKAGE = "MatherR")
    }  else if(func == "K") {
      tmp <- .Call("Algorithm5K",
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
                   PACKAGE = "MatherR")
    }
    
  })[3]

  cat("\n  C++ code execution took ", tt, "\n")
  
  pops.by.time <- t(tmp$outNS)

  ## with K, this would work
  ## pops.by.time <- reshape(tmp$longOut, timevar = "index",
  ##                         direction = "wide", idvar = "genot")
  ## but that is crazy
  
  ## if(sampling) {
  ##   TotalPopSize <- tmp$TotalPopSize
  ##   if(TotalPopSize == 0) return(NA)

  ##   if(is.null(timeSample))
  ##     timeSample <- nrows(pops.by.time)
  ##   else
  ##     timeSample <- 
  ##   pops <- pops.by.time[ , ]
    
  ##   index <- sample()
    
  ## } else {

  ## FIXME00: most of this should be done in C++!!!
  ## this could come right from C++
    NumMutations <- apply(tmp$Genotypes, 1, sum)
    
    ## just for debugging
    if(is.matrix(tmp$Params)) {
      colnames(tmp$Params) <- c("Flag",
                                "birth", "popSize",
                                "timeLastUpdate", "W", "R",
                                "tis", "nextMutationTime",
                                "mutation")
    }
    colnames(pops.by.time) <- rep("", ncol(pops.by.time))
    colnames(pops.by.time) <- c("Iteration", "Time",
                                paste("Sp_", 1:tmp$NumSpecies, sep = ""))
    
    
    if(tmp$NumSpecies > 1) {
      muts.by.time <- cbind(pops.by.time[, c(1, 2), drop = FALSE] ,
                            t(apply(pops.by.time[, -c(1, 2), drop = FALSE], 1,
                                    function(x) tapply(x, NumMutations, sum))))
    } else {
      muts.by.time <- pops.by.time
    }
    
    CountNumDrivers <- apply(tmp$Genotypes[, 1:numDrivers, drop = FALSE], 1, sum)
    if(tmp$NumSpecies > 1) {
      if(length(unique(CountNumDrivers )) > 1) {
        drivers.by.time <- cbind(pops.by.time[, c(1, 2), drop = FALSE] ,
                                 t(apply(pops.by.time[, -c(1, 2), drop = FALSE], 1,
                                         function(x) tapply(x, CountNumDrivers, sum)))) 
      } else {
        drivers.by.time <- cbind(pops.by.time[, c(1, 2), drop = FALSE] ,
                                 rowSums(pops.by.time[, -c(1, 2), drop = FALSE]))
      }
      colnames(drivers.by.time) <- c("Iteration", "Time",
                                   paste("dr_", colnames(drivers.by.time)[-c(1,2)],
                                         sep = ""))
    } else {
      drivers.by.time <- NULL
    }
    
    colnames(muts.by.time)[c(1, 2)] <- c("Iteration", "Time")
    
    return(list(pops.by.time = pops.by.time,
                muts.by.time = muts.by.time,
                drivers.by.time = drivers.by.time,
                max.num.drivers = max(CountNumDrivers),
                rest = tmp))
##  }
}


sampleZZ <- function(zz, seed = "auto"){
  if(seed == "auto") {
    ## paste as numeric(hostname()) y el segundo from time
  }

  fname <- paste(fileroot, randomstgring)
  save(...., compress = FALSE)

  
}


plotPop <- function(z, na.subs = TRUE, log = "y", type = "l",
                    lty = 1:8, col = 1:9, ...) {
  y <- z$pops.by.time[, 3:ncol(z$pops.by.time)]
  if(na.subs){
    y[y == 0] <- NA
  }
    
  matplot(x = z$pops.by.time[, 2],
          y = y,
          log = log, type = type,
          ...)
}
## for plotting, substitute 0s by NAs

plotDrivers <- function(z, na.subs = TRUE, log = "y", type = "l",
                        lty = 1:9, col = c(8, "orange", 6:1),
                        lwd = 2, ...) {
  y <- z$drivers.by.time[, 3:ncol(z$drivers.by.time)]
  if(na.subs){
    y[y == 0] <- NA
  }

  matplot(x = z$drivers.by.time[, 2],
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


selectSB <- function(x, threshold = 5, maxGenes = 14, weighted = TRUE) {
  ## x: the genotypes, with columns as genotypes, rows as genes
  ## threshold is percentage here.
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
  plot(g1)
  ## print(xx)
  lc <- largest.cliques(g1)
  if(length(lc) > 1) {
    cat("\n WARNING: more than one largest clique\n")
  }
  lcx <- lc[[1]]
  if(length(lcx) > maxGenes)
    lcx <- lcx[order(diag(xx)[lcx], decreasing = TRUE)][1:maxGenes]
  return(lcx)
}




selectPercent <- function(x, threshold = 5, maxGenes = 14) {
  xx <- rowSums(x)
  threshold <- ceiling((threshold/100) * ncol(x))
  lcx <- which(xx >= threshold)
  if(length(lcx) > maxGenes)
    lcx <- lcx[order(xx[lcx], decreasing = TRUE)][1:maxGenes]
  return(lcx)
}

selectGenes <- function(x, threshold, method, maxGenes = 14) {
  if(method == "SB")
    return(selectSB(x, threshold, maxGenes = maxGenes))
  if(method == "Percent")
    return(selectPercent(x, threshold, maxGenes = maxGenes))
}


veltcs <- function(x) {
  ## return a vector from edge list of transitive closure
  if(length(x) > 1) {
    if(class(x[[1]]) == "graphNEL")
      tc0 <- lapply(x, transitive.closure)
    else ## adjacency matrix
      tc0 <- lapply(x, transClos)
    el <- unlist(lapply(tc0, edgeList), recursive = FALSE)

  } else{
    if(class(x) == "graphNEL")
      tc0 <- transitive.closure(x)
    else ## adjacency matrix
      tc0 <- transClos(x)
    el <- edgeList(tc0)
  }
  ell <- sapply(el, function(x) paste(x[1], x[2], sep = "_"))
  return(unique(ell)) ## important with mixtures
}

metrics1 <- function(x, y) {
  ## My PFD and PND
  ## x is model
  ## y is true
  ex <- veltcs(x)
  ey <- veltcs(y)
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



call.external.cbn <- function(data, file = "testcbn", eparam = 0.05,
                              temp = 10, steps = 200, silent = FALSE,
                              init.poset = TRUE) {
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
  zzz <- system(paste("export OMP_NUM_THREADS=", detectCores(), sep = ""),
                intern = FALSE)
  dir.create(file)
  
  if(init.poset) {
    write.linear.poset(data, file)
  } else { ## Use ct-cbn to search and create starting poset;
    ## possibly eternal. NOT RECOMMENDED
    writeLines(as.character(c(ncol(data), 0)),
               con = paste(file, ".poset", sep = ""))
    ## First create the lambda file
    zzz <- system(paste("h-cbn -f",  file, "-w"), intern = silent)
    cat("\n\n")
    ## this call requires a lambda file
    zzz <- system(paste("h-cbn -f",  file, "-e", eparam,
                 "-w -m"), intern = silent)
    cat("\n\n")
  }
  zzz <- system(paste("h-cbn -f",  file, "-s", 
                      "-T", temp,  "-N", steps,
                      "-m -w"), intern = silent)
  cat("\n\n")
  ## the final poset in file/00000.poset
}

read.poset <- function(dirname) {
  tmp <- scan(paste(dirname, "/00000.poset", sep = ""))
  tmp <- tmp[-(length(tmp))]
  tmp <- tmp[-1]
  tmp <- matrix(tmp, ncol = 2, byrow = TRUE)
  return(tmp)
}

poset.to.graph <- function(x, names, type = "graphNEL") {
  m <- length(names) 
  m1 <- matrix(0, nrow = m, ncol = m)
  colnames(m1) <- names
  rownames(m1) <- names
  ## a debugging check
  if(nrow(x) > 0) {
    if(max(x) != (m - 1))
      warning("\n in poset, max(x) != (m - 1)")
  }
  if(nrow(x) > 0) {
    m1[x + 1] <- 1
  }
  if(length(names) > 1) {
    no.ancestor <- which(apply(m1, 2, function(x) all(x == 0)))[-1]
    m1[cbind(1, no.ancestor)] <- 1
  } ## o.w. do nothing
  if(type == "adjmat") return(m1)
  else if (type == "graphNEL") return(as(m1, "graphNEL"))
  ## does not show the labels
  ## else if (type == "igraph") return(graph.adjacency(m1))
}



## Later, we will want much more info recovered, such as the probs, etc.
run.cbn <- function(x, file,
                    temp = 10, steps = 200,
                    silent = FALSE,
                    type.out = "graphNEL",
                    init.poset = TRUE,
                    eparam = 0.05) {
  zzz <- call.external.cbn(x, file = file, eparam = eparam,
                           temp = temp, steps = steps,
                           silent = silent, init.poset = init.poset)
  cnames <- colnames(x)
  poset <- read.poset(file)
  gr <- poset.to.graph(poset, names = c("Root", cnames), type = type.out)
  return(gr)
}


run.oncotree <- function(x, type.out = "graphNEL") {
  m <- ncol(x)
  onco.fit <- oncotree.fit(x)
  ## root gets number 1
  onco.fit$parent$parent.num
  ## convert to adjacency matrix
  ## but remember 1 is actually the root, so we add a 1, and later rename
  p.to.child <- cbind( onco.fit$parent$parent.num[-1], 2:(ncol(x) + 1))
  m1 <- matrix(0, nrow = m + 1, ncol = m + 1)
  m1[p.to.child] <- 1
  cnames <- c("Root", colnames(x))
  colnames(m1) <- cnames
  rownames(m1) <- cnames
  if(type.out == "adjmat") return(m1)
  else if (type.out == "graphNEL") return(as(m1, "graphNEL"))
}



run.rtreemix <- function(x, K, noise = FALSE, only.graphnel = TRUE,
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



graph.to.poset <- function(x) {
  ## FIXME: this are characters, not numeric
  return(matrix(as.numeric(unlist(edgeList(x))), ncol = 2,
                byrow = TRUE))
}


adjmat.to.restrictTable <- function(x) {
  ## we have the zero
  x <- x[-1, -1]
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
  x1 <- poset.to.graph(x, c("nn", 1:max(x)), "adjmat")
  adjmat.to.restrictTable(x1)
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
