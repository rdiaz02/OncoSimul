## Copyright 2013, 2014 Ramon Diaz-Uriarte

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


samplePop <- function(x, timeSample = "last", typeSample = "whole",
                      thresholdWhole = 0.5) {
    if(inherits(x, "oncosimulpop"))
        z <- do.call(rbind,
                     lapply(x,
                            get.mut.vector,
                            timeSample = timeSample,
                            typeSample = typeSample,
                            thresholdWhole = thresholdWhole))
    else {
        z <- get.mut.vector(x,
                            timeSample = timeSample,
                            typeSample = typeSample,
                            thresholdWhole = thresholdWhole)
        dim(z) <- c(1, length(z))
    }
    cat("\n Subjects by Genes matrix of ",
        nrow(z), " subjects and ",
        ncol(z), " genes:\n")
    colnames(z) <- paste0("G.", seq_len(ncol(z)))
    return(z)
}


oncoSimulPop <- function(Nindiv,
                         poset,
                         model = "Bozic",
                         numPassengers = 30,
                         mu = 5e-7,
                         detectionSize = 1e6,
                         detectionDrivers = 4,
                         sampleEvery = 5,
                         initSize = 2000,
                         s = 0.1,
                         sh = -1,
                         K = initSize/(exp(1) - 1),
                         keepEvery = sampleEvery,
                         finalTime = 0.25 * 25 * 365,
                         onlyCancer = TRUE,
                         max.memory = 2000,
                         max.wall.time = 200,
#                         endTimeEvery = -9,
                         silent = TRUE,
                         mc.cores = detectCores()) {

    if(.Platform$OS.type == "windows") {
        if(mc.cores != 1)
            message("You are running Windows. Setting mc.cores = 1")
        mc.cores <- 1
    }
    pop <- mclapply(seq.int(Nindiv),
                    function(x)
                    oncoSimulIndiv(
                        poset = poset,
                        model = model,
                        numPassengers = numPassengers,
                        mu = mu,
                        detectionSize = detectionSize,
                        detectionDrivers = detectionDrivers,
                        sampleEvery = sampleEvery,
                        initSize = initSize,
                        s = s,
                        sh = sh,
                        K = K,
                        keepEvery = keepEvery,
                        finalTime = finalTime,
                        onlyCancer = onlyCancer,
                        max.memory = max.memory,
                        max.wall.time = max.wall.time,
##                        endTimeEvery = endTimeEvery,
                        silent = silent),
                    mc.cores = mc.cores
                    )
    class(pop) <- "oncosimulpop"
    attributes(pop)$call <- match.call()
    return(pop)
}

## log( (K+N)/K  ) = 1; k + n = k * exp(1); k(exp - 1) = n; k = n/(exp - 1)

## pass an adj matrix.
## offer convert poset to adj matrix.

oncoSimulIndiv <- function(poset,
                           model = "Bozic",
                           numPassengers = 30,
                           mu = 5e-7,
                           detectionSize = 1e6,
                           detectionDrivers = 4,
                           sampleEvery = 5,
                           initSize = 2000,
                           s = 0.1,
                           sh = -1,
                           K = initSize/(exp(1) - 1),
                           keepEvery = sampleEvery,
                           finalTime = 0.25 * 25 * 365,
                           onlyCancer = TRUE,
                           max.memory = 2000,
                           max.wall.time = 200,
##                           endTimeEvery = -9,
                           silent = TRUE
                           ) {
    call <- match.call()
    rt <- poset.to.restrictTable(poset)


    numDrivers <- nrow(rt)
    numGenes <- (numDrivers + numPassengers)
    
    if(numGenes > 64)
        stop("Largest possible number of genes is 64")

    if(keepEvery < sampleEvery)
        warning("setting keepEvery to sampleEvery")


    ## legacies from poor name choices
    typeFitness <- switch(model,
                          "Bozic" = "bozic1",
                          "Exp" = "exp",
                          "McFarlandLog" = "mcfarlandlog",
                          "McFL" = "mcfarlandlog",
                          stop("No valid value for model")
                          )

    if(typeFitness == "exp") {
        death <- 1
        mutatorGenotype <- 1
    } else {
        death <- -99
        mutatorGenotype <- 0
    }
    birth <- -99

    if( (typeFitness == "mcfarlandlog") &&
       (sampleEvery > 0.05)) {
        warning("With the McFarland model you generally want smaller sampleEvery")
    }
    
    if(typeFitness == "mcfarlandlog") {
        endTimeEvery <- keepEvery
    } else {
        endTimeEvery <- -9
    }
    ## endTimeEvery <- -9


    ## A simulation stops if cancer or finalTime appear, the first
    ## one. But if we set onlyCnacer = FALSE, we also accept simuls
    ## without cancer (or without anything)
    
    doneSimuls <- FALSE
    while(!doneSimuls) {
        op <- try(oncoSimul.internal(restrict.table = rt,
                                 numGenes = numGenes,
                                 typeFitness = typeFitness,
                                 typeCBN = "-99",
                                 birth = birth,
                                 s = s,
                                 sh = sh,
                                 death = death,  
                                 mu =  mu,  
                                 initSize =  initSize, 
                                 sampleEvery =  sampleEvery,  
                                 detectionSize =  detectionSize, 
                                 mutatorGenotype = mutatorGenotype,  
                                 finalTime = finalTime, 
                                 initSize_species = 2000, 
                                 initSize_iter = 500, 
                                 seed_gsl = NULL, 
                                 verbosity = 1, 
                                 initMutant = -1, 
                                 speciesFS = 40000,  
                                 ratioForce = 2,  
                                 max.memory = max.memory, 
                                 max.wall.time = max.wall.time, 
                                 keepEvery = keepEvery,  
                                 alpha = 0.0015,  
                                 K = K, 
                                 endTimeEvery = endTimeEvery, 
                                 finalDrivers = detectionDrivers, 
                                 silent = silent ),
                  silent = silent)

        if(!inherits(op, "try-error")) {
            if(!silent) {
                cat("\n ... finished this run:")
                cat("\n       Total Pop Size = ", op$TotalPopSize)
                cat("\n       Drivers Last = ", op$MaxDriversLast)
                cat("\n       Final Time = ", op$FinalTime)
                ## cat("\n       Numerical Issuses?", op$ti_dbl_min)
                
                ##cat("\n")
            }
            ## browser()
            if(onlyCancer) {
                doneSimuls <- reachCancer(op, ndr = detectionDrivers,
                                             detectionSize = detectionSize,
                                             maxPopSize = 1e15)
            } else {
                doneSimuls <- TRUE
            }
            if(!silent) {
                if(doneSimuls)
                    cat("\n ... Keeping this one\n")
                else
                    cat("\n ... Cancer not reached\n")
            }
        } else {
            if(length(grep("BAIL OUT NOW", op)))
                stop("Unrecoverable error")
            if(!silent)
                cat("\n Simulation aborted because of numerical or other problems.",
                    "Proceeding to next one.\n")
        }
    }
    class(op) <- "oncosimul"
    attributes(op)$call <- call
    return(op)
}


summary.oncosimul <- function(object, ...) {
    tmp <- object[c("NumClones", "TotalPopSize", "LargestClone",
               "MaxNumDrivers", "MaxDriversLast",
               "NumDriversLargestPop", "TotalPresentDrivers",
               "FinalTime", "NumIter", "HittedWallTime")]
    tmp$errorMF <- object$other$errorMF
    if(tmp$errorMF == -99) tmp$errorMF <- NA
    tmp$OccurringDrivers <- object$OccurringDrivers
    return(as.data.frame(tmp))
}

print.oncosimul <- function(x, ...) {
    cat("\nIndividual OncoSimul trajectory with call:\n ")
    print(attributes(x)$call)
    cat("\n")
    print(summary(x))
}

## I want this to return things storable
summary.oncosimulpop <- function(object, ...) {
    as.data.frame(rbindlist(lapply(object, summary)))
}

print.oncosimulpop <- function(x, ...) {
    cat("\nPopulation of OncoSimul trajectories of ",
        length(x), " individuals. Call :\n")
    print(attributes(x)$call)
    cat("\n")
    print(summary(x))
}


plot.oncosimulpop <- function(x, ask = TRUE,
                              col = c(8, "orange", 6:1),
                              log = "y",
                              ltyClone = 2:6,
                              lwdClone = 0.2,
                              ltyDrivers = 1,
                              lwdDrivers = 3,
                              xlab = "Time units",
                              ylab = "Number of cells",
                              plotClones = TRUE,
                              plotDrivers = TRUE,
                              addtot = FALSE,
                              addtotlwd = 0.5,
                              yl = NULL,
                              thinData = FALSE,
                              thinData.keep = 0.1,
                              thinData.min = 2,
                              ...
                              ) {
    op <- par(ask = ask)
    on.exit(par(op))
    null <- lapply(x, function(z)
                   plot.oncosimul(z,
                          col = col,
                          log = log,
                          ltyClone = ltyClone,
                          lwdClone = lwdClone,
                          ltyDrivers = ltyDrivers,
                          lwdDrivers = lwdDrivers,
                          xlab = xlab,
                          ylab = ylab,
                          plotClones = plotClones,
                          plotDrivers = plotDrivers,
                          addtot = addtot,
                          addtotlwd = addtotlwd,
                          yl = yl,
                          thinData = thinData,
                          thinData.keep = thinData.keep,
                          thinData.min = thinData.min,
                          ...))
}
                            


plot.oncosimul <- function(x, col = c(8, "orange", 6:1),
                           log = "y",
                           ltyClone = 2:6,
                           lwdClone = 0.2,
                           ltyDrivers = 1,
                           lwdDrivers = 3,
                           xlab = "Time units",
                           ylab = "Number of cells",
                           plotClones = TRUE,
                           plotDrivers = TRUE,
                           addtot = FALSE,
                           addtotlwd = 0.5,
                           yl = NULL,
                           thinData = FALSE,
                           thinData.keep = 0.1,
                           thinData.min = 2,
                           ...
                           ) {

    if(thinData)
        x <- thin.pop.data(x, keep = thinData.keep, min.keep = thinData.min)
    
    ndr <- apply(x$Genotypes[1:x$NumDrivers, , drop = FALSE], 2, sum)

    if(is.null(yl)) {
        if(log %in% c("y", "xy", "yx") )
            yl <- c(1, max(apply(x$pops.by.time[, -1, drop = FALSE], 1, sum)))
        else
            yl <- c(0, max(apply(x$pops.by.time[, -1, drop = FALSE], 1, sum)))
    }
    if(plotClones) {
        plotClones(x,
                   ndr = ndr, 
                   xlab = xlab,
                   ylab = ylab,
                   lty = ltyClone,
                   col = col, 
                   ylim = yl,
                   lwd = lwdClone,
                   axes = FALSE,
                   log = log,
                   ...)
    }

    if(plotClones && plotDrivers)
        par(new = TRUE)
    
    if(plotDrivers){
        plotDrivers0(x,
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
}


poset2AdjMat <- function(x) {
    return(poset.to.graph(x, names = 1:max(x), addroot = FALSE,
                          type = "adjmat"))
}


## adjMatNoDeps <- function(ngenes = 11) {
##     posetToAdj(cbind(ngenes, 0))
## }

## posetNoDeps <- function(ngenes = 11){
##     cbind(0, ngenes)
## }


plotPoset <- function(x, names = NULL, addroot = FALSE,
                       box = FALSE, ...) {
  if(is.null(names)) {
    if(addroot) names <- c("Root", 1:max(x))
    else names <- 1:max(x)
  }
  plot(poset.to.graph(x, names, addroot), ...)
  if(box)
    box()
}

## plotAdjMat <- function(x) {
##     plot(as(x, "graphNEL"))
## }


############# The rest are internal functions


get.mut.vector.whole <- function(tmp, timeSample = "last", threshold = 0.5) {
    ## Obtain, from  results from a simulation run, the vector
    ## of 0/1 corresponding to each gene.
    
    ## threshold is the min. proportion for a mutation to be detected
    ## We are doing whole tumor sampling here, as in Sprouffske

    ## timeSample: do we sample at end, or at a time point, chosen
    ## randomly, from all those with at least one driver?
    
        
    if(timeSample == "last") {
        return(as.numeric((tcrossprod(tmp$pops.by.time[nrow(tmp$pops.by.time), -1],
                                      tmp$Genotypes)/tmp$TotalPopSize) > threshold))
    } else if (timeSample %in% c("uniform", "unif")) {
        the.time <- sample(which(tmp$PerSampleStats[, 4] > 0), 1)
        pop <- tmp$pops.by.time[the.time, -1]
        popSize <- tmp$PerSampleStats[the.time, 1]
        return( as.numeric((tcrossprod(pop, tmp$Genotypes)/popSize) > threshold) )
    }
}


get.mut.vector.singlecell <- function(tmp, timeSample = "last") {
    ## No threshold, as single cell.

    ## timeSample: do we sample at end, or at a time point, chosen
    ## randomly, from all those with at least one driver?
    
    if(timeSample == "last") {
        the.time <- nrow(tmp$pops.by.time)
    } else if (timeSample %in% c("uniform", "unif")) {
        the.time <- sample(which(tmp$PerSampleStats[, 4] > 0), 1)
    }
    pop <- tmp$pops.by.time[the.time, -1]
    ##       popSize <- tmp$PerSampleStats[the.time, 1]
    ## genot <- sample(seq_along(pop), 1, prob = pop)
    return(tmp$Genotypes[, sample(seq_along(pop), 1, prob = pop)])
}


get.mut.vector <- function(x, timeSample = "whole", typeSample = "last",
                           thresholdWhole = 0.5) {
    if(typeSample %in% c("wholeTumor", "whole")) {
        get.mut.vector.whole(x, timeSample = timeSample,
                             threshold = thresholdWhole)
    } else if(typeSample %in%  c("singleCell", "single")) {
        get.mut.vector.singlecell(x, timeSample = timeSample)
    }
}







reachCancer <- function(x, ndr = 0, detectionSize = 0,
                        maxPopSize = 1e15) {
    return(
        ( ((x$TotalPopSize >= detectionSize) ||
           (x$MaxDriversLast >= ndr)) &&
         ## (x$ti_dbl_min == 0) && ## silly, since now impossible
         (x$TotalPopSize < maxPopSize) ## numerical issues here
         ))
}


## adjM.or.poset.to.restrictTable <- function(x) {
##     if(nrow(x) == 1) {
##         is.poset <- TRUE 
##     }
##     if(nrow(x) == ncol(x))
    
##     if(all(colnames(x) == rownames(x)) ||
##        (ncol(x) > 2) ) {
##         if(any(! (x %in% c(0, 1) )))
##             stop("This looks like an adjacency matrix, but entries are not 0 or 1")
        
##     }
    
##     if(ncol(x) == 2) {
##         if(colnames(x) == c("Ancestor"))

##     }

## }


oncoSimul.internal <- function(restrict.table,
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
                      speciesFS = 40000,
                      ratioForce = 2,
                      max.memory = 20000,
                      max.wall.time = 3600,
                      keepEvery = 20,
                      alpha = 0.0015,
                      K = 1000,
                      endTimeEvery = NULL,
                      finalDrivers = 1000,
                      silent = TRUE) {
  ## the value of 20000, in megabytes, for max.memory sets a limit of ~ 20 GB
  
  ## FIXME: check argument types for typeFitness 

    ## FIXME: keepEvery not a multiple of sampleEvery

  if(initSize_species < 10) {
    warning("initSize_species too small?")
  }
  if(initSize_iter < 100) {
    warning("initSize_iter too small?")
  }
  if(keepEvery < sampleEvery)
    warning("setting keepEvery to sampleEvery")
  if(is.null(seed_gsl)) {## passing a null creates a random seed
    seed_gsl <- as.integer(round(runif(1, min = 0, max = 2^16)))
    if(!silent)
      cat(paste("\n Using ", seed_gsl, " as seed for GSL\n"))
  }

  numDrivers <- nrow(restrict.table)
  if(length(unique(restrict.table[, 1])) != numDrivers)
    stop("BAIL OUT NOW: EH??!! length(unique(restrict.table[, 1])) != numDrivers)")
  ddr <- restrict.table[, 1]
  if(any(diff(ddr) != 1))
    stop("BAIL OUT NOW:  any(diff(ddr) != 1")
  ## sanity checks
  if(max(restrict.table[, 1]) != numDrivers)
    stop("BAIL OUT NOW: max(restrict.table[, 1]) != numDrivers")
  if(numDrivers > numGenes)
    stop("BAIL OUT NOW: numDrivers > numGenes")
  
  non.dep.drivers <- restrict.table[which(restrict.table[, 2] == 0), 1]


  if( (typeFitness == "bozic1") && (mutatorGenotype) )
    warning("Using fitness bozic1 with mutatorGenotype;",
            "this will have no effect.")

  if( (typeFitness == "exp") && (death != 1) )
    warning("Using fitness exp with death != 1")


  if( (is.null(endTimeEvery) || (endTimeEvery > 0)) &&
      (typeFitness %in% c("bozic1", "exp") )) {
          warning(paste("endTimeEvery will take a positive value. ",
                        "This will make simulations not stop until the next ",
                        "endTimeEvery has been reached. Thus, in simulations ",
                        "with very fast growth, simulations can take a long ",
                        "time to finish, or can hit the wall time limit. "))
      }
  if(is.null(endTimeEvery))
    endTimeEvery <- keepEvery
  if( (endTimeEvery > 0) && (endTimeEvery %% keepEvery) )
    warning("!(endTimeEvery %% keepEvery)")
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
    stop("BAIL OUT NOW: Negative dependencies in restriction table")

  ## transpose the table
  rtC <- convertRestrictTable(restrict.table)

     
  ## return the matching call? call <- match.call()
  ## and then return(c(.Call(), call))
  call <- match.call()
  return(c(.Call("C_Algorithm5",
      ##"Algorithm5",
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
                 keepEvery,
                 alpha,
                 sh,
                 K,
                 endTimeEvery,
                 finalDrivers),
##               PACKAGE = "OncoSimulR"),
##           call = call,
           NumDrivers = numDrivers
##         ,  initMutant = initMutant
           ))
}


## colnames.to.pops.by.time <- function(pops.by.time) {
##   if(prod(dim(pops.by.time)) > 1) {  
##     ## colnames(pops.by.time) <- rep("", ncol(pops.by.time))
##     colnames(pops.by.time) <- c("Time",
##                                 paste("Clone_", 1:tmp$NumClones, sep = ""))
##   }
## }


create.muts.by.time <- function(tmp) { ## tmp is the output from Algorithm5
  if(tmp$NumClones > 1) {
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
  if(tmp$NumClones > 1) {
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



## For plotting, this helps decrease huge file sizes, while still showing
## the start of each clone, if it was originally recorded.

thin.pop.data <- function(x, keep = 0.1, min.keep = 3) {
    norig <- nrow(x$pops.by.time)
    keep1 <- round(seq.int(from = 1, to = norig,
                           length.out = round(norig * keep)))
    keep2 <- apply(x$pops.by.time[, -1, drop = FALSE],
                   1, function(x) any((x[x > 0] < min.keep)))
    keep <- sort(union(keep1, keep2))
    x$pops.by.time <- x$pops.by.time[keep, , drop = FALSE]
    return(x)
}



plotClones <- function(z, ndr = NULL, na.subs = TRUE,
                    log = "y", type = "l",
                    lty = 1:8, col = 1:9, ...) {

    ## if given ndr, we order columns based on ndr, so clones with more
    ## drivers are plotted last

    y <- z$pops.by.time[, 2:ncol(z$pops.by.time), drop = FALSE]
    
    if(na.subs){
        y[y == 0] <- NA
    }
  if(!is.null(ndr)) {
      ## could be done above, to avoid creating
      ## more copies
      oo <- order(ndr)
      y <- y[, oo, drop = FALSE]
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


rtNoDep <- function(numdrivers) {
  ## create a restriction table with no dependencies
  x <- matrix(nrow = numdrivers, ncol = 3)
  x[, 1] <- 1:numdrivers
  x[, 2] <- 0
  x[, 3] <- -9
  return(x)
}



convertRestrictTable <- function(x) {
  t.restrictTable <- matrix(as.integer(x),
                            ncol = nrow(x), byrow = TRUE)

  t.restrictTable[-2, ] <- t.restrictTable[-2, ] - 1
  return(t.restrictTable)
}





adjmat.to.restrictTable <- function(x) {
    ## we have the zero
    ## x <- x[-1, -1]
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











## get.mut.vector.whole <- function(filename, timeSample, threshold = 0.5,
##                                  remove.offending = TRUE) {
##     ## FIXME: make remove.offending = FALSE. It is extremely dangerous!!!
##     ## Obtain, from a file with results from a simulation run, the vector
##     ## of 0/1 corresponding to each gene.
    
##     ## threshold is the min. proportion for a mutation to be detected
##     ## We are doing whole tumor sampling here, as in Sprouffske

##     ## timeSample: do we sample at end, or at a time point, chosen
##     ## randomly, from all those with at least one driver?
    
##     ## We can be using rds and RData

##     ## This is fragile: if things fail below, but there is a tmp object in
##     ## global env, we can get garbage.

##     ## Can be made more robust by assigning to another name on load,
##     ## and in the future I will only use rds
##     ## or if RData, do as in ADaCGH2 with named RDatas
##     rt <- try({
##         if(length(grep("RData$", filename))) {
##             load(filename) ## we expect the object to be called tmp
##         } else {
##             tmp <- readRDS(filename)
##         }
##     })
##     if(inherits(rt, "try-error")) {
##         if(remove.offending) {
##             warning("Removing offending filename ", filename)
##             file.remove(filename)
##             return(NA)
##         }
##     } else {
        
##         if(timeSample == "last") {
##             return(as.numeric((tcrossprod(tmp$pops.by.time[nrow(tmp$pops.by.time), -1],
##                                           tmp$Genotypes)/tmp$TotalPopSize) > threshold))
##         } else if (timeSample == "uniform") {
##             the.time <- sample(which(tmp$PerSampleStats[, 4] > 0), 1)
##             pop <- tmp$pops.by.time[the.time, -1]
##             popSize <- tmp$PerSampleStats[the.time, 1]
##             return( as.numeric((tcrossprod(pop, tmp$Genotypes)/popSize) > threshold) )
##         }
##   }
## }





## get.mut.vector.singlecell <- function(filename, timeSample,
##                                       remove.offending = TRUE) {
##     ## FIXME: make remove.offending = FALSE. It is extremely dangerous!!!
##     ## Obtain, from a file with results from a simulation run, the vector
##     ## of 0/1 corresponding to each gene.
    
##     ## No threshold, as single cell.

##     ## timeSample: do we sample at end, or at a time point, chosen
##     ## randomly, from all those with at least one driver?
    
##     ## We can be using rds and RData

##     ## This is fragile: if things fail below, but there is a tmp object in
##     ## global env, we can get garbage.

##     ## Can be made more robust by assigning to another name on load,
##     ## and in the future I will only use rds
##     ## or if RData, do as in ADaCGH2 with named RDatas
##     rt <- try({
##         if(length(grep("RData$", filename))) {
##             load(filename) ## we expect the object to be called tmp
##         } else {
##             tmp <- readRDS(filename)
##         }
##     })
##     if(inherits(rt, "try-error")) {
##         if(remove.offending) {
##             warning("Removing offending filename ", filename)
##             file.remove(filename)
##             return(NA)
##         }
##     } else {
##         if(timeSample == "last") {
##             the.time <- nrow(tmp$pops.by.time)
##         } else if (timeSample == "uniform") {
##             the.time <- sample(which(tmp$PerSampleStats[, 4] > 0), 1)
##         }
##         pop <- tmp$pops.by.time[the.time, -1]
##         ##       popSize <- tmp$PerSampleStats[the.time, 1]
##         ## genot <- sample(seq_along(pop), 1, prob = pop)
##         return(tmp$Genotypes[, sample(seq_along(pop), 1, prob = pop)])
##   }
## }






