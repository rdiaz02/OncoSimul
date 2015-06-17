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





oncoSimulSample <- function(Nindiv,
                            fp,
                            model = "Exp",
                            numPassengers = 0,
                            mu = 1e-6,
                            detectionSize = round(runif(Nindiv, 1e5, 1e8)),
                            detectionDrivers = {
                                if(inherits(fp, "fitnessEffects")) {
                                    if(length(fp$drv)) {
                                        nd <- (2: round(0.75 * length(fp$drv)))
                                    } else {
                                        nd <- 0
                                    }
                                } else {
                                    nd <- (2 : round(0.75 * max(fp)))
                                }
                                sample(nd, Nindiv,
                                       replace = TRUE)
                            },
                            sampleEvery = ifelse(model %in% c("Bozic", "Exp"), 1,
                                0.025),
                            initSize = 500,
                            s = 0.1,
                            sh = -1,
                            K = initSize/(exp(1) - 1),
                            minDetectDrvCloneSz = "auto",
                            extraTime = 0,
                            finalTime = 0.25 * 25 * 365,
                            onlyCancer = TRUE,
                            max.memory = 2000,
                            max.wall.time.total = 600,
                            max.num.tries.total = 500 * Nindiv,
                            ## well, obviously they are errors
                            ## errorHitWallTime = TRUE,
                            ## errorHitMaxTries = TRUE,
                            verbosity  = 1,
                            typeSample = "whole",
                            thresholdWhole = 0.5){
    ## No longer using mclapply, because of the way we use the limit on
    ## the number of tries.
    
    ## leaving detectionSize and detectionDrivers as they are, produces
    ## the equivalente of uniform sampling. For last, fix a single number

    if(max.num.tries.total < Nindiv)
        stop(paste("You have requested something impossible: ",
                   "max.num.tries.total < Nindiv"))
    attemptsLeft <- max.num.tries.total
    attemptsUsed <- 0
    numToRun <- Nindiv
    pop <- vector(mode = "list", length = Nindiv)
    indiv <- 1
    ## FIXME! really pass params such as extraTime and minDDr as vectors,
    ## or give a warning
    params <- data.frame(seq.int(Nindiv),
                         extraTime = extraTime,
                         minDetectDrvCloneSz = minDetectDrvCloneSz,
                         detectionSize = detectionSize,
                         detectionDrivers = detectionDrivers)[, -1, drop = FALSE]

    f.out.attempts <- function() {
        message("Run out of attempts")
        return(list(
            popSummary = NA,
            popSample = NA,
            HittedMaxTries = TRUE,
            hittedWallTime = FALSE,
            UnrecoverExcept = FALSE
        ))    
    }    

    f.out.time <- function() {
        message("Run out of time")
        return(list(
            popSummary = NA,
            popSample = NA,
            HittedMaxTries = FALSE,
            HittedWallTime = TRUE,
            UnrecoverExcept = FALSE
        ))    
    }    

    f.out.attempts.cpp <- function() {
        message("Run out of attempts (in C++)")
        return(list(
            popSummary = NA,
            popSample = NA,
            HittedMaxTries = TRUE,
            hittedWallTime = FALSE,
            UnrecoverExcept = FALSE
        ))    
    }    

    f.out.time.cpp <- function() {
        message("Run out of time (in C++)")
        return(list(
            popSummary = NA,
            popSample = NA,
            HittedMaxTries = FALSE,
            HittedWallTime = TRUE,
            UnrecoverExcept = FALSE
        ))    
    }    

    f.out.unrecover.except <- function(x) {
        message("Unrecoverable exception (in C++)")
        return(list(
            popSummary = NA,
            popSample = NA,
            HittedMaxTries = NA,
            HittedWallTime = NA,
            UnrecoverExcept = TRUE,
            ExceptionMessage = x$other$ExceptionMessage
        ))    
    }    
   
    
    startTime <- Sys.time()
    while(TRUE) {
        
        possibleAttempts <- attemptsLeft - (numToRun - 1)
        ## I think I do not want a try here.
        tmp <-  oncoSimulIndiv(fp = fp,
                               model = model,
                               numPassengers = numPassengers,
                               mu = mu,
                               detectionSize = params[indiv, "detectionSize"],
                               detectionDrivers = params[indiv, "detectionDrivers"],
                               sampleEvery = sampleEvery,
                               initSize = initSize,
                               s = s,
                               sh = sh,
                               K = K,
                               minDetectDrvCloneSz = params[indiv, "minDetectDrvCloneSz"],
                               extraTime = params[indiv, "extraTime"],
                               finalTime = finalTime,
                               max.memory = max.memory,
                               max.wall.time = max.wall.time.total,
                               max.num.tries = possibleAttempts,
                               verbosity = verbosity,
                               keepEvery = -9,
                               onlyCancer = onlyCancer,
                               errorHitWallTime = TRUE,
                               errorHitMaxTries = TRUE)
        
        if(tmp$other$UnrecoverExcept) {
            return(f.out.unrecover.except(tmp))
        }
        
        pop[[indiv]] <- tmp
        numToRun <- (numToRun - 1)
        attemptsUsed <- attemptsUsed + tmp$other$attemptsUsed
        attemptsLeft <- (max.num.tries.total - attemptsUsed)
        indiv <- indiv + 1
        
        ## We need to check in exactly this order. Attempts left only
        ## matters if no remaining individuals to run. But C++ might bail
        ## out in exactly the last individual

        if(  
            (exists("HittedMaxTries", where = tmp) &&
                 tmp[["HittedMaxTries"]])  ) {
            ## in C++ code
            return(f.out.attempts.cpp())
        } else if(  
            (exists("HittedWallTime", where = tmp) &&
                 tmp[["HittedWallTime"]])  ) {
            ## in C++ code
            return(f.out.time.cpp())
        } else if( indiv > Nindiv ) {
            if(verbosity > 0)
                message(paste("Successfully sampled ", Nindiv, " individuals"))
            class(pop) <- "oncosimulpop"
            if(inherits(fp, "fitnessEffects")) {
                geneNames <- names(getNamesID(fp))
            } else {
                geneNames <- NULL
            }
            return(list(
                popSummary = summary(pop),
                popSample = samplePop(pop, typeSample = typeSample,
                    thresholdWhole = thresholdWhole, geneNames = geneNames),
                attemptsUsed = attemptsUsed,
                probCancer = Nindiv/attemptsUsed,
                HittedMaxTries = FALSE,
                HittedWallTime = FALSE,
                UnrecoverExcept = FALSE
            ))
        } else if( attemptsLeft <= 0 ) {
            return(f.out.attempts())
        } else  if( as.double(difftime(Sys.time(), startTime, units = "secs"))
                   > max.wall.time.total ) {
            return(f.out.time())
        } 
    }
}


samplePop <- function(x, timeSample = "last", typeSample = "whole",
                      thresholdWhole = 0.5,
                      geneNames = NULL) {
    gN <- geneNames
    if(inherits(x, "oncosimulpop")) {
        z <- do.call(rbind,
                     lapply(x,
                            get.mut.vector,
                            timeSample = timeSample,
                            typeSample = typeSample,
                            thresholdWhole = thresholdWhole))
        ## We need to check if the object is coming from v.2., to avoid
        ## having to force passing a vector of names
        if(is.null(gN) && (!is.null(x[[1]]$geneNames)))
            gN <- x[[1]]$geneNames
    } else {
        z <- get.mut.vector(x,
                            timeSample = timeSample,
                            typeSample = typeSample,
                            thresholdWhole = thresholdWhole)
        dim(z) <- c(1, length(z))
        if(is.null(gN) && (!is.null(x$geneNames)))
            gN <- geneNames
    }
    message("\n Subjects by Genes matrix of ",
        nrow(z), " subjects and ",
            ncol(z), " genes.\n")

    if(!is.null(gN)) {
        colnames(z) <- gN
    } else {
        colnames(z) <- seq_len(ncol(z))
    }
    return(z)
}


oncoSimulPop <- function(Nindiv,
                         fp,
                         model = "Exp",
                         numPassengers = 30,
                         mu = 1e-6,
                         detectionSize = 1e8,
                         detectionDrivers = 4,
                         sampleEvery = ifelse(model %in% c("Bozic", "Exp"), 1,
                             0.025),
                         initSize = 500,
                         s = 0.1,
                         sh = -1,
                         K = initSize/(exp(1) - 1),
                         keepEvery = sampleEvery, 
                         minDetectDrvCloneSz = "auto",
                         extraTime = 0,
                         ## used to be this
                         ## ifelse(model \%in\% c("Bozic", "Exp"), -9,
                         ##                       5 * sampleEvery),
                         finalTime = 0.25 * 25 * 365,
                         onlyCancer = TRUE,
                         max.memory = 2000,
                         max.wall.time = 200,
                         max.num.tries = 500,
                         errorHitWallTime = TRUE,
                         errorHitMaxTries = TRUE,
                         verbosity  = 0,
                         mc.cores = detectCores()) {

    if(.Platform$OS.type == "windows") {
        if(mc.cores != 1)
            message("You are running Windows. Setting mc.cores = 1")
        mc.cores <- 1
    }
    pop <- mclapply(seq.int(Nindiv),
                    function(x)
                    oncoSimulIndiv(
                        fp = fp,
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
                        minDetectDrvCloneSz = minDetectDrvCloneSz,
                        extraTime = extraTime,
                        finalTime = finalTime,
                        onlyCancer = onlyCancer,
                        max.memory = max.memory,
                        max.wall.time = max.wall.time,
                        max.num.tries = max.num.tries,
                        errorHitWallTime = errorHitWallTime,
                        errorHitMaxTries = errorHitMaxTries,
                        verbosity = verbosity),
                    mc.cores = mc.cores
                    )
    class(pop) <- "oncosimulpop"
    attributes(pop)$call <- match.call()
    return(pop)
}

## where is the default K coming from? Here:
## log( (K+N)/K  ) = 1; k + n = k * exp(1); k(exp - 1) = n; k = n/(exp - 1)

oncoSimulIndiv <- function(fp,
                           model = "Exp",
                           numPassengers = 30,
                           mu = 1e-6,
                           detectionSize = 1e8,
                           detectionDrivers = 4,
                           sampleEvery = ifelse(model %in% c("Bozic", "Exp"), 1,
                               0.025),
                           initSize = 500,
                           s = 0.1,
                           sh = -1,
                           K = initSize/(exp(1) - 1),
                           keepEvery = sampleEvery,
                           minDetectDrvCloneSz = "auto",
                           extraTime = 0,
                           ## used to be this
                           ## ifelse(model \%in\% c("Bozic", "Exp"), -9,
                           ##                     5 * sampleEvery),
                           finalTime = 0.25 * 25 * 365,
                           onlyCancer = TRUE,
                           max.memory = 2000,
                           max.wall.time = 200,
                           max.num.tries = 500,
                           errorHitWallTime = TRUE,
                           errorHitMaxTries = TRUE,
                           verbosity = 0,
                           initMutant = NULL,
                           seed = NULL
                           ) {
    call <- match.call()
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
        warning("With the McFarland model you often want smaller sampleEvery")
    }

    if(minDetectDrvCloneSz == "auto") {
        if(model %in% c("Bozic", "Exp") )
            minDetectDrvCloneSz <- 0
        else if (model %in% c("McFL", "McFarlandLog"))
            minDetectDrvCloneSz <- eFinalMf(initSize, s, detectionDrivers)
        else
            stop("Unknown model")
    }

    if(mu < 0) {
        stop("mutation rate (mu) is negative")
    }

    if(is.null(seed)) {## passing a null creates a random seed
        ## name is a legacy. This is really the seed for the C++ generator.
        seed <- as.integer(round(runif(1, min = 0, max = 2^16)))
    }
    if(verbosity >= 2)
        cat(paste("\n Using ", seed, " as seed for C++ generator\n"))

    if( (keepEvery > 0) & (keepEvery < sampleEvery)) {
        keepEvery <- sampleEvery
        warning("setting keepEvery <- sampleEvery")
    }

    if( (typeFitness == "bozic1") && (mutatorGenotype) )
        warning("Using fitness bozic1 with mutatorGenotype;",
                "this will have no effect.")

    if( (typeFitness == "exp") && (death != 1) )
        warning("Using fitness exp with death != 1")

    
    if(!inherits(fp, "fitnessEffects")) {
        if(any(unlist(lapply(list(fp, 
                                  numPassengers,
                                  s, sh), is.null)))) {
            m <- paste("You are using the old poset format.",
                       "You must specify all of poset, numPassengers",
                       "s, and sh.")
            stop(m)
            if(length(initMutant) > 1)
                stop("With the old poset, initMutant can only take a single value.")
        }
        ## if(message.v1)
        ##     message("You are using the old poset format. Consider using the new one.")
   
    
        ## A simulation stops if cancer or finalTime appear, the first
        ## one. But if we set onlyCnacer = FALSE, we also accept simuls
        ## without cancer (or without anything)
        
        op <- try(oncoSimul.internal(poset = fp, ## restrict.table = rt,
                                     ## numGenes = numGenes,
                                     numPassengers = numPassengers,
                                     typeCBN = "CBN",
                                     birth = birth,
                                     s = s,
                                     death = death,  
                                     mu =  mu,  
                                     initSize =  initSize, 
                                     sampleEvery =  sampleEvery,  
                                     detectionSize =  detectionSize, 
                                     finalTime = finalTime, 
                                     initSize_species = 2000, 
                                     initSize_iter = 500, 
                                     seed = seed, 
                                     verbosity = verbosity, 
                                     speciesFS = 40000,  
                                     ratioForce = 2,
                                     typeFitness = typeFitness,
                                     max.memory = max.memory,
                                     mutatorGenotype = mutatorGenotype,                                   
                                     initMutant = -1, 
                                     max.wall.time = max.wall.time,
                                     max.num.tries = max.num.tries,
                                     keepEvery = keepEvery,  
                                     alpha = 0.0015,  
                                     sh = sh,
                                     K = K, 
                                     minDetectDrvCloneSz = minDetectDrvCloneSz,
                                     extraTime = extraTime,
                                     detectionDrivers = detectionDrivers,
                                     onlyCancer = onlyCancer,
                                     errorHitWallTime = errorHitWallTime,
                                     errorHitMaxTries = errorHitMaxTries),
                  silent = !verbosity)
        objClass <- "oncosimul"
    } else {
        op <- try(nr_oncoSimul.internal(rFE = fp, 
                                        birth = birth,
                                        death = death,  
                                        mu =  mu,  
                                        initSize =  initSize, 
                                        sampleEvery =  sampleEvery,  
                                        detectionSize =  detectionSize, 
                                        finalTime = finalTime, 
                                        initSize_species = 2000, 
                                        initSize_iter = 500, 
                                        seed = seed, 
                                        verbosity = verbosity, 
                                        speciesFS = 40000,  
                                        ratioForce = 2,
                                        typeFitness = typeFitness,
                                        max.memory = max.memory,
                                        mutatorGenotype = mutatorGenotype,                                   
                                        initMutant = initMutant, 
                                        max.wall.time = max.wall.time,
                                        max.num.tries = max.num.tries,
                                        keepEvery = keepEvery,  
                                        alpha = 0.0015,  
                                        K = K, 
                                        minDetectDrvCloneSz = minDetectDrvCloneSz,
                                        extraTime = extraTime,
                                        detectionDrivers = detectionDrivers,
                                        onlyCancer = onlyCancer,
                                        errorHitWallTime = errorHitWallTime,
                                        errorHitMaxTries = errorHitMaxTries),
                  silent = !verbosity)
        objClass <- c("oncosimul", "oncosimul2")
    }
    if(inherits(op, "try-error")) {
        ##         if(length(grep("BAIL OUT NOW", op)))
        stop(paste("Unrecoverable error:", op ))
    }
    if(verbosity >= 2) {
        cat("\n ... finished this run:")
        cat("\n       Total Pop Size = ", op$TotalPopSize)
        cat("\n       Drivers Last = ", op$MaxDriversLast)
        cat("\n       Final Time = ", op$FinalTime, "\n")
    }

    class(op) <- objClass
    attributes(op)$call <- call
    return(op)
}

summary.oncosimul <- function(object, ...) {

    if(object$HittedWallTime || object$HittedMaxTries ||
       object$other$UnrecoverExcept)
        return(NA)
    else {
        tmp <- object[c("NumClones", "TotalPopSize", "LargestClone",
                        "MaxNumDrivers", "MaxDriversLast",
                        "NumDriversLargestPop", "TotalPresentDrivers",
                        "FinalTime", "NumIter", "HittedWallTime")]
 
        tmp$errorMF <- object$other$errorMF
        tmp$minDMratio <- object$other$minDMratio
        tmp$minBMratio <- object$other$minBMratio
        if( (tmp$errorMF == -99)) tmp$errorMF <- NA
        if( (tmp$minDMratio == -99)) tmp$minDMratio <- NA
        if( (tmp$minBMratio == -99)) tmp$minBMratio <- NA
        tmp$OccurringDrivers <- object$OccurringDrivers
        return(as.data.frame(tmp))
    }
}




print.oncosimul <- function(x, ...) {
    cat("\nIndividual OncoSimul trajectory with call:\n ")
    print(attributes(x)$call)
    cat("\n")
    print(summary(x))

    if(inherits(x, "oncosimul2")) {
        ## we know small object
        cat("\n")
        cat("Final population composition:\n")
        df <- data.frame(Genotype = x$GenotypesLabels,
                         N = x$pops.by.time[nrow(x$pops.by.time), -1])
        print(df)
    }
}

## I want this to return things that are storable
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
                              lwdClone = 0.9,
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
                           lwdClone = 0.9,
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

    ## uvx
    if(!inherits(x, "oncosimul2"))
        ndr <- colSums(x$Genotypes[1:x$NumDrivers, , drop = FALSE])
    else {
        ndr <- colSums(x$Genotypes[x$Drivers, , drop = FALSE])
    }
    
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
                     ndr,
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



plotPoset <- function(x, names = NULL, addroot = FALSE,
                      box = FALSE, ...) {
    if(is.null(names)) {
        if(addroot) names <- c("Root", 1:max(x))
        else names <- 1:max(x)
    }
    plot(posetToGraph(x, names =  names, addroot = addroot,
                      type = "graphNEL",
                      strictAdjMat = FALSE), ...)
    if(box)
        box()
}

plotAdjMat <- function(adjmat) {
    plot(as(adjmat, "graphNEL"))
}



############# The rest are internal functions

get.mut.vector.whole <- function(tmp, timeSample = "last", threshold = 0.5) {
    ## Obtain, from  results from a simulation run, the vector
    ## of 0/1 corresponding to each gene.
    
    ## threshold is the min. proportion for a mutation to be detected
    ## We are doing whole tumor sampling here, as in Sprouffske

    ## timeSample: do we sample at end, or at a time point, chosen
    ## randomly, from all those with at least one driver?
    
    
    if(timeSample == "last") {
        return(as.numeric(
            (tcrossprod(tmp$pops.by.time[nrow(tmp$pops.by.time), -1],
                        tmp$Genotypes)/tmp$TotalPopSize) > threshold))
    } else if (timeSample %in% c("uniform", "unif")) {
        the.time <- sample(which(tmp$PerSampleStats[, 4] > 0), 1)
        pop <- tmp$pops.by.time[the.time, -1]
        popSize <- tmp$PerSampleStats[the.time, 1]
        return( as.numeric((tcrossprod(pop,
                                       tmp$Genotypes)/popSize) > threshold) )
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


oncoSimul.internal <- function(poset, ## restrict.table,
                               numPassengers, 
                               ## numGenes,
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
                               seed,
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
                               ## endTimeEvery,
                               detectionDrivers,
                               onlyCancer,
                               errorHitWallTime,
                               max.num.tries,
                               errorHitMaxTries,
                               minDetectDrvCloneSz,
                               extraTime) {

    ## the value of 20000, in megabytes, for max.memory sets a limit of ~ 20 GB
  

    ## if(keepEvery < sampleEvery)
    ##     warning("setting keepEvery to sampleEvery")

    ## a backdoor to allow passing restrictionTables directly
    if(inherits(poset, "restrictionTable"))
        restrict.table <- poset
    else
        restrict.table <- poset.to.restrictTable(poset)
    numDrivers <- nrow(restrict.table)
    numGenes <- (numDrivers + numPassengers)
    
    if(numGenes > 64)
        stop("Largest possible number of genes is 64")

    
    if(initSize_species < 10) {
        warning("initSize_species too small?")
    }
    if(initSize_iter < 100) {
        warning("initSize_iter too small?")
    }

    ## numDrivers <- nrow(restrict.table)
    if(length(unique(restrict.table[, 1])) != numDrivers)
        stop("BAIL OUT NOW: length(unique(restrict.table[, 1])) != numDrivers)")
    ddr <- restrict.table[, 1]
    if(any(diff(ddr) != 1))
        stop("BAIL OUT NOW:  any(diff(ddr) != 1")
    ## sanity checks
    if(max(restrict.table[, 1]) != numDrivers)
        stop("BAIL OUT NOW: max(restrict.table[, 1]) != numDrivers")
    if(numDrivers > numGenes)
        stop("BAIL OUT NOW: numDrivers > numGenes")
    
    non.dep.drivers <- restrict.table[which(restrict.table[, 2] == 0), 1]




    ## if( (is.null(endTimeEvery) || (endTimeEvery > 0)) &&
    ##    (typeFitness %in% c("bozic1", "exp") )) {
    ##     warning(paste("endTimeEvery will take a positive value. ",
    ##                   "This will make simulations not stop until the next ",
    ##                   "endTimeEvery has been reached. Thus, in simulations ",
    ##                   "with very fast growth, simulations can take a long ",
    ##                   "time to finish, or can hit the wall time limit. "))
    ## }
    ## if(is.null(endTimeEvery))
    ##     endTimeEvery <- keepEvery
    ## if( (endTimeEvery > 0) && (endTimeEvery %% keepEvery) )
    ##     warning("!(endTimeEvery %% keepEvery)")
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

    return(c(
        BNB_Algo5(rtC,
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
        seed,
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
        # endTimeEvery,
        detectionDrivers,
        onlyCancer,
        errorHitWallTime,
        max.num.tries,
        errorHitMaxTries,
        minDetectDrvCloneSz,
        extraTime
    ),
             NumDrivers = numDrivers
             ))

}






create.muts.by.time <- function(tmp) { ## tmp is the output from Algorithm5
    if(tmp$NumClones > 1) {
        NumMutations <- apply(tmp$Genotypes, 2, sum)
        muts.by.time <- cbind(tmp$pops.by.time[, c(1), drop = FALSE],
                              t(apply(tmp$pops.by.time[, -c(1),
                                                       drop = FALSE], 1,
                                      function(x) tapply(x,
                                                         NumMutations, sum))))
        colnames(muts.by.time)[c(1)] <- "Time"
    } else {
        muts.by.time <- tmp$pops.by.time
    }
    return(muts.by.time)
} 


create.drivers.by.time <- function(tmp, ndr) {
    ## CountNumDrivers <- apply(tmp$Genotypes[1:numDrivers, ,drop = FALSE], 2, sum)
    CountNumDrivers <- ndr
    if(tmp$NumClones >= 1) {
        if(tmp$NumClones == 1) {
            if(ncol(tmp$pops.by.time) != 2)
                stop("This is impossible!")
            drivers.by.time <- tmp$pops.by.time
        } else {
            if(length(unique(CountNumDrivers )) > 1) {
                drivers.by.time <- cbind(tmp$pops.by.time[, c(1),
                                                          drop = FALSE] ,
                                         t(apply(tmp$pops.by.time[, -c(1),
                                                                  drop = FALSE],
                                                 1,
                                                 function(x)
                                                 tapply(x,
                                                        CountNumDrivers,
                                                        sum)))) 
            } else {
                drivers.by.time <- cbind(tmp$pops.by.time[, c(1),
                                                          drop = FALSE] ,
                                         rowSums(tmp$pops.by.time[, -c(1),
                                                                  drop =FALSE]))
            }
        }
        colnames(drivers.by.time) <- c("Time",
                                       paste("dr_",
                                             colnames(drivers.by.time)[-c(1)],
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
                         ndr,
                         timescale = 4,
                         trim.no.drivers = TRUE,
                         addtot = TRUE,
                         addtotlwd = 2,
                         na.subs = TRUE, log = "y", type = "l",
                         lty = 1:9, col = c(8, "orange", 6:1),
                         lwd = 2,
                         ...) {
    ## z <- create.drivers.by.time(x, numDrivers)
    z <- create.drivers.by.time(x, ndr)
    if(trim.no.drivers && x$MaxDriversLast) {
        fi <- which(apply(z[, -c(1, 2), drop = FALSE], 1,
                          function(x) sum(x) > 0))[1]
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


## simulate from generative model. This might not be fully correct!!!

simposet <- function(poset, p) {
    ## if (length(parent.nodes) != length (child.nodes)){
    ##     print("An Error Occurred")
    ## }
    ##    else {
    num.genes <- max(poset) - 1 ## as root is not a gene
    genotype <-t(c(1, rep(NA, num.genes)))
    colnames(genotype) <- as.character(0:num.genes)
    
    
    poset$runif <- runif(nrow(poset))
    ## this.relation.prob.OK could be done outside, but having it inside
    ## the loop would allow to use different thresholds for different
    ## relationships
    for (i in (1:nrow(poset))) {
        child <- poset[i, 2]
        this.relation.prob.OK <- as.numeric(poset[i, "runif"] > p)
        the.parent <- genotype[ poset[i, 1] ] ## it's the value of parent in genotype. 
        if (is.na(genotype[child])){
            genotype[child] <- this.relation.prob.OK * the.parent  
        }
        else
            genotype[child] <- genotype[child]*(this.relation.prob.OK * the.parent)
    }
    ##    }
    
    return(genotype)
}




eFinalMf <- function(initSize, s, j) {
    ## Expected final sizes for McF, when K is set to the default.
    # j is number of drivers
    ## as it says, with no passengers
    ## Set B(d) = D(N)
    K <- initSize/(exp(1) - 1)
    return(K * (exp( (1 + s)^j) - 1))
}


## to plot and adjacency matrix in this context can do
## plotPoset(intAdjMatToPoset(adjMat))
## where intAdjMatToPoset is from best oncotree code: generate-random-trees.
## No! the above is simpler
