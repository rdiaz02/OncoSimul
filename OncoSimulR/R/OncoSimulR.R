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


## Does it even make sense to keepPhylog with oncoSimulSample?
## Yes, but think what we store.


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
                                if(length(nd) == 1) ## for sample
                                       nd <- c(nd, nd)
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
                            keepPhylog = FALSE, 
                            max.memory = 2000,
                            max.wall.time.total = 600,
                            max.num.tries.total = 500 * Nindiv,
                            ## well, obviously they are errors
                            ## errorHitWallTime = TRUE,
                            ## errorHitMaxTries = TRUE,
                            typeSample = "whole",
                            thresholdWhole = 0.5,
                            initMutant = NULL,
                            verbosity  = 1,
                            seed = "auto"){
    ## No longer using mclapply, because of the way we use the limit on
    ## the number of tries.
    
    ## leaving detectionSize and detectionDrivers as they are, produces
    ## the equivalente of uniform sampling. For last, fix a single number

    if(keepPhylog)
        warning(paste("oncoSimulSample does not return the phylogeny",
                      "for now, so there is little point in storing it."))

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

    ## FIXME: we are not triggering an error, just a message. This is on
    ## purpose, since some of these conditions DO provide useful
    ## output. Do we want these to be errors?
    f.out.attempts <- function() {
        message("Run out of attempts")
        return(list(
            popSummary = NA,
            popSample = NA,
            HittedMaxTries = TRUE,
            HittedWallTime = FALSE,
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
            HittedWallTime = FALSE,
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
                               errorHitMaxTries = TRUE,
                               seed = seed,
                               initMutant = initMutant,
                               keepPhylog = keepPhylog)
        
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
                message(paste0("Successfully sampled ", Nindiv, " individuals"))
            class(pop) <- "oncosimulpop"
            if(inherits(fp, "fitnessEffects")) {
                geneNames <- names(getNamesID(fp))
            } else {
                geneNames <- NULL
            }
            return(list(
                popSummary = summary(pop),
                popSample = samplePop(pop, timeSample = "last",
                                      typeSample = typeSample,
                                      thresholdWhole = thresholdWhole,
                                      geneNames = geneNames),
                attemptsUsed = attemptsUsed,
                probCancer = Nindiv/attemptsUsed,
                HittedMaxTries = FALSE,
                HittedWallTime = FALSE,
                UnrecoverExcept = FALSE
            ))
        } else if( attemptsLeft <= 0 ) {
              ## it is very unlikely this will ever happen. 
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
                         keepPhylog = FALSE,
                         max.memory = 2000,
                         max.wall.time = 200,
                         max.num.tries = 500,
                         errorHitWallTime = TRUE,
                         errorHitMaxTries = TRUE,
                         initMutant = NULL,
                         verbosity  = 0,
                         mc.cores = detectCores(),
                         seed = "auto") {

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
                        verbosity = verbosity,
                        seed = seed, keepPhylog = keepPhylog,
                        initMutant = initMutant
                    ),
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
                           keepPhylog = FALSE,
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
    ## We do not test for equality to 0. That might be a weird but
    ## legitimate case?

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
           
        }
        if(length(initMutant) > 1)
            stop("With the old poset, initMutant can only take a single value.")
        ## Seeding C++ is now much better in new version
        if(is.null(seed) || (seed == "auto")) {## passing a null creates a random seed
            ## name is a legacy. This is really the seed for the C++ generator.
            ## Nope, we cannot use 2^32, because as.integer will fail.
            seed <- as.integer(round(runif(1, min = 0, max = 2^16)))
        }
        if(verbosity >= 2)
            cat(paste("\n Using ", seed, " as seed for C++ generator\n"))

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
        if(is.null(seed)) {## Passing a null creates a random seed
            ## We use a double, to be able to pass in range > 2^16.
            ## Do not use 0, as that is our way of signaling to C++ to
            ## generate the seed.
            seed <- round(runif(1, min = 1, max = 2^32))
            if(verbosity >= 2)
                cat(paste("\n Using ", seed, " as seed for C++ generator\n"))
        } else if(seed == "auto") {
            seed <- 0.0
            if(verbosity >= 2)
                cat("\n A (high quality) random seed will be generated in C++\n")
        }
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
                                        errorHitMaxTries = errorHitMaxTries,
                                        keepPhylog = keepPhylog),
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
    
    if(object$other$UnrecoverExcept) ## yes, when bailing out from
                                     ## except. can have just minimal
                                     ## content
        return(NA)
    else if (object$HittedWallTime || object$HittedMaxTries)
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
    cat("\nPopulation of OncoSimul trajectories of",
        length(x), "individuals. Call :\n")
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
                              plotDiversity = FALSE,
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
                                  plotDiversity = plotDiversity,
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
                           plotDiversity = FALSE,
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
    if(plotDiversity) {
        par(fig = c(0, 1, 0.8, 1))
        m1 <- par()$mar
        m <- m1
        m[c(1, 3)] <- c(0, 0.7)
        op <- par(mar = m )
        plotShannon(x)
        par(op)
        m1[c(3)] <- 0.2
        op <- par(mar = m1)
        par(fig = c(0, 1, 0, 0.8), new = TRUE)  
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

## this function seems to never be used
## plotAdjMat <- function(adjmat) {
##     plot(as(adjmat, "graphNEL"))
## }



which_N_at_T <- function(x, N = 1, t = "last") {
    if((length(t) == 1) && (t == "last"))
        T <- nrow(x$pops.by.time)
    else if(length(t) == 2) {
        if(t[1] < 0)
            warning("smallest t must be >= 0; setting to 0")
        if(t[1] > t[2])
            stop("t[1] must be <= t[2]")
        if(t[2] > max(x$pops.by.time[, 1]))
            message("t[2] > largest time; setting it to the max")
        T <- which(
            (x$pops.by.time[, 1] >= t[1]) &
                (x$pops.by.time[, 1] <= t[2]))
    }
    else
        stop("t must be either 'last' or a vector of length 2")
    z <- which(x$pops.by.time[T, -1, drop = FALSE] >= N, arr.ind = TRUE)[, 2]
    z <- unique(z) ## recall we removed first column but we index from first.
    return(z)
}

phylogClone <- function(x, N = 1, t = "last", keepEvents = TRUE) {
    if(!inherits(x, "oncosimul2"))
        stop("Phylogenetic information is only stored with v >= 2")
    z <- which_N_at_T(x, N, t)
    tG <- x$GenotypesLabels[z] ## only for GenotypesLabels we keep all
                               ## sample size info at each period
    df <- x$other$PhylogDF
    if(!keepEvents) { ## is this just a graphical thing? or not?
        df <- df[!duplicated(df[, c(1, 2)]), ]
    }
    g <- igraph::graph.data.frame(df[, c(1, 2)])
    ## nodes <- match(tG, V(g)$name)
    nodesInP <- unique(unlist(igraph::neighborhood(g, order = 1e9,
                                                   nodes = tG,
                                                   mode = "in")))
    ## Remember that the phylog info can contain clones that are
    ## not in pops.by.time, as they go extinct between creation
    ## and sampling.
    allLabels <- unique(as.character(unlist(df[, c(1, 2)])))
    nodesRm <- setdiff(allLabels, V(g)$name[nodesInP])
    g <- igraph::delete.vertices(g, nodesRm)
    tmp <- list(graph = g, df = df)
    class(tmp) <- c(class(tmp), "phylogClone")
    return(tmp)
    ## trivial to return an adjacency matrix if needed. The keepEvents = FALSE
}



plotClonePhylog <- function(x, N = 1, t = "last",
                            timeEvents = FALSE,
                            keepEvents = FALSE,
                            fixOverlap = TRUE,
                            returnGraph = FALSE, ...) {
    if(!inherits(x, "oncosimul2"))
        stop("Phylogenetic information is only stored with v >=2")
    if(nrow(x$other$PhylogDF) == 0)
        stop("It seems you run the simulation with keepPhylog= FALSE")
    pc <- phylogClone(x, N, t, keepEvents)
    l0 <- igraph::layout.reingold.tilford(pc$g)
    if(!timeEvents) {
        plot(pc$g, layout = l0)
    } else {
        l1 <- l0
        indexAppear <- match(V(pc$g)$name, as.character(pc$df[, 2]))
        firstAppear <- pc$df$time[indexAppear]
        firstAppear[1] <- 0
        l1[, 2] <- (max(firstAppear) - firstAppear)
        if(fixOverlap) {
            dx <- which(duplicated(l1[, 1]))
            if(length(dx)) {
                ra <- range(l1[, 1])
                l1[dx, 1] <- runif(length(dx), ra[1], ra[2])
            }
        }
        plot(pc$g, layout = l1)         
    }
    if(returnGraph)
        return(pc$g)
}








############# The rest are internal functions


get.the.time.for.sample <- function(tmp, timeSample) {
    if(timeSample == "last") {
        if(tmp$TotalPopSize == 0) {
            warning(paste("Final population size is 0.",
                          "Thus, there is nothing to sample with",
                          "sampling last. You will get NAs"))
            the.time <- -99
        } else {
              the.time <- nrow(tmp$pops.by.time)
          }
    } else if (timeSample %in% c("uniform", "unif")) {
          candidate.time <- which(tmp$PerSampleStats[, 4] > 0)
          
          if (length(candidate.time) == 0) {
              warning(paste("There is not a single sampled time",
                            "at which there are any mutants.",
                            "Thus, no uniform sampling possible.",
                            "You will get NAs"))
              the.time <- -99
              ## return(rep(NA, nrow(tmp$Genotypes)))
          } else if (length(candidate.time) == 1) {
                message("Only one sampled time period with mutants.")
                the.time <- candidate.time
            } else {
                  the.time <- sample(candidate.time, 1)
              }
      } else {
            stop("Unknown timeSample option")
        }
    return(the.time)
}


get.mut.vector <- function(x, timeSample, typeSample,
                           thresholdWhole) {
    the.time <- get.the.time.for.sample(x, timeSample)
    if(the.time < 0) { 
        return(rep(NA, nrow(x$Genotypes)))
    } 
    pop <- x$pops.by.time[the.time, -1]
    
    if(all(pop == 0)) {
        stop("You found a bug: this should never happen")
    }
    
    if(typeSample %in% c("wholeTumor", "whole")) {
        popSize <- x$PerSampleStats[the.time, 1]
        return( as.numeric((tcrossprod(pop,
                                       x$Genotypes)/popSize) >= thresholdWhole) )
    } else if (typeSample %in%  c("singleCell", "single")) {

          return(x$Genotypes[, sample(seq_along(pop), 1, prob = pop)])
      } else {
            stop("Unknown typeSample option")
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

    ## These can never be set by the user
    ## if(initSize_species < 10) {
    ##     warning("initSize_species too small?")
    ## }
    ## if(initSize_iter < 100) {
    ##     warning("initSize_iter too small?")
    ## }

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

eFinalMf <- function(initSize, s, j) {
    ## Expected final sizes for McF, when K is set to the default.
    # j is number of drivers
    ## as it says, with no passengers
    ## Set B(d) = D(N)
    K <- initSize/(exp(1) - 1)
    return(K * (exp( (1 + s)^j) - 1))
}



## We are not using this anymore
## create.muts.by.time <- function(tmp) { ## tmp is the output from Algorithm5
##     if(tmp$NumClones > 1) {
##         NumMutations <- apply(tmp$Genotypes, 2, sum)
##         muts.by.time <- cbind(tmp$pops.by.time[, c(1), drop = FALSE],
##                               t(apply(tmp$pops.by.time[, -c(1),
##                                                        drop = FALSE], 1,
##                                       function(x) tapply(x,
##                                                          NumMutations, sum))))
##         colnames(muts.by.time)[c(1)] <- "Time"
##     } else {
##         muts.by.time <- tmp$pops.by.time
##     }
##     return(muts.by.time)
## } 


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

shannonI <- function(x) {
    sx <- sum(x)
    p <- x/sx
    p <- p[p > 0]
    return(-sum(p * log(p)))
}

plotShannon <- function(z) {
    h <- apply(z$pops.by.time[, 2:ncol(z$pops.by.time), drop = FALSE],
               1, shannonI)
    plot(x = z$pops.by.time[, 1],
         y = h, type = "l", xlab = "", ylab = "H", axes = FALSE)
    box()
    axis(2)
}

## simpsonI <- function(x) {
##     sx <- sum(x)
##     p <- x/sx
##     p <- p[p > 0]
##     return(sum(p^2)))
## }

## plotSimpson <- function(z) {
    
##     h <- apply(z$pops.by.time[, 2:ncol(z$pops.by.time), drop = FALSE],
##                1, shannonI)
##     plot(x = z$pops.by.time[, 1],
##          y = h, lty = "l", xlab = "", ylab = "H")
## }


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
    box()
}


plotDrivers0 <- function(x,
                         ndr,
                         timescale = 1,
                         trim.no.drivers = FALSE,
                         addtot = TRUE,
                         addtotlwd = 2,
                         na.subs = TRUE, log = "y", type = "l",
                         lty = 1:9, col = c(8, "orange", 6:1),
                         lwd = 2,
                         ...) {
    ## z <- create.drivers.by.time(x, numDrivers)
    z <- create.drivers.by.time(x, ndr)
    ## we can now never enter here because trim.no.drivers is always FALSE
    ## in call.
    if(trim.no.drivers && x$MaxDriversLast) {
        fi <- which(apply(z[, -c(1, 2), drop = FALSE], 1,
                          function(x) sum(x) > 0))[1]
        z <- z[fi:nrow(z), , drop = FALSE]
    }
    y <- z[, 2:ncol(z), drop = FALSE]
    if(na.subs){
        y[y == 0] <- NA
    }

    ## Likewise, we can never enter here now as timescale fixed at 1. And
    ## this is silly.
    ## if(timescale != 1) {
    ##     time <- timescale * z[, 1]
    ## } else {
    ##     time <- z[, 1]
    ## }
    time <- timescale * z[, 1]
    
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




## No longer used
## rtNoDep <- function(numdrivers) {
##     ## create a restriction table with no dependencies
##     x <- matrix(nrow = numdrivers, ncol = 3)
##     x[, 1] <- 1:numdrivers
##     x[, 2] <- 0
##     x[, 3] <- -9
##     return(x)
## }


## Simulate from generative model. This is a quick function, and is most
## likely wrong! Never used for anything.

## simposet <- function(poset, p) {
##     ## if (length(parent.nodes) != length (child.nodes)){
##     ##     print("An Error Occurred")
##     ## }
##     ##    else {
##     num.genes <- max(poset) - 1 ## as root is not a gene
##     genotype <-t(c(1, rep(NA, num.genes)))
##     colnames(genotype) <- as.character(0:num.genes)
    
    
##     poset$runif <- runif(nrow(poset))
##     ## this.relation.prob.OK could be done outside, but having it inside
##     ## the loop would allow to use different thresholds for different
##     ## relationships
##     for (i in (1:nrow(poset))) {
##         child <- poset[i, 2]
##         this.relation.prob.OK <- as.numeric(poset[i, "runif"] > p)
##         the.parent <- genotype[ poset[i, 1] ] ## it's the value of parent in genotype. 
##         if (is.na(genotype[child])){
##             genotype[child] <- this.relation.prob.OK * the.parent  
##         }
##         else
##             genotype[child] <- genotype[child]*(this.relation.prob.OK * the.parent)
##     }
##     ##    }
    
##     return(genotype)
## }


## to plot and adjacency matrix in this context can do
## plotPoset(intAdjMatToPoset(adjMat))
## where intAdjMatToPoset is from best oncotree code: generate-random-trees.
## No! the above is simpler




## get.mut.vector.whole <- function(tmp, timeSample = "last", threshold = 0.5) {
##     ## Obtain, from  results from a simulation run, the vector
##     ## of 0/1 corresponding to each gene.
    
##     ## threshold is the min. proportion for a mutation to be detected
##     ## We are doing whole tumor sampling here, as in Sprouffske

##     ## timeSample: do we sample at end, or at a time point, chosen
##     ## randomly, from all those with at least one driver?
    
##     if(timeSample == "last") {
##         if(tmp$TotalPopSize == 0)
##             warning(paste("Final population size is 0.",
##                           "Thus, there is nothing to sample with ",
##                           "sampling last. You will get NAs"))
##         return(as.numeric(
##             (tcrossprod(tmp$pops.by.time[nrow(tmp$pops.by.time), -1],
##                         tmp$Genotypes)/tmp$TotalPopSize) > threshold))
##     } else if (timeSample %in% c("uniform", "unif")) {
##           candidate.time <- which(tmp$PerSampleStats[, 4] > 0)
          
##           if (length(candidate.time) == 0) {
##               warning(paste("There is not a single sampled time",
##                             "at which there are any mutants.",
##                             "Thus, no uniform sampling possible.",
##                             "You will get NAs"))
##               return(rep(NA, nrow(tmp$Genotypes)))
##           } else if (length(candidate.time) == 1) {
##                 the.time <- candidate.time
##             } else {
##                   the.time <- sample(candidate.time, 1)
##               }
##           pop <- tmp$pops.by.time[the.time, -1]
##           popSize <- tmp$PerSampleStats[the.time, 1]
##           ## if(popSize == 0)
##           ##     warning(paste("Population size at this time is 0.",
##           ##                   "Thus, there is nothing to sample at this time point.",
##           ##                   "You will get NAs"))
##           return( as.numeric((tcrossprod(pop,
##                                        tmp$Genotypes)/popSize) > threshold) )
##       }
## }



##           the.time <- sample(which(tmp$PerSampleStats[, 4] > 0), 1)
##           if(length(the.time) == 0) {
##               warning(paste("There are no clones with drivers at any time point.",
##                             "No uniform sampling possible.",
##                             "You will get a vector of NAs."))
##             return(rep(NA, nrow(tmp$Genotypes)))  
##           }
## get.mut.vector.singlecell <- function(tmp, timeSample = "last") {
##     ## No threshold, as single cell.

##     ## timeSample: do we sample at end, or at a time point, chosen
##     ## randomly, from all those with at least one driver?
    
##     if(timeSample == "last") {
##         the.time <- nrow(tmp$pops.by.time)
##     } else if (timeSample %in% c("uniform", "unif")) {
##          candidate.time <- which(tmp$PerSampleStats[, 4] > 0)
         
##          if (length(candidate.time) == 0) {
##              warning(paste("There is not a single sampled time",
##                            "at which there are any mutants.",
##                            "Thus, no uniform sampling possible.",
##                            "You will get NAs"))
##              return(rep(NA, nrow(tmp$Genotypes)))
##          } else if (length(candidate.time) == 1) {
##                the.time <- candidate.time
##            } else {
##                  the.time <- sample(candidate.time, 1)
##              }

##      }
##     pop <- tmp$pops.by.time[the.time, -1]
##     ##       popSize <- tmp$PerSampleStats[the.time, 1]
##     ## genot <- sample(seq_along(pop), 1, prob = pop)
##     if(all(pop == 0)) {
##         warning(paste("All clones have a population size of 0",
##                       "at the chosen time. Nothing to sample.",
##                       "You will get NAs"))
##         return(rep(NA, nrow(tmp$Genotypes)))
##     } else {
##           return(tmp$Genotypes[, sample(seq_along(pop), 1, prob = pop)])
##       }
## }


## get.mut.vector <- function(x, timeSample = "whole", typeSample = "last",
##                            thresholdWhole = 0.5) {
##     if(typeSample %in% c("wholeTumor", "whole")) {
##         get.mut.vector.whole(x, timeSample = timeSample,
##                              threshold = thresholdWhole)
##     } else if(typeSample %in%  c("singleCell", "single")) {
##         get.mut.vector.singlecell(x, timeSample = timeSample)
##     }
## }





## plotClonePhylog <- function(x, timeEvent = FALSE,
##                             showEvents = TRUE,
##                             fixOverlap = TRUE) {
##     if(!inherits(x, "oncosimul2"))
##         stop("Phylogenetic information is only stored with v >=2")
##     if(nrow(x$other$PhylogDF) == 0)
##         stop("It seems you run the simulation with keepPhylog= FALSE")
##     ## requireNamespace("igraph")
##     df <- x$other$PhylogDF
##     if(!showEvents) {
##         df <- df[!duplicated(df[, c(1, 2)]), ]
##     }
##     g <- igraph::graph.data.frame(df)
##     l0 <- igraph::layout.reingold.tilford(g)
##     if(!timeEvent) {
##         plot(g, layout = l0)
##     } else {
##         l1 <- l0
##         indexAppear <- match(V(g)$name, as.character(df[, 2]))
##         firstAppear <- df$time[indexAppear]
##         firstAppear[1] <- 0
##         l1[, 2] <- (max(firstAppear) - firstAppear)
##         if(fixOverlap) {
##             dx <- which(duplicated(l1[, 1]))
##             if(length(dx)) {
##                 ra <- range(l1[, 1])
##                 l1[dx, 1] <- runif(length(dx), ra[1], ra[2])
##             }
##         }
##         plot(g, layout = l1)         
##     }
## }
