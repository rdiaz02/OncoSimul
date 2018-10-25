## Copyright 2013, 2014, 2015, 2016, 2017 Ramon Diaz-Uriarte

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
                            muEF = NULL,
                            detectionSize = round(runif(Nindiv, 1e5, 1e8)),
                            detectionDrivers = {
                                if(inherits(fp, "fitnessEffects")) {
                                    if(length(fp$drv)) {
                                        nd <- (2: round(0.75 * length(fp$drv)))
                                    } else {
                                        nd <- 9e6 
                                    }
                                } else {
                                    nd <- (2 : round(0.75 * max(fp)))
                                }
                                if(length(nd) == 1) ## for sample
                                       nd <- c(nd, nd)
                                sample(nd, Nindiv,
                                       replace = TRUE)
                            },
                            detectionProb = NA,
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
                            mutationPropGrowth = ifelse(model == "Bozic",
                                                        FALSE, TRUE),
                            max.memory = 2000,
                            max.wall.time.total = 600,
                            max.num.tries.total = 500 * Nindiv,
                            ## well, obviously they are errors
                            ## errorHitWallTime = TRUE,
                            ## errorHitMaxTries = TRUE,
                            typeSample = "whole",
                            thresholdWhole = 0.5,
                            initMutant = NULL,
                            AND_DrvProbExit = FALSE,
                            fixation = NULL,
                            verbosity  = 1,
                            showProgress = FALSE,
                            seed = "auto"){
    ## No longer using mclapply, because of the way we use the limit on
    ## the number of tries.
    
    ## leaving detectionSize and detectionDrivers as they are, produces
    ## the equivalente of uniform sampling. For last, fix a single number

    ## detectionDrivers when there are none: had we left it at 0, then
    ## when there are no drivers we would stop at the first sampling
    ## period.
    
    if(Nindiv < 1)
        stop("Nindiv must be >= 1")
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
                               muEF = muEF,
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
                               keepPhylog = keepPhylog,
                               mutationPropGrowth = mutationPropGrowth,
                               detectionProb = detectionProb,
                               AND_DrvProbExit = AND_DrvProbExit,
                               fixation = fixation)        
        if(tmp$other$UnrecoverExcept) {
            return(f.out.unrecover.except(tmp))
        }
        
        pop[[indiv]] <- tmp
        numToRun <- (numToRun - 1)
        attemptsUsed <- attemptsUsed + tmp$other$attemptsUsed
        attemptsLeft <- (max.num.tries.total - attemptsUsed)
        if(showProgress) {
            cat("......  Done individual ", indiv,
                ". Used ", attemptsUsed, "attempts.",
                ". Running for ", as.double(difftime(Sys.time(),
                                                    startTime, units = "secs")),
                " secs.\n"
                )
        }
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

add_noise <- function(x, properr) {
    if(properr <= 0) {
        return(x)
    }
    else {
        if(properr > 1)
            stop("Proportion with error cannot be > 1")
        nn <- prod(dim(x))
        flipped <- sample(nn, round(nn * properr))
        x[flipped] <- as.integer(!x[flipped])
        return(x)
    }
}

samplePop <- function(x, timeSample = "last",
                      typeSample = "whole",
                      thresholdWhole = 0.5,
                      geneNames = NULL,
                      popSizeSample = NULL,
                      propError = 0) {
    ## timeSample <- match.arg(timeSample)
    gN <- geneNames
    
    if(!is.null(popSizeSample) && (length(popSizeSample) > 1) &&
       (length(popSizeSample) != length(x))) {
        message("length popSizeSample != number of subjects")
        popSizeSample <- rep(popSizeSample, length.out = length(x))
        }
    ## A hack to prevent Map from crashing with a NULL
    if(is.null(popSizeSample)) popSizeSample <- -99
    if(inherits(x, "oncosimulpop")) {
        z <- do.call(rbind,
                     Map(get.mut.vector,
                         x = x,
                         timeSample = timeSample,
                         typeSample = typeSample,
                         thresholdWhole = thresholdWhole,
                         popSizeSample = popSizeSample
                         ))
        ## We need to check if the object is coming from v.2., to avoid
        ## having to force passing a vector of names
        if(is.null(gN) && (!is.null(x[[1]]$geneNames)))
            gN <- x[[1]]$geneNames
    } else {
        z <- get.mut.vector(x,
                            timeSample = timeSample,
                            typeSample = typeSample,
                            thresholdWhole = thresholdWhole,
                            popSizeSample = popSizeSample)
        dim(z) <- c(1, length(z))
        if(is.null(gN) && (!is.null(x$geneNames)))
            gN <- x$geneNames
    }
    message("\n Subjects by Genes matrix of ",
        nrow(z), " subjects and ",
            ncol(z), " genes.\n")

    if(propError > 0) {
        z <- add_noise(z, propError)
    }
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
                         numPassengers = 0,
                         mu = 1e-6,
                         muEF = NULL,
                         detectionSize = 1e8,
                         detectionDrivers = 4,
                         detectionProb = NA,
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
                         mutationPropGrowth = ifelse(model == "Bozic",
                                                     FALSE, TRUE),
                         max.memory = 2000,
                         max.wall.time = 200,
                         max.num.tries = 500,
                         errorHitWallTime = TRUE,
                         errorHitMaxTries = TRUE,
                         initMutant = NULL,
                         AND_DrvProbExit = FALSE,
                         fixation = NULL,
                         verbosity  = 0,
                         mc.cores = detectCores(),
                         seed = "auto") {

    if(Nindiv < 1)
        stop("Nindiv must be >= 1")
    
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
                        muEF = muEF,
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
                        initMutant = initMutant,
                        mutationPropGrowth = mutationPropGrowth,
                        detectionProb = detectionProb,
                        AND_DrvProbExit = AND_DrvProbExit,
                        fixation = fixation),
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
                           numPassengers = 0,
                           mu = 1e-6,
                           muEF = NULL,
                           detectionSize = 1e8,
                           detectionDrivers = 4,
                           detectionProb = NA,
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
                           mutationPropGrowth = ifelse(model == "Bozic",
                                                       FALSE, TRUE),
                           max.memory = 2000,
                           max.wall.time = 200,
                           max.num.tries = 500,
                           errorHitWallTime = TRUE,
                           errorHitMaxTries = TRUE,
                           verbosity = 0,
                           initMutant = NULL,
                           AND_DrvProbExit = FALSE,
                           fixation = NULL,
                           seed = NULL
                           ) {
    call <- match.call()
    if(all(c(is_null_na(detectionProb),
             is_null_na(detectionSize),
             is_null_na(detectionDrivers),
             is_null_na(finalTime),
             is_null_na(fixation)
             )))
        stop("At least one stopping condition should be given.",
             " At least one of detectionProb, detectionSize, detectionDrivers,",
             " finalTime. Otherwise, we'll run until aborted by max.wall.time,",
             " max.num.tries, and the like.")

    if(AND_DrvProbExit && (is_null_na(detectionProb) || is_null_na(detectionDrivers)))
        stop("AND_DrvProbExit is TRUE: both of detectionProb and detectionDrivers",
             " must be non NA.")
    if(AND_DrvProbExit && !is_null_na(detectionSize)) {
        warning("With AND_DrvProbExit = TRUE, detectionSize is ignored.")
        detectionSize <- NA
    }
    if(inherits(fp, "fitnessEffects")) {
        s <- sh <- NULL ## force it soon!
    }

    ## legacies from poor name choices
    typeFitness <- switch(model,
                          "Bozic" = "bozic1",
                          "Exp" = "exp",
                          "McFarlandLog" = "mcfarlandlog",
                          "McFL" = "mcfarlandlog",
                          stop("No valid value for model")
                          )
    if(initSize < 1)
        stop("initSize < 1")
    
    if( (K < 1) && (model %in% c("McFL", "McFarlandLog") )) {
        stop("Using McFarland's model: K cannot be < 1")
    }       ##  if ( !(model %in% c("McFL", "McFarlandLog") )) {
            ## K <- 1 ## K is ONLY used for McFarland; set it to 1, to avoid
            ##        ## C++ blowing.

    if(typeFitness == "exp") {
        death <- 1
        ## mutationPropGrowth <- 1
    } else {
        death <- -99
        ## mutationPropGrowth <- 0
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
            minDetectDrvCloneSz <- initSize
        ## minDetectDrvCloneSz <- eFinalMf(initSize, s, detectionDrivers)
        else
            stop("Unknown model")
    }

    if( (length(mu) > 1) && !inherits(fp, "fitnessEffects"))
        stop("Per-gene mutation rates cannot be used with the old poset format")

    if(any(mu < 0)) {
        stop("(at least one) mutation rate (mu) is negative")
    }
    ## We do not test for equality to 0. That might be a weird but
    ## legitimate case?

    ## No user-visible magic numbers
    ## if(is.null(keepEvery))
    ##     keepEvery <- -9
    if(is_null_na(keepEvery)) keepEvery <- -9

    
    if( (keepEvery > 0) & (keepEvery < sampleEvery)) {
        keepEvery <- sampleEvery
        warning("setting keepEvery <- sampleEvery")
    }

    if( (typeFitness == "bozic1") && (mutationPropGrowth) )
        warning("Using fitness Bozic (bozic1) with mutationPropGrowth = TRUE;",
                "this will have no effect.")

    if( (typeFitness == "exp") && (death != 1) )
        warning("Using fitness exp with death != 1")

    if(!is_null_na(detectionDrivers) && (detectionDrivers >= 1e9))
        stop("detectionDrivers > 1e9; this doesn't seem reasonable")
    if(is_null_na(detectionDrivers)) detectionDrivers <- (2^31) - 1
    if(is_null_na(detectionSize)) detectionSize <- Inf
    if(is_null_na(finalTime)) finalTime <- Inf

    if(is_null_na(sampleEvery)) stop("sampleEvery cannot be NULL or NA")
    
    if(!inherits(fp, "fitnessEffects")) {
        if(any(unlist(lapply(list(fp, 
                                  numPassengers,
                                  s, sh), is.null)))) {
            m <- paste("You are using the old poset format.",
                       "You must specify all of poset, numPassengers",
                       "s, and sh.")
            stop(m)
           
        }
        if(AND_DrvProbExit) {
            stop("The AND_DrvProbExit = TRUE setting is invalid",
                 " with the old poset format.")
        }
        if(!is.null(muEF))
            stop("Mutator effects cannot be specified with the old poset format.")
        if( length(initMutant) > 0)  
            warning("With the old poset format you can no longer use initMutant.",
                    " The initMutant you passed will be ignored.")
        ## stop("With the old poset, initMutant can only take a single value.")
        if(!is_null_na(fixation))
            stop("'fixation' cannot be specified with the old poset format.")
        ## Seeding C++ is now much better in new version
        if(is.null(seed) || (seed == "auto")) {## passing a null creates a random seed
            ## name is a legacy. This is really the seed for the C++ generator.
            ## Nope, we cannot use 2^32, because as.integer will fail.
            seed <- as.integer(round(runif(1, min = 0, max = 2^16)))
        }
        if(verbosity >= 2)
            cat(paste("\n Using ", seed, " as seed for C++ generator\n"))

        if(!is_null_na(detectionProb)) stop("detectionProb cannot be used in v.1 objects")
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
                                     speciesFS = 10000,  
                                     ratioForce = 2,
                                     typeFitness = typeFitness,
                                     max.memory = max.memory,
                                     mutationPropGrowth = mutationPropGrowth,                                   
                                     initMutant = -1, 
                                     max.wall.time = max.wall.time,
                                     max.num.tries = max.num.tries,
                                     keepEvery = keepEvery,  
                                     ## alpha = 0.0015,  
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
        s <- sh <- NULL ## force it.
        if(numPassengers != 0)
            warning(paste("Specifying numPassengers has no effect",
                          " when using fitnessEffects objects. ",
                          " The fitnessEffects objects are much more ",
                          "flexible and you can use, for example,",
                          "the noIntGenes component for passengers."))
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
        if(!is_null_na(fixation)) {
            if( (!is.list(fixation)) && (!is.vector(fixation))  )
                stop("'fixation' must be a list or a vector.")
            if(!(all(unlist(lapply(fixation, is.vector)))))
                stop("Each element of 'fixation' must be a single element character vector.")
            if(!(all(unlist(lapply(fixation, class)) == "character")))
                stop("Each element of 'fixation' must be a single element character vector.")
            if(!(all( unlist(lapply(fixation, length)) == 1)))
                stop("Each element of 'fixation' must be a single element character vector.")
            if(AND_DrvProbExit)
                stop("It makes no sense to pass AND_DrvProbExit and a fixation list.")
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
                                        speciesFS = 10000,  
                                        ratioForce = 2,
                                        typeFitness = typeFitness,
                                        max.memory = max.memory,
                                        mutationPropGrowth = mutationPropGrowth,                                   
                                        initMutant = initMutant, 
                                        max.wall.time = max.wall.time,
                                        max.num.tries = max.num.tries,
                                        keepEvery = keepEvery,  
                                        ## alpha = 0.0015,  
                                        K = K, 
                                        minDetectDrvCloneSz = minDetectDrvCloneSz,
                                        extraTime = extraTime,
                                        detectionDrivers = detectionDrivers,
                                        onlyCancer = onlyCancer,
                                        errorHitWallTime = errorHitWallTime,
                                        errorHitMaxTries = errorHitMaxTries,
                                        keepPhylog = keepPhylog,
                                        MMUEF = muEF,
                                        detectionProb = detectionProb,
                                        AND_DrvProbExit = AND_DrvProbExit,
                                        fixation = fixation),
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
    ## This should be present even in HittedWallTime and HittedMaxTries
    ## if those are not regarded as errors
    pbp <- ("pops.by.time" %in% names(object) )
    
    if(object$other$UnrecoverExcept) { ## yes, when bailing out from
                                     ## except. can have just minimal
                                     ## content
        return(NA)
    } else if( !pbp & (object$HittedWallTime || object$HittedMaxTries) ) {
        return(NA)
    } else if ( !pbp & !(object$HittedWallTime || object$HittedMaxTries) ) {
        stop("Eh? No pops.by.time but did not hit wall time or max tries? BUG!")
    } else {
        tmp <- object[c("NumClones", "TotalPopSize", "LargestClone",
                        "MaxNumDrivers", "MaxDriversLast",
                        "NumDriversLargestPop", "TotalPresentDrivers",
                        "FinalTime", "NumIter", "HittedWallTime",
                        "HittedMaxTries")]
 
        tmp$errorMF <- object$other$errorMF
        tmp$minDMratio <- object$other$minDMratio
        tmp$minBMratio <- object$other$minBMratio
        if( (tmp$errorMF == -99)) tmp$errorMF <- NA
        if( (tmp$minDMratio == -99)) tmp$minDMratio <- NA
        if( (tmp$minBMratio == -99)) tmp$minBMratio <- NA
        tmp$OccurringDrivers <- object$OccurringDrivers
        return(as.data.frame(tmp, stringsAsFactors = FALSE))
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
## summary.oncosimulpop <- function(object, ...) {
##     as.data.frame(rbindlist(lapply(object, summary)))
## }

summary.oncosimulpop <- function(object, ...) {
    tmp <- lapply(object, summary)
    rm <- which(unlist(lapply(tmp, function(x) (length(x) == 1) && (is.na(x)))))
    if(length(rm) > 0)
        if(length(rm) < length(object)) {
        warning("Some simulations seem to have failed and will be removed",
                " from the summary. The failed runs are ",
                paste(rm, collapse = ", "),
                ".")
        tmp <- tmp[-rm]
        } else {
            warning("All simulations failed.")
            return(NA)
        }
    as.data.frame(rbindlist(tmp))
}



print.oncosimulpop <- function(x, ...) {
    cat("\nPopulation of OncoSimul trajectories of",
        length(x), "individuals. Call :\n")
    print(attributes(x)$call)
    cat("\n")
    print(summary(x))
}


plot.oncosimulpop <- function(x, ask = TRUE,
                              show = "drivers", 
                              type = ifelse(show == "genotypes",
                                            "stacked", "line"),
                              col = "auto",
                              log = ifelse(type == "line", "y", ""),
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
                              ylim = NULL,
                              xlim = NULL,
                              thinData = FALSE,
                              thinData.keep = 0.1,
                              thinData.min = 2,
                              plotDiversity = FALSE,
                              order.method = "as.is",
                              stream.center = TRUE,
                              stream.frac.rand = 0.01,
                              stream.spar = 0.2,
                              border = NULL,
                              lwdStackedStream = 1,
                              srange = c(0.4, 1),
                              vrange = c(0.8, 1),
                              breakSortColors = "oe",
                              legend.ncols = "auto",
                              ...
                              ) {
    op <- par(ask = ask)
    on.exit(par(op))
    null <- lapply(x, function(z)
        plot.oncosimul(z,
                       show = show,
                       type = type,
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
                       ylim = ylim,
                       xlim = xlim,
                       thinData = thinData,
                       thinData.keep = thinData.keep,
                       thinData.min = thinData.min,
                       plotDiversity = plotDiversity,
                       order.method = order.method,
                       stream.center = stream.center,
                       stream.frac.rand = stream.frac.rand,
                       stream.spar = stream.spar,
                       border = border,
                       lwdStackedStream = lwdStackedStream,
                       srange = srange,
                       vrange = vrange,
                       breakSortColors = breakSortColors,
                       legend.ncols = legend.ncols,
                       ...))
}


## plot.oncosimul <- function(x, col = c(8, "orange", 6:1),
##                            log = "y",
##                            ltyClone = 2:6,
##                            lwdClone = 0.9,
##                            ltyDrivers = 1,
##                            lwdDrivers = 3,
##                            xlab = "Time units",
##                            ylab = "Number of cells",
##                            plotClones = TRUE,
##                            plotDrivers = TRUE,
##                            addtot = FALSE,
##                            addtotlwd = 0.5,
##                            yl = NULL,
##                            thinData = FALSE,
##                            thinData.keep = 0.1,
##                            thinData.min = 2,
##                            plotDiversity = FALSE,
##                            ...
##                            ) {

##     if(thinData)
##         x <- thin.pop.data(x, keep = thinData.keep, min.keep = thinData.min)

##     ## uvx
##     if(!inherits(x, "oncosimul2"))
##         ndr <- colSums(x$Genotypes[1:x$NumDrivers, , drop = FALSE])
##     else {
##         ndr <- colSums(x$Genotypes[x$Drivers, , drop = FALSE])
##     }
    
##     if(is.null(yl)) {
##         if(log %in% c("y", "xy", "yx") )
##             yl <- c(1, max(apply(x$pops.by.time[, -1, drop = FALSE], 1, sum)))
##         else
##             yl <- c(0, max(apply(x$pops.by.time[, -1, drop = FALSE], 1, sum)))
##     }
##     if(plotDiversity) {
##         par(fig = c(0, 1, 0.8, 1))
##         m1 <- par()$mar
##         m <- m1
##         m[c(1, 3)] <- c(0, 0.7)
##         op <- par(mar = m )
##         plotShannon(x)
##         par(op)
##         m1[c(3)] <- 0.2
##         op <- par(mar = m1)
##         par(fig = c(0, 1, 0, 0.8), new = TRUE)  
##     }
##     if(plotClones) {
##         plotClones(x,
##                    ndr = ndr, 
##                    xlab = xlab,
##                    ylab = ylab,
##                    lty = ltyClone,
##                    col = col, 
##                    ylim = yl,
##                    lwd = lwdClone,
##                    axes = FALSE,
##                    log = log,
##                    ...)
##     }

##     if(plotClones && plotDrivers)
##         par(new = TRUE)
    
##     if(plotDrivers){
##         plotDrivers0(x,
##                      ndr,
##                      timescale = 1,
##                      trim.no.drivers = FALSE,
##                      xlab = "", ylab = "",
##                      lwd = lwdDrivers,
##                      lty = ltyDrivers,
##                      col = col, 
##                      addtot = addtot,
##                      addtotlwd = addtotlwd,
##                      log = log, ylim = yl,
##                      ...)
##     }
    
## }


plot.oncosimul <- function(x,
                           show = "drivers", 
                           type = ifelse(show == "genotypes",
                                         "stacked", "line"),
                           col = "auto",
                           log = ifelse(type == "line", "y", ""),
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
                           ylim = NULL,
                           xlim = NULL,
                           thinData = FALSE,
                           thinData.keep = 0.1,
                           thinData.min = 2,
                           plotDiversity = FALSE,
                           order.method = "as.is",
                           stream.center = TRUE,
                           stream.frac.rand = 0.01,
                           stream.spar = 0.2,
                           border = NULL,
                           lwdStackedStream = 1,
                           srange = c(0.4, 1),
                           vrange = c(0.8, 1),
                           breakSortColors = "oe",
                           legend.ncols = "auto",
                           ...
                           ) {


    if(!(type %in% c("stacked", "stream", "line")))
        stop("Type of plot unknown: it must be one of",
             "stacked, stream or line")

    if(!(show %in% c("genotypes", "drivers")))
        stop("show must be one of ",
             "genotypes or drivers")

    if(!(breakSortColors %in% c("oe", "distave", "random")))
        stop("breakSortColors must be one of ",
             "oe, distave, or random")

    

    colauto <- FALSE
    if(col == "auto" && (type == "line") && (show == "drivers"))
        col <- c(8, "orange", 6:1)
    if(col == "auto" && (show == "genotypes")) {
        ## For categorical data, I find Dark2, Paired, or Set1 to work best.
        col <- colorRampPalette(brewer.pal(8, "Dark2"))(ncol(x$pops.by.time) - 1)
        colauto <- TRUE
    }
    
    if(show == "genotypes") {
        plotDrivers <- FALSE
        plotClones <- TRUE
    }
    
    if(thinData)
        x <- thin.pop.data(x, keep = thinData.keep, min.keep = thinData.min)

    if(!is.null(xlim))
        x <- xlim.pop.data(x, xlim)
    
    ## For genotypes, ndr is now the genotypes.  Actually, ndr is now just
    ## a sequence 1:(ncol(y) - 1)

    ## The user will want to change the colors, like a colorRamp, etc. Or
    ## rainbow.

    ## genotypes and line, always call plotDrivers0
    if(show == "drivers") {
        if(!inherits(x, "oncosimul2"))
            ndr <- colSums(x$Genotypes[1:x$NumDrivers, , drop = FALSE])
        else {
            ndr <- colSums(x$Genotypes[x$Drivers, , drop = FALSE])
        }
    } else { ## show we are showing genotypes
        ndr <- 1:(ncol(x$pops.by.time) - 1)
    }
    
    if((type == "line") && is.null(ylim)) {
        if(log %in% c("y", "xy", "yx") )
            ylim <- c(1, max(apply(x$pops.by.time[, -1, drop = FALSE], 1, sum)))
        else
            ylim <- c(0, max(apply(x$pops.by.time[, -1, drop = FALSE], 1, sum)))
    }
    if(plotDiversity) {
        oppd <- par(fig = c(0, 1, 0.8, 1))
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

    ## Shows its history: plotClones makes plotDrivers0 unneeded with
    ## stacked and stream plots. But now so with line plot.
    ## When showing genotypes, plotDrivers0 with line only used for
    ## showing the legend.
    if(plotClones) {
        plotClonesSt(x,
                     ndr = ndr,
                     show = show,
                     na.subs = TRUE,
                     log = log,
                     lwd = lwdClone,
                     lty = ifelse(show == "drivers", ltyClone, ltyDrivers),
                     col = col, 
                     order.method = order.method,
                     stream.center = stream.center,
                     stream.frac.rand = stream.frac.rand,
                     stream.spar = stream.spar,
                     border = border,
                     srange = srange,
                     vrange = vrange,
                     type = type,
                     breakSortColors = breakSortColors,
                     colauto = colauto,
                     legend.ncols = legend.ncols,
                     lwdStackedStream = lwdStackedStream,
                     xlab = xlab,
                     ylab = ylab,
                     ylim = ylim,
                     xlim = xlim,
                     ...)
    }

    if(plotClones && plotDrivers && (type == "line"))
        par(new = TRUE)
    
    if( plotDrivers && (type == "line") ) {
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
                     log = log, ylim = ylim,
                     xlim = xlim,
                     legend.ncols = legend.ncols,
                     ...)
    }
    if(plotDiversity) {
        par(oppd)
    }
    
}

plotClonesSt <- function(z,
                         ndr,
                         show = "drivers",
                         na.subs = TRUE,
                         log = "y",
                         lwd = 1,
                         ## type = "l",
                         lty = 1:8, col = 1:9,
                         order.method = "as.is",
                         stream.center = TRUE,
                         stream.frac.rand = 0.01,
                         stream.spar = 0.2,
                         border = NULL,
                         srange = c(0.4, 1),
                         vrange = c(0.8, 1),
                         type = "stacked",
                         breakSortColors = "oe",
                         colauto = TRUE,
                         legend.ncols = "auto",
                         lwdStackedStream = 1,
                         xlab = "Time units",
                         ylab = "Number of cells",
                         ylim = NULL,
                         xlim = NULL,
                         ...) {

    ## if given ndr, we order columns based on ndr, so clones with more
    ## drivers are plotted last

    y <- z$pops.by.time[, 2:ncol(z$pops.by.time), drop = FALSE]

    ## Code in stacked and stream plots relies on there being no NAs. Could
    ## change it, but it does not seem reasonable.
    ##  But my original plotting code runs faster and is simpler if 0 are
    ##  dealt as NAs (which also makes log transformations simpler).
    
    if(type %in% c("stacked", "stream") )
        na.subs <- FALSE
    
    if(na.subs){
        y[y == 0] <- NA
    }
    ## if(is.null(ndr))
    ##     stop("Should never have null ndr")
    ## if(!is.null(ndr)) {
        ## could be done above, to avoid creating
        ## more copies
    oo <- order(ndr)
    y <- y[, oo, drop = FALSE]
    ndr <- ndr[oo]
    if(show == "drivers") {
        col <- rep(col, length.out = (1 + max(ndr)))[ndr + 1]
        lty <- rep(lty, length.out = ncol(y))
    } else {
        if(length(col) < max(ndr))
            warning("Repeating colors; you might want to",
                    "pass a col vector of more elements")
        col <- rep(col, length.out = (max(ndr)))[ndr]
    }
    ## }
    if(type == "line") {
        matplot(x = z$pops.by.time[, 1],
                y = y,
                log = log, type = "l",
                col = col, lty = lty,
                lwd = lwd,
                xlab = xlab,
                ylab = ylab,
                ylim = ylim,
                xlim = xlim,
                ...)
        box()
        if(show == "genotypes") {
            if(!inherits(z, "oncosimul2")) {
                ldrv <- genotypeLabel(z)
            } else {
                ldrv <- z$GenotypesLabels
            }
            ldrv[ldrv == ""] <- "WT"
            ldrv[ldrv == " _ "] <- "WT"
            if(legend.ncols == "auto") {
                if(length(ldrv) > 6) legend.ncols <- 2
                else legend.ncols <- 1
            }
            legend(x = "topleft",
                   title = "Genotypes",
                   lty = lty,
                   col = col, lwd = lwd,
                   legend = ldrv,
                   ncol = legend.ncols)
        }
    } else {
        ymax <- colSums(y)
        if((show == "drivers") || ((show == "genotypes") && (colauto))) {
            cll <- myhsvcols(ndr, ymax, srange = srange, vrange = vrange,
                             breakSortColors = breakSortColors)
        } else {
            cll <- list(colors = col)
        }
        x <- z$pops.by.time[, 1]
        if(grepl("y", log)) {
            stop("It makes little sense to do a stacked/stream",
                 "plot after taking the log of the y data.")
        }
        if(grepl("x", log)) {
            x <- log10(x + 1)
        }
        
        if (type == "stacked") {
            plot.stacked2(x = x,
                          y = y,
                          order.method = order.method,
                          border = border,
                          lwd = lwdStackedStream,
                          col = cll$colors,
                          log = log,
                          xlab = xlab,
                          ylab = ylab,
                          ylim = ylim,
                          xlim = xlim,
                          ...) 
        } else if (type == "stream") {
            plot.stream2(x = x,
                         y = y,
                         order.method = order.method,
                         border = border,
                         lwd = lwdStackedStream,
                         col = cll$colors,
                         frac.rand = stream.frac.rand,
                         spar = stream.spar,
                         center = stream.center,
                         log = log,
                         xlab = xlab,
                         ylab = ylab,
                         ylim = ylim,
                         xlim = xlim,
                         ...)
        } else if(type = "fish") {
            fish <- createFishObject(frac.table = t(y),
                                     parents = parentfish,
                                     timepoints = x)

        }
        if(show == "drivers") {
            if(legend.ncols == "auto") {
                if(length(cll$colorsLegend$Drivers) > 6) legend.ncols <- 2
                else legend.ncols <- 1
            }
            legend(x = "topleft",
                   title = "Number of drivers",
                   pch = 15,
                   ## lty = 1,
                   ## lwd = 2,
                   col = cll$colorsLegend$Color,
                   legend = cll$colorsLegend$Drivers,
                   ncol = legend.ncols)
        } else if (show == "genotypes") {
            if(!inherits(z, "oncosimul2")) {
                ldrv <- genotypeLabel(z)
            } else {
                ldrv <- z$GenotypesLabels
            }
            ldrv[ldrv == ""] <- "WT"
            ldrv[ldrv == " _ "] <- "WT"            
            if(legend.ncols == "auto") {
                if(length(ldrv) > 6) legend.ncols <- 2
                else legend.ncols <- 1
            }
            legend(x = "topleft",
                   title = "Genotypes",
                   pch = 15,
                   col = cll$colors,
                   legend = ldrv,
                   ncol = legend.ncols)
        }
    }
}





relabelLogaxis <- function(side = 2) {
    ## we do a plot but pass the already transformed data; relabel
    ## afterwards.
    po <- axis( side = side, labels = FALSE, tick = FALSE, lwd = 0)
    axis(side = side, labels = 10^po, at = po, tick = TRUE)
}




myhsvcols <- function(ndr, ymax, srange = c(0.4, 1),
                      vrange = c(0.8, 1),
                      breakSortColors = "oe") {
    ## Generate a set of colors so that:
    ##  - easy to tell when we increase number of drivers
    ##  - reasonably easy to tell in a legend
    ##  - different clones with same number of drivers have "similar" colors

    ## I use hsv color specification as this seems the most reasonable.
    
    minor <- table(ndr)
    major <- length(unique(ndr)) ## yeah same as length(minor), but least
                                 ## surprise
    
    h <- seq(from = 0, to = 1, length.out = major + 1)[-1]
    ## do not keep similar hues next to each other
    if(breakSortColors == "oe") {
        oe <- seq_along(h) %% 2
        h <- h[order(oe, h)]
    } else if(breakSortColors == "distave"){
        sl <- seq_along(h)
        h <- h[order(-abs(mean(sl) - sl))]
    } else if(breakSortColors == "random") {
        rr <- order(runif(length(h)))
        h <- h[rr]
    } 
    
    hh <- rep(h, minor)
    
    sr <- unlist(lapply(minor, function(x) 
        seq(from = srange[1], to = srange[2], length.out = x)))
    sv <- unlist(lapply(minor, function(x) 
        seq(from = vrange[1], to = vrange[2], length.out = x))
        )

    colors <- hsv(hh, sr, sv)

    ## This gives "average" or "median" color for legend
    ## colorsLegend <- aggregate(list(Color = colors), list(Drivers = ndr),
    ##                           function(x)
    ##                               as.character(x[((length(x) %/% 2) + 1 )]))

    ## Give the most abundant class color as the legend. Simpler to read
    colorsLegend <- by(data.frame(Color = colors, maxnum = ymax),
                       list(Drivers = ndr),
                       function(x) as.character(x$Color[which.max(x$maxnum)]))
    colorsLegend <- data.frame(Drivers = as.integer(row.names(colorsLegend)),
                               Color = cbind(colorsLegend)[, 1],
                               stringsAsFactors = FALSE)
    ## To show what it would look like
    ## plot(1:(sum(minor)), col = colors, pch = 16, cex = 3)
    ## legend(1, length(ndr), col = colorsLegend$Color, legend = names(minor),
    ##        pch = 16)

    return(list(colors = colors,
                colorsLegend = colorsLegend))
}

genotypeLabel <- function(x) {
    ## For oncosimul objects, that do not have the GenotypesLabels object
    apply(x$Genotypes, 2,
          function(x) paste(which(x == 1), sep = "", collapse = ", "))
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
                         legend.ncols = "auto",
                         ...) {
    ## z <- create.drivers.by.time(x, numDrivers)
    z <- create.drivers.by.time(x, ndr)
    ## We can now never enter here because trim.no.drivers is always FALSE
    ## in call.
    ## if(trim.no.drivers && x$MaxDriversLast) {
    ##     fi <- which(apply(z[, -c(1, 2), drop = FALSE], 1,
    ##                       function(x) sum(x) > 0))[1]
    ##     z <- z[fi:nrow(z), , drop = FALSE]
    ## }
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
    ## Allow for the weird case of possibly missing drivers
    col2 <- rep(col, length.out = (1 + max(ndr)))
    col2 <- col2[sort(unique(ndr)) + 1]
    matplot(x = time,
            y = y,
            type = type, log = log, lty = lty, col = col2, lwd = lwd,
            ...)
    if(addtot) {
        tot <- rowSums(y, na.rm = TRUE)
        lines(time, tot, col = "black", lty = 1, lwd = addtotlwd)
    }
    
    ## This will work even if under the weird case of a driver missing
    ldrv <- unlist(lapply(strsplit(colnames(y), "dr_", fixed = TRUE),
                          function(x) x[2]))
    if(legend.ncols == "auto") {
        if(length(ldrv) > 6) legend.ncols <- 2
        else legend.ncols <- 1
    }
    legend(x = "topleft",
           title = "Number of drivers",
           lty = lty, col = col2, lwd = lwd,
           legend = ldrv,
           ncol = legend.ncols)
    ## legend = (1:ncol(y)) - 1)
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

    if( (length(tG) == 1) && (tG == "")) {
        warning("There never was a descendant of WT")
    }
    
    df <- x$other$PhylogDF
    if(nrow(df) == 0) {
        warning("PhylogDF has 0 rows: no descendants of initMutant ever appeared. ",
                "This also happens if you did not set 'keepPhylog = TRUE'.")
        return(NA)
    }
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
        stop("It seems you run the simulation with keepPhylog= FALSE. ",
             "This error can also appear if your simulation exited ",
             "very fast, before any clones beyond the initial were ",
             "generated.")
    pc <- phylogClone(x, N, t, keepEvents)
    ## if(is.na(pc)) {
    ##     ## This should not be reachable, as caught before
    ##     ## where we check for nrow of PhylogDF   
    ##     warning("No clone phylogeny available. Exiting without plotting.")
    ##     return(NULL)
    ## }
        
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

closest_time <- function(x, size) {
    ## Find the first time when a given size is reached
    sizes <- rowSums(x$pops.by.time[, -1, drop = FALSE])
    candidates <- which(sizes >= size)
    if(length(candidates) == 0) {
        warning(paste("Pop size never >= requested size.",
                      "Thus, there is nothing to sample. You will get NAs"))
        return(-99)
    } else {
        return(candidates[1])
    }
}


get.the.time.for.sample <- function(tmp, timeSample, popSizeSample) {
    if( !is.null(popSizeSample) && (popSizeSample >= 0) )  {
        the.time <- closest_time(tmp, popSizeSample)
    } else if(timeSample == "last") {
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
                           thresholdWhole, popSizeSample) {
    if(is.null(x$TotalPopSize)) {
        warning(paste("It looks like this simulation never completed.",
                      " Maybe it reached maximum allowed limits.",
                      " You will get NAs"))
        return(rep(NA, length(x$geneNames)))
    }
    the.time <- get.the.time.for.sample(x, timeSample, popSizeSample)
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
                               mutationPropGrowth, ## make it explicit
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
    if(numGenes < 2)
        stop("There must be at least two genes (loci) in the fitness effects.",
             "If you only care about a degenerate case with just one,",
             "you can enter a second gene",
             "with fitness effect of zero.")
    if(numGenes > 64)
        stop("Largest possible number of genes (loci) is 64 for version 1.",
             "You are strongly encouraged to use the new specification",
             "as in version 2.")

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
        BNB_Algo5(restrictTable = rtC,
        numDrivers = numDrivers,
        numGenes = numGenes,
        typeCBN_= typeCBN,
        s = s, 
        death = death,
        mu = mu,
        initSize = initSize,
        sampleEvery = sampleEvery,
        detectionSize = detectionSize,
        finalTime = finalTime,
        initSp = initSize_species,
        initIt = initSize_iter,
        seed = seed,
        verbosity = verbosity,
        speciesFS = speciesFS,
        ratioForce = ratioForce,
        typeFitness_ = typeFitness,
        maxram = max.memory,
        mutationPropGrowth = as.integer(mutationPropGrowth),
        initMutant = initMutant,
        maxWallTime = max.wall.time,
        keepEvery = keepEvery,
        sh = sh,
        K = K,
        detectionDrivers = detectionDrivers,
        onlyCancer = onlyCancer,
        errorHitWallTime = errorHitWallTime,
        maxNumTries = max.num.tries,
        errorHitMaxTries = errorHitMaxTries,
        minDetectDrvCloneSz = minDetectDrvCloneSz,
        extraTime = extraTime
    ),
    NumDrivers = numDrivers
             ))

}

OncoSimulWide2Long <- function(x) {
    ## Put data in long format, for ggplot et al
    
    if(!inherits(x, "oncosimul2")) {
        ndr <- colSums(x$Genotypes[1:x$NumDrivers, , drop = FALSE])
        genotLabels <- genotypeLabel(x)
    } else {
        ndr <- colSums(x$Genotypes[x$Drivers, , drop = FALSE])
        genotLabels <- x$GenotypesLabels
    }
    genotLabels[genotLabels == ""] <- "WT"
    genotLabels[genotLabels == " _ "] <- "WT"
    y <- x$pops.by.time[, 2:ncol(x$pops.by.time), drop = FALSE]
    y[y == 0] <- NA
    
    oo <- order(ndr)
    y <- y[, oo, drop = FALSE]
    ndr <- ndr[oo]

    nc <- ncol(y)
    nr <- nrow(y)
    y <- as.vector(y)
    return(data.frame(Time = rep(x$pops.by.time[, 1], nc),
                      Y = y,
                      Drivers = factor(rep(ndr, rep(nr, nc))),
                      Genotype = rep(genotLabels, rep(nr, nc))))
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

xlim.pop.data <- function(x, xlim) {
    x$pops.by.time <- x$pops.by.time[
    (x$pops.by.time[, 1] >= xlim[1]) &
    (x$pops.by.time[, 1] <= xlim[2]),   ]
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



is_null_na <- function(x) {
    ## For arguments, if user passes "undefined" or "not set"
    ## See also http://stackoverflow.com/a/19655909
    if(is.function(x)) return(FALSE)
    if( is.null(x) ||
        ( (length(x) == 1) && (is.na(x)) ) ||
        ( (length(x) == 1) && (x == "") ) ## careful here
       )  {
        return(TRUE)
    } else {
        return(FALSE)
    }
}


## Not used anymore, but left here in case they become useful.
## Expected numbers at equilibrium under McFarland's
## eFinalMf <- function(initSize, s, j) {
##     ## Expected final sizes for McF, when K is set to the default.
##     # j is number of drivers
##     ## as it says, with no passengers
##     ## Set B(d) = D(N)
##     K <- initSize/(exp(1) - 1)
##     return(K * (exp( (1 + s)^j) - 1))
## }

## mcflE <- function(p, s, initSize) {
##     K <- initSize/(exp(1) - 1)
##     ## Expected number at equilibrium
##     return( K * (exp((1 + s)^p) - 1))
## }

## mcflEv <- function(p, s, initSize) {
##     ## expects vectors for p and s
##     K <- initSize/(exp(1) - 1)
    
##     ## Expected number at equilibrium
##     return( K * (exp(prod((1 + s)^p)) - 1))
## }












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


## plotClones <- function(z, ndr = NULL, na.subs = TRUE,
##                        log = "y", type = "l",
##                        lty = 1:8, col = 1:9, ...) {

##     ## if given ndr, we order columns based on ndr, so clones with more
##     ## drivers are plotted last

##     y <- z$pops.by.time[, 2:ncol(z$pops.by.time), drop = FALSE]
    
##     if(na.subs){
##         y[y == 0] <- NA
##     }
##     if(!is.null(ndr)) {
##         ## could be done above, to avoid creating
##         ## more copies
##         oo <- order(ndr)
##         y <- y[, oo, drop = FALSE]
##         ndr <- ndr[oo]
##         col <- col[ndr + 1]
##     }
##     matplot(x = z$pops.by.time[, 1],
##             y = y,
##             log = log, type = type,
##             col = col, lty = lty,
##             ...)
##     box()
## }





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
