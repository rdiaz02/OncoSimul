## Run this when you change code. Load the libraries, different versions,
## in two different R.


library(OncoSimulT)
library(help = OncoSimulT)$info[[1]][c(1, 4, 5)]
data(example_trees)
max.wall.time <- 60
max.memory <- 20000

set.seed(1)
options(width = 200)

printout <- function(tmp22, rt, sh, namef) {
  np <- tmp22$TotalPopSize
  ft <- tmp22$FinalTime
  hitWallTime <- tmp22$HittedWallTime
  totalPresentD <- tmp22$TotalPresentDrivers
  MaxDriversLast <- tmp22$MaxDriversLast
  NumDrLgPopLast <- tmp22$NumDriversLargestPopLast
  PropLgPopLast <- tmp22$PropLargestPopLast
  Samples <- tmp22$outi
  ErrorMF <- tmp22$other$errorMF
  oo <- data.frame(rt, sh, np, ft, hitWallTime, totalPresentD,
                   MaxDriversLast, NumDrLgPopLast, PropLgPopLast,
                   Samples, ErrorMF, namef)

  names(oo) <- c("rt", "sh", "PopSize", "FinalTime", "hitWallTime",
                 "totalPresentDrivers", "MaxDriversLast", "NumDrLgPopLast",
                 "PropLgPopLast", "numSamples", "ErrorMF" ,"filename")
  print(oo)
  rm(oo)
}


params <- data.frame(rtmember = "rt1101", sh = -1)


runit <- function() {
    if(model == "bozic") {
        typeFitness <- "bozic1"
        detectionSize <- 1e9
        sampleEvery <- keep.every <- 0.1
        initSize <- 500
        mu <- 1e-6
        s <- 0.1
        finalTime <- (1/4) * 25 * 365
        mutatorGenotype <- 0
        birth <- -99
        death <- -99
        K <- 9999999
        ndr <- 999999
    } else if (model == "exp") {
        typeFitness <- "exp"
        detectionSize <- 1e9
        sampleEvery <- keep.every <- 0.1
        initSize <- 1000
        mu <- 1e-7
        s <- 0.1
        finalTime <- (1/4) * 25 * 365
        birth <- -99
        death <- 1
        mutatorGenotype <- 1
        K <- 9999999999
        ndr <- 99999999
    } else if(model == "mcf4"){
        typeFitness <- "mcfarlandlog"
        ndr <- 4
        sampleEvery <- 0.05
        mu <- 5e-7
        s <- 0.1
        initSize <- K <- 2000
        keep.every <- 0.1
        endTimeEvery <- 15
        detectionSize <- 1e6 ## irrelevant
        finalTime <- 15000
        mutatorGenotype <- 1 ## irrelevant
        birth <- -99
        death <- 1
    } else if(model == "mcf6"){
        typeFitness <- "mcfarlandlog"
        ndr <- 6
        sampleEvery <- 0.05
        mu <- 5e-7
        s <- 0.1
        initSize <- K <- 2000
        keep.every <- 0.1
        endTimeEvery <- 15
        detectionSize <- 1e6 ## irrelevant
        finalTime <- 15000
        mutatorGenotype <- 1 ## irrelevant
        birth <- -99
        death <- 1
    }
    if(model %in% c("mcf4", "mcf6")) {
        meetCancer <- function(x, ndr = NULL, detectionSize = NULL,
                               initSize = NULL,
                               maxPopSize = 1e15) {
            if(is.null(ndr) || is.null(initSize))
                stop("ndr and initSize must have values")
            return(
                (## (x$TotalPopSize >= detectionSize) ||
                 ( (x$MaxDriversLast >= ndr) && (x$TotalPopSize > initSize)  ) &&
                 (x$ti_dbl_min == 0) &&
                 (x$TotalPopSize < maxPopSize) ## numerical issues here
                 ))
        }
    } else if(model %in% c("exp", "bozic") ) {
        meetCancer <- function(x, ndr = NULL, detectionSize = NULL,
                               initSize = NULL, ## yes, we do nothing with this
                               maxPopSize = 1e15) {
            if(is.null(detectionSize))
                stop("detectionSize must have a value")
            return(
                (x$TotalPopSize >= detectionSize) &
                (x$ti_dbl_min == 0) &
                (x$TotalPopSize < maxPopSize) ## numerical issues here
                )
        }
    }

    for(i in 1:nrow(params)) {
        cat(paste0("\n Doing param i = ", i, " = ", params[i, 1], ", ",
                   params[i, 2]))
        n.sims <- 0
        while(n.sims < n.sims.total) {
            ## cat(paste0("\n         Doing rep ", nr))
            ## rs <<- .Random.seed ## rerun with .Random.seed <- rs
            ## cat(paste("\n rt = ", as.character(params[i, 1])))
            ## cat(paste("\n sh = ", params[i, 2], "\n"))
            tmp <- try(Algo5(restrict.table = get(as.character(params[i, 1])),
                             sh = params[i, 2],
                             initSize = initSize,
                             mu = mu,
                             sampleEvery = sampleEvery,
                             detectionSize = detectionSize,
                             max.wall.time = max.wall.time,
                             max.memory = max.memory,
                             finalTime = finalTime, 
                             s = s,
                             numGenes = 60,
                             typeFitness = typeFitness,
                             birth = birth,  death = death,
                             initMutant = -1,
                             mutatorGenotype = mutatorGenotype,
                             K = K,
                             finalDrivers = ndr,
                             typeCBN = "CBN",
                             initSize_species = 20000,
                             initSize_iter = 600,
                             ratioForce = 0.5,
                             seed_gsl = NULL,
                             speciesFS = 40000,
                             keep.every = keep.every,
                             endTimeEvery = -9,
                             alpha = -99,
                             verbosity = 0,
                             silent = FALSE))
            if(!inherits(tmp, "try-error")) {
                mtc <- meetCancer(tmp, detectionSize = detectionSize,
                                  ndr = ndr, initSize = initSize)
                if(mtc) {
                    n.sims <- n.sims + 1
                    namef <- paste(model, "_", params[i, 1],
                                   ".sh_", params[i, 2],
                                   ".MDr_", tmp$MaxDriversLast,
                                   ".LDr_", tmp$NumDriversLargestPopLast,
                                   ".Time_", round(tmp$FinalTime),
                                   ".Pop_", format(tmp$TotalPopSize, scientific = TRUE),
                                   ".ErrorMF_", round(tmp$other$errorMF, 4),
                                   "_", sep = "")
                    rnn <- paste(sample(c(letters, LETTERS, 0:9), 12, replace = TRUE),
                                 collapse = "")
                    namef <- paste(namef, rnn, ".rds", sep = "")
                    ## saveRDS(tmp, file = namef, compress = TRUE)
                    printout(tmp, params[i, 1], params[i, 2], namef)
                }
            }
        }
        ## params[i, "n.cancer"] <- n.cancer
        ## params[i, "ntot"] <- reps
        ## cat(paste0("\n    n.cancer = ", n.cancer, "\n"))
    }
}

n.sims.total <- 1
for(model in c("exp", "bozic", "mcf4", "mcf6")) {
    runit()
}


set.seed(34)
for(model in c("exp", "bozic", "mcf4", "mcf6")) {
    runit()
}


set.seed(20149854)
for(model in c("exp", "bozic", "mcf4", "mcf6")) {
    runit()
}
