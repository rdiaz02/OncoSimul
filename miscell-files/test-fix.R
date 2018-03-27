library(OncoSimulR)


## source("~/Proyectos/predictability-local-maxima/simulations-2018-03/ruggify-functions.R")
i <- 3
ng <- 7
RNGkind("Mersenne-Twister")

## set.seed(i)
## x <- single_ruggified_DAG(ng, TRUE, 50,  0.25)
## fll <- x$fitness_landscape
## fll[fll[, "Fitness"] < 1e-8, "Fitness"] <- 0
## fee <- allFitnessEffects(genotFitness = fll, drvNames = LETTERS[1:ng])
## save(file = "fee.RData", fee, x)


load("fee.RData")
set.seed(i)
r3 <- oncoSimulIndiv( fp = fee,
                     model = "McFL",
                     initSize = 2000,
                     mu = 1e-4,
                     detectionSize = NA,
                     sampleEvery = .03,
                     keepEvery = 1,
                     finalTime = 50000,
                     fixation = unlist(x[["labelled_peaks"]]), 
                     detectionDrivers = NA,
                     detectionProb = NA,
                     onlyCancer = TRUE,
                     max.num.tries = 500,
                     max.wall.time = 20, 
                     errorHitMaxTries = TRUE,
                     keepPhylog = FALSE)

summary(r3)

## stopping at ABCG, which is not a maximum, not a labelled peak

r3$pops.by.time

x$labelled_peaks
r3$GenotypesLabels[c(4, 5, 6, 8)]
r3$pops.by.time[, c(4, 5, 6, 8) + 1]




library(OncoSimulR)
## source("~/Proyectos/predictability-local-maxima/simulations-2018-03/ruggify-functions.R")
i <- 3
ng <- 7
RNGkind("Mersenne-Twister")
load("fee.RData")
set.seed(i)
r4 <- oncoSimulIndiv( fp = fee,
                     model = "McFL",
                     initSize = 2000,
                     mu = 1e-4,
                     detectionSize = NA,
                     sampleEvery = .03,
                     keepEvery = 1,
                     finalTime = 50000,
                     fixation = paste0("_,", unlist(x[["labelled_peaks"]])),
                     detectionDrivers = NA,
                     detectionProb = NA,
                     onlyCancer = TRUE,
                     max.num.tries = 500,
                     max.wall.time = 20, 
                     errorHitMaxTries = TRUE,
                     keepPhylog = FALSE)

summary(r4)
r4$pops.by.time
x$labelled_peaks


set.seed(i)
r5 <- oncoSimulIndiv( fp = fee,
                     model = "McFL",
                     initSize = 2000,
                     mu = 1e-4,
                     detectionSize = NA,
                     sampleEvery = .03,
                     keepEvery = 1,
                     finalTime = 50000,
                     fixation = c(paste0("_,", unlist(x[["labelled_peaks"]])),
                                 fixation_tolerance = 0.05),
                     detectionDrivers = NA,
                     detectionProb = NA,
                     onlyCancer = TRUE,
                     max.num.tries = 500,
                     max.wall.time = 20, 
                     errorHitMaxTries = TRUE,
                     keepPhylog = FALSE)
summary(r5)

54$pops.by.time
x$labelled_peaks















r4$GenotypesLabels[c(4, 5, 6, 8, 12, 13)]
r4$pops.by.time[160, c(4, 5, 6, 8, 12, 13) + 1]
sum(r4$pops.by.time[160, c(4, 5, 6, 8, 12, 13) + 1])



## Add tests
## Add tolerance
## Add to help files
##   - no sense with exponential
##   - list can be used
##   examples of getting gene combs not genotypes
