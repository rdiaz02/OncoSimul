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



## Add tests of min_successive_fixation
## Add tolerance
## Add to help files
##   - no sense with exponential
##   - list can be used
##   examples of getting gene combs not genotypes

## documentation: change cPDetect?


## min_succesive_fixation

library(OncoSimulR)

        initS <- 2000
    r1 <- matrix(0, ncol = 6, nrow = 9)
    colnames(r1) <- c(LETTERS[1:5], "Fitness")
    r1[1, 6] <- 1
    r1[cbind((2:4), c(1:3))] <- 1
    r1[2, 6] <- 1.4
    r1[3, 6] <- 1.32
    r1[4, 6] <- 1.32
    r1[5, ] <- c(0, 1, 0, 0, 1, 1.5)
    r1[6, ] <- c(0, 0, 1, 1, 0, 1.54)
    r1[7, ] <- c(1, 0, 1, 1, 0, 1.65)
    r1[8, ] <- c(1, 1, 1, 1, 0, 1.75)
    r1[9, ] <- c(1, 1, 1, 1, 1, 1.85)
    class(r1) <- c("matrix", "genotype_fitness_matrix")
    ## plot(r1) ## to see the fitness landscape
    local_max_g <- c("A", "B, E", "A, B, C, D, E")
    local_max <- paste0("_,", local_max_g)
    fr1 <- allFitnessEffects(genotFitness = r1, drvNames = LETTERS[1:5])
    ## we pass sets of genes, so not stopping at genotypes



r0 <- oncoSimulIndiv(
                       fp = fr1,
                       model = "McFL",
                       initSize = initS,
                       mu = 1e-4,
                       detectionSize = NA,
                       sampleEvery = .03,
                       keepEvery = 1, 
                       finalTime = 50000,
                       fixation = local_max_g, 
                       detectionDrivers = NA,
                       detectionProb = NA,
                       onlyCancer = TRUE,
                       max.num.tries = 500,
                       max.wall.time = 20, 
                       errorHitMaxTries = TRUE,
    keepPhylog = FALSE)

## if you run with the C++ with the DP1 spitting out,
## you will see the differences

set.seed(1)
r00 <-      oncoSimulIndiv(
    fp = fr1,
    model = "McFL",
    initSize = initS,
    mu = 1e-4,
    detectionSize = NA,
    sampleEvery = .03,
    keepEvery = 1, 
    finalTime = 50000,
    fixation = local_max, 
    detectionDrivers = NA,
    detectionProb = NA,
    onlyCancer = TRUE,
    max.num.tries = 500,
    max.wall.time = 20, 
    errorHitMaxTries = TRUE,
    keepPhylog = FALSE)

set.seed(1)
r000 <-
    oncoSimulIndiv(
    fp = fr1,
    model = "McFL",
    initSize = initS,
    mu = 1e-4,
    detectionSize = NA,
    sampleEvery = .03,
    keepEvery = 1, 
    finalTime = 50000,
    fixation = c(local_max, min_successive_fixation = 300),
    detectionDrivers = NA,
    detectionProb = NA,
    onlyCancer = TRUE,
    max.num.tries = 500,
    max.wall.time = 20, 
    errorHitMaxTries = TRUE,
    keepPhylog = FALSE)




r1 <- oncoSimulPop(50,
                       fp = fr1,
                       model = "McFL",
                       initSize = initS,
                       mu = 1e-4,
                       detectionSize = NA,
                       sampleEvery = .03,
                       keepEvery = 1, 
                       finalTime = 50000,
                       fixation = local_max_g, 
                       detectionDrivers = NA,
                       detectionProb = NA,
                       onlyCancer = TRUE,
                       max.num.tries = 500,
                       max.wall.time = 20, 
                       errorHitMaxTries = TRUE,
                       keepPhylog = FALSE,
                       mc.cores = 2)
    sr1 <- summary(r1)
    expect_true(!all(sr1$TotalPopSize == sr1$LargestClone))
    sp1 <- samplePop(r1, "last", "singleCell")
    sgsp1 <- sampledGenotypes(sp1)
    expect_true(!all(sgsp1$Genotype %in% local_max_g))
    ## genotypes, exactly
    mm <- rep(1e-4, 5)
    mm[3] <- 2e-4
    mm[1] <- 1e-5
    names(mm) <- LETTERS[1:5]
    r2 <- oncoSimulPop(50,
                       fp = fr1,
                       model = "McFL",
                       initSize = initS,
                       mu = 1e-4,
                       detectionSize = NA,
                       sampleEvery = .03,
                       keepEvery = 1, 
                       finalTime = 50000,
                       fixation = local_max, 
                       detectionDrivers = NA,
                       detectionProb = NA,
                       onlyCancer = TRUE,
                       max.num.tries = 500,
                       max.wall.time = 20, 
                       errorHitMaxTries = TRUE,
                       keepPhylog = FALSE,
                       mc.cores = 2)
    sr2 <- summary(r2)
    expect_true(all(sr2$TotalPopSize == sr2$LargestClone))
    sp2 <- samplePop(r2, "last", "singleCell")
    sgsp2 <- sampledGenotypes(sp2)
    expect_true(all(sgsp2$Genotype %in% local_max_g))
    ## tolerance

r3 <- oncoSimulPop(50,
                       fp = fr1,
                       model = "McFL",
                       initSize = initS,
                       mu = 1e-4,
                       detectionSize = NA,
                       sampleEvery = .03,
                       keepEvery = 1, 
                       finalTime = 50000,
                       fixation = c(local_max, fixation_tolerance = 0.1),
                       detectionDrivers = NA,
                       detectionProb = NA,
                       onlyCancer = TRUE,
                       max.num.tries = 500,
                       max.wall.time = 20, 
                       errorHitMaxTries = TRUE,
                       keepPhylog = FALSE,
                       mc.cores = 2)

sr3 <- summary(r3)
    expect_true(!all(sr3$TotalPopSize == sr3$LargestClone))
    sp3 <- samplePop(r3, "last", "singleCell")
    sgsp3 <- sampledGenotypes(sp3)
    expect_true(!all(sgsp3$Genotype %in% local_max_g))


#### Can we do "fixation for anything?"
####  In C++ code, ensure it is the same genotype?

library(OncoSimulR)

        initS <- 2000
    r1 <- matrix(0, ncol = 6, nrow = 9)
    colnames(r1) <- c(LETTERS[1:5], "Fitness")
    r1[1, 6] <- 1
    r1[cbind((2:4), c(1:3))] <- 1
    r1[2, 6] <- 1.4
    r1[3, 6] <- 1.32
    r1[4, 6] <- 1.32
    r1[5, ] <- c(0, 1, 0, 0, 1, 1.5)
    r1[6, ] <- c(0, 0, 1, 1, 0, 1.54)
    r1[7, ] <- c(1, 0, 1, 1, 0, 1.65)
    r1[8, ] <- c(1, 1, 1, 1, 0, 1.75)
    r1[9, ] <- c(1, 1, 1, 1, 1, 1.85)
    class(r1) <- c("matrix", "genotype_fitness_matrix")
    ## plot(r1) ## to see the fitness landscape
    local_max_g <- c("A", "B, E", "A, B, C, D, E")
    local_max <- paste0("_,", local_max_g)
    fr1 <- allFitnessEffects(genotFitness = r1, drvNames = LETTERS[1:5])
    ## we pass sets of genes, so not stopping at genotypes

anything <- LETTERS[1:5]


r2 <- oncoSimulPop(50,
                   fp = fr1,
                   model = "McFL",
                   initSize = initS,
                   mu = 1e-4,
                   detectionSize = NA,
                   sampleEvery = .03,
                   keepEvery = 1, 
                   finalTime = 50000,
                   fixation = c(anything,
                                fixation_tolerance = 1e-4,
                                min_successive_fixation = 500),
                   detectionDrivers = NA,
                   detectionProb = NA,
                   onlyCancer = TRUE,
                   max.num.tries = 500,
                   max.wall.time = 20, 
                   errorHitMaxTries = TRUE,
                   keepPhylog = FALSE,
                   mc.cores = 2)

(sr2 <- summary(r2))



r3 <- oncoSimulIndiv(fp = fr1,
                   model = "McFL",
                   initSize = initS,
                   mu = 1e-4,
                   detectionSize = NA,
                   sampleEvery = .03,
                   keepEvery = 1, 
                   finalTime = 50000,
                   fixation = c(anything,
                                fixation_tolerance = 1e-4,
                                min_successive_fixation = 500),
                   detectionDrivers = NA,
                   detectionProb = NA,
                   onlyCancer = TRUE,
                   max.num.tries = 500,
                   max.wall.time = 20, 
                   errorHitMaxTries = TRUE,
                   keepPhylog = FALSE)

sr3 <- summary(r3)




r5 <- oncoSimulIndiv(fp = fr1,
                   model = "McFL",
                   initSize = initS,
                   mu = 1e-4,
                   detectionSize = NA,
                   sampleEvery = .03,
                   keepEvery = 1, 
                   finalTime = 50000,
                   fixation = c(local_max,
                                fixation_tolerance = 1e-4,
                                min_successive_fixation = 50,
                                fixation_min_size = 6000),
                   detectionDrivers = NA,
                   detectionProb = NA,
                   onlyCancer = TRUE,
                   max.num.tries = 500,
                   max.wall.time = 20, 
                   errorHitMaxTries = TRUE,
                   keepPhylog = FALSE)
summary(r5)





r3 <- oncoSimulIndiv(
                       fp = fr1,
                       model = "McFL",
                       initSize = initS,
                       mu = 1e-4,
                       detectionSize = NA,
                       sampleEvery = .03,
                       keepEvery = 1, 
                       finalTime = 50000,
                       fixation = c(local_max, fixation_tolerance = 0.1),
                       detectionDrivers = NA,
                       detectionProb = NA,
                       onlyCancer = TRUE,
                       max.num.tries = 500,
                       max.wall.time = 20, 
                       errorHitMaxTries = TRUE,
                       keepPhylog = FALSE)






r1 <- matrix(0, ncol = 6, nrow = 9)
colnames(r1) <- c(LETTERS[1:5], "Fitness")
r1[1, 6] <- 1
r1[cbind((2:4), c(1:3))] <- 1
r1[2, 6] <- 1.4
r1[3, 6] <- 1.32
r1[4, 6] <- 1.32
r1[5, ] <- c(0, 1, 0, 0, 1, 1.5)
r1[6, ] <- c(0, 0, 1, 1, 0, 1.54)
r1[7, ] <- c(1, 0, 1, 1, 0, 1.65)
r1[8, ] <- c(1, 1, 1, 1, 0, 1.75)
r1[9, ] <- c(1, 1, 1, 1, 1, 1.85)
class(r1) <- c("matrix", "genotype_fitness_matrix")
## plot(r1) ## to see the fitness landscape

## The local fitness maxima are
## c("A", "B, E", "A, B, C, D, E")

fr1 <- allFitnessEffects(genotFitness = r1, drvNames = LETTERS[1:5])
initS <- 2000
r3 <- oncoSimulPop(10,
                  fp = fr1,
                  model = "McFL",
                  initSize = initS,
                  mu = 1e-4,
                  detectionSize = NA,
                  sampleEvery = .03,
                  keepEvery = 1, 
                  finalTime = 50000,
                  fixation = c(paste0("_,",
                  c("A", "B, E", "A, B, C, D, E")),
                  fixation_tolerance = 0.1,
                  min_successive_fixation = 200,
                  fixation_min_size = 2300),
                  detectionDrivers = NA,
                  detectionProb = NA,
                  onlyCancer = TRUE,
                  max.num.tries = 500,
                  max.wall.time = 20, 
                  errorHitMaxTries = TRUE,
                  keepPhylog = FALSE,
                  mc.cores = 2)
summary(r3)
## All final genotypes should be local maxima                       
sp3 <- samplePop(r3, "last", "singleCell")
sgsp3 <- sampledGenotypes(sp3)


sgsp3$Genotype %in% local_max_g

