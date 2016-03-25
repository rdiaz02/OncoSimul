RNGkind("L'Ecuyer-CMRG") ## for the mclapplies
## If crashes I want to see where: thus output seed.
## The tests below can occasionally fail (but that probability decreases
## as we increase number of pops), as they should.

cat("\n", date(), "\n") ## whole file takes about 9 seconds
date()
test_that("Ordering of number of clones with mutpropgrowth", {
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n omp1: the seed is", pseed, "\n")
    pops <- 100
    lni <- 200
    no <- 5e3
    ni <- c(5, 2, rep(0, lni))
    names(ni) <- c("a", "b", paste0("n", seq.int(lni)))
    fe <- allFitnessEffects(noIntGenes = ni)
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n omp1a: the seed is", pseed, "\n")
    nca <- oncoSimulPop(pops, fe, finalTime =1,
                        mutationPropGrowth = TRUE,
                        initSize = no,
                        initMutant = "a",
                        onlyCancer = FALSE, seed = NULL, mc.cores = 2)
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n omp1b: the seed is", pseed, "\n")
    ncb <- oncoSimulPop(pops, fe, finalTime =1,
                        mutationPropGrowth = TRUE,
                        initSize = no,
                        initMutant = "b",
                        onlyCancer = FALSE, seed = NULL, mc.cores = 2)
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n omp1c: the seed is", pseed, "\n")
    nca2 <- oncoSimulPop(pops, fe, finalTime =1,
                        mutationPropGrowth = FALSE,
                        initSize = no,
                        initMutant = "a",
                        onlyCancer = FALSE, seed = NULL, mc.cores = 2)
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n omp1d: the seed is", pseed, "\n")
    ncb2 <- oncoSimulPop(pops, fe, finalTime =1,
                        mutationPropGrowth = FALSE,
                        initSize = no,
                        initMutant = "b",
                        onlyCancer = FALSE, seed = NULL, mc.cores = 2)
    ## I once saw a weird thing
    expect_true(var(summary(nca)$NumClones) > 1e-4)
    expect_true(var(summary(ncb)$NumClones) > 1e-4)
    expect_true(var(summary(nca2)$NumClones) > 1e-4)
    expect_true(var(summary(ncb2)$NumClones) > 1e-4)
    ## The real comparison
    expect_true( median(summary(nca)$NumClones) >
                 median(summary(nca2)$NumClones))
    expect_true( median(summary(nca)$NumClones) >
                 median(summary(ncb)$NumClones))
    expect_true( median(summary(ncb)$NumClones) >
                 median(summary(ncb2)$NumClones))
})
cat("\n", date(), "\n")

cat("\n", date(), "\n")
test_that("Ordering of number of clones with mutpropgrowth, McFL", {
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n omp2: the seed is", pseed, "\n")
    pops <- 100
    lni <- 200
    no <- 5e3
    ni <- c(5, 2, rep(0, lni))
    names(ni) <- c("a", "b", paste0("n", seq.int(lni)))
    fe <- allFitnessEffects(noIntGenes = ni)
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n omp2a: the seed is", pseed, "\n")
    nca <- oncoSimulPop(pops, fe, finalTime =1,
                        mutationPropGrowth = TRUE,
                        initSize = no, model = "McFL",
                        initMutant = "a",
                        onlyCancer = FALSE, seed = NULL, mc.cores = 2)
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n omp2b: the seed is", pseed, "\n")
    ncb <- oncoSimulPop(pops, fe, finalTime =1,
                        mutationPropGrowth = TRUE,
                        initSize = no, model = "McFL",
                        initMutant = "b",
                        onlyCancer = FALSE, seed = NULL, mc.cores = 2)
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n omp2c: the seed is", pseed, "\n")
    nca2 <- oncoSimulPop(pops, fe, finalTime =1,
                        mutationPropGrowth = FALSE,
                        initSize = no, model = "McFL",
                        initMutant = "a",
                        onlyCancer = FALSE, seed = NULL, mc.cores = 2)
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n omp2d: the seed is", pseed, "\n")
    ncb2 <- oncoSimulPop(pops, fe, finalTime =1,
                        mutationPropGrowth = FALSE,
                        initSize = no, model = "McFL",
                        initMutant = "b",
                        onlyCancer = FALSE, seed = NULL, mc.cores = 2)
    ## I once saw a weird thing
    expect_true(var(summary(nca)$NumClones) > 1e-4)
    expect_true(var(summary(ncb)$NumClones) > 1e-4)
    expect_true(var(summary(nca2)$NumClones) > 1e-4)
    expect_true(var(summary(ncb2)$NumClones) > 1e-4)
    ## The real comparison
    expect_true( median(summary(nca)$NumClones) >
                 median(summary(nca2)$NumClones))
    expect_true( median(summary(nca)$NumClones) >
                 median(summary(ncb)$NumClones))
    expect_true( median(summary(ncb)$NumClones) >
                 median(summary(ncb2)$NumClones))
})

cat("\n", date(), "\n")
test_that("oncoSimulSample Without initmutant and modules", {
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n osS: the seed is", pseed, "\n")
    pops <- 60
    lni <- 1 ## no fitness effects genes
    fni <- 50 ## fitness effects genes
    no <- 1e4 
    ft <- 4 
    s3 <- 2.5 
    mu <- 1e-5 
    ## noInt have no fitness effects, but can accumulate mutations
    ni <- rep(0, lni)
    ## Those with fitness effects in one module, so
    ## neither fitness nor mut. rate blow up
    gn <- paste(paste0("a", 1:fni), collapse = ", ")
    f3 <- allFitnessEffects(epistasis = c("A" = s3),
                            geneToModule = c("A" = gn),
                            noIntGenes = ni)
    x <- 1e-9 ## so basically anything that appears once
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n osSa: the seed is", pseed, "\n")
    b1 <- oncoSimulSample(pops,
                          f3,
                          mu = mu,
                          mutationPropGrowth = FALSE,
                          finalTime =ft,
                          initSize = no,
                          onlyCancer = FALSE,
                          sampleEvery = 0.01,
                          detectionSize = 1e9,
                          detectionDrivers = 99,
                          seed =NULL,
                          thresholdWhole = x)
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n osSb: the seed is", pseed, "\n")
    b2 <- oncoSimulSample(pops,
                         f3,
                         mu = mu,
                         mutationPropGrowth = TRUE,
                         finalTime =ft,
                         initSize = no,
                         onlyCancer = FALSE,
                         sampleEvery = 0.01,
                          detectionSize = 1e9,
                          detectionDrivers = 99,
                          seed =NULL,
                         thresholdWhole = x)
    b1$popSummary[1:5, c(1:3, 8:9)]
    summary(b1$popSummary[, "NumClones"])
    summary(b1$popSummary[, "TotalPopSize"])
    b2$popSummary[1:5, c(1:3, 8:9)]
    summary(b2$popSummary[, "NumClones"])
    summary(b2$popSummary[, "TotalPopSize"])
    ## cc1 and cc2 should all be smaller than pops, or you are maxing
    ## things and not seeing patterns
    (cc1 <- colSums(b1$popSample))
    (cc2 <- colSums(b2$popSample))
    ## Of course, these are NOT really mutationsPerClone: we collapse over
    ## whole population.
    (mutsPerClone1 <- rowSums(b1$popSample))
    (mutsPerClone2 <- rowSums(b2$popSample))
    summary(mutsPerClone1)
    summary(mutsPerClone2)
    expect_true( mean(mutsPerClone2) >
                 mean(mutsPerClone1))
    expect_true( median(b2$popSummary[, "NumClones"]) >
                 median(b1$popSummary[, "NumClones"]))
})
cat("\n", date(), "\n")


cat("\n", date(), "\n")
test_that("oncoSimulSample Without initmutant and modules, McFL", {
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n osSMcFL: the seed is", pseed, "\n")
    pops <- 60
    lni <- 1 ## no fitness effects genes
    fni <- 50 ## fitness effects genes
    no <- 1e4 ## note we use only 10 in the other example below
    ft <- 4 
    s3 <- 2.5 
    mu <- 1e-5 
    ## noInt have no fitness effects, but can accumulate mutations
    ni <- rep(0, lni)
    ## Those with fitness effects in one module, so
    ## neither fitness nor mut. rate blow up
    gn <- paste(paste0("a", 1:fni), collapse = ", ")
    f3 <- allFitnessEffects(epistasis = c("A" = s3),
                            geneToModule = c("A" = gn),
                            noIntGenes = ni)
    x <- 1e-9 ## 1/no
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n osSMcFLa: the seed is", pseed, "\n")
    b1 <- oncoSimulSample(pops,
                          f3,
                          mu = mu,
                          mutationPropGrowth = FALSE,
                          finalTime =ft,
                          initSize = no,
                          onlyCancer = FALSE,
                          sampleEvery = 0.01,
                          detectionSize = 1e9,
                          detectionDrivers = 99,
                          seed =NULL,
                          model = "McFL",
                          thresholdWhole = x)
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n osSMcFLb: the seed is", pseed, "\n")
    b2 <- oncoSimulSample(pops,
                         f3,
                         mu = mu,
                         mutationPropGrowth = TRUE,
                         finalTime =ft,
                         initSize = no,
                         onlyCancer = FALSE,
                         sampleEvery = 0.01,
                          detectionSize = 1e9,
                          detectionDrivers = 99,
                         seed =NULL,
                         model = "McFL",
                         thresholdWhole = x)
    b1$popSummary[1:5, c(1:3, 8:9)]
    summary(b1$popSummary[, "NumClones"])
    summary(b1$popSummary[, "TotalPopSize"])
    b2$popSummary[1:5, c(1:3, 8:9)]
    summary(b2$popSummary[, "NumClones"])
    summary(b2$popSummary[, "TotalPopSize"])
    ## cc1 and cc2 should all be smaller than pops, or you are maxing
    ## things and not seeing patterns
    (cc1 <- colSums(b1$popSample))
    (cc2 <- colSums(b2$popSample))
    (mutsPerClone1 <- rowSums(b1$popSample))
    (mutsPerClone2 <- rowSums(b2$popSample))
    summary(mutsPerClone1)
    summary(mutsPerClone2)
    expect_true( mean(mutsPerClone2) >
                 mean(mutsPerClone1))
    expect_true( median(b2$popSummary[, "NumClones"]) >
                 median(b1$popSummary[, "NumClones"]))
})
cat("\n", date(), "\n")






cat("\n Ended test.mutPropGrowth: ", date(), "\n")


## ##     A way to check is to see the output from the C++ code with the
## ##     verbosity option.

## RNGkind("Mersenne-Twister")
## ni <- rep(0.4, 20)
## names(ni) <- c("a", "b", "c", "d", paste0("n", 1:16))
## fe <- allFitnessEffects(noIntGenes = ni)
## set.seed(5) 
## oncoSimulIndiv(fe, finalTime =30,
##                mutationPropGrowth = TRUE,
##                initSize = 1e4,
##                mu = 1e-06,
##                verbosity = 6,
##                onlyCancer = FALSE)

## ###### Iteration 30.
## ## mutation
## ## child
## 1.4 * 1e-06 * 19

## ni <- rep(0.4, 20)
## names(ni) <- c("a", "b", "c", "d", paste0("n", 1:16))
## fe <- allFitnessEffects(noIntGenes = ni)
## set.seed(25) 
## oncoSimulIndiv(fe, finalTime =40,
##                mutationPropGrowth = TRUE,
##                initSize = 1e4,
##                mu = 1e-06,
##                verbosity = 6,
##                onlyCancer = FALSE)
## ## Iteration 48.
## ## Birth of child:
## 1.4 * 1.4
## ## Mutation of child
## 1.96 * 1e-06 * 18



