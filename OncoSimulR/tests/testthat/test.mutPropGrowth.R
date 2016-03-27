cat(paste("\n Starting at mutPropGrowth ", date(), "\n"))

RNGkind("L'Ecuyer-CMRG") ## for the mclapplies
## If crashes I want to see where: thus output seed.
## The tests below can occasionally fail (but that probability decreases
## as we increase number of pops), as they should.

cat("\n", date(), "\n") ## whole file takes about 16 seconds

mutsPerClone <- function(x, per.pop.mean = TRUE) {
    perCl <- function(z)
        unlist(lapply(z$GenotypesWDistinctOrderEff, length))
    perCl2 <- function(z)
        mean(unlist(lapply(z$GenotypesWDistinctOrderEff, length)))

    if(per.pop.mean)    
        unlist(lapply(x, function(u) perCl2(u)))
    else
        lapply(x, function(u) perCl(u))
}

## ## this tests takes 10 seconds. Moved to long.
## date()
## test_that("mutPropGrowth diffs with s> 0", {
    
##     pseed <- sample(9999999, 1)
##     set.seed(pseed)
##     cat("\n mgp1: the seed is", pseed, "\n")
##     ft <- 4 ## 2.7
##     pops <- 50
##     lni <- 150 ## 100
##     no <- 1e3 ## 5e1 
##     ni <- c(2, rep(0, lni)) ## 2 ## 4 ## 5
##     mu <- 1e-5 ## 1e-6
##     names(ni) <- c("a", paste0("n", seq.int(lni)))
##     ni <- sample(ni) ## scramble 
##     fe <- allFitnessEffects(noIntGenes = ni)
##     pseed <- sample(9999999, 1)
##     set.seed(pseed)
##     cat("\n mpg1a: the seed is", pseed, "\n")
##     nca <- oncoSimulPop(pops, fe, finalTime = ft,
##                         mu = mu,
##                         mutationPropGrowth = TRUE,
##                         initSize = no, sampleEvery = 0.1,
##                         initMutant = "a", keepEvery = 1,
##                         onlyCancer = FALSE, seed = NULL, mc.cores = 2)
##     pseed <- sample(9999999, 1)
##     set.seed(pseed)
##     cat("\n mpg1c: the seed is", pseed, "\n")
##     nca2 <- oncoSimulPop(pops, fe, finalTime = ft,
##                          mu = mu,
##                         mutationPropGrowth = FALSE,
##                         initSize = no, sampleEvery = 0.1,
##                         initMutant = "a", keepEvery = 1,
##                         onlyCancer = FALSE, seed = NULL, mc.cores = 2)
##     ## summary(nca)[1:20, c(1, 2, 3, 8, 9)]
##     ## summary(nca2)[1:20, c(1, 2, 3, 8, 9)]
##     ## I once saw a weird thing
##     expect_true(var(summary(nca)$NumClones) > 1e-4)
##     expect_true(var(summary(nca2)$NumClones) > 1e-4)
##     ## summary(summary(nca)$NumClones)
##     ## summary(summary(nca2)$NumClones)
##     ## summary(mutsPerClone(nca))
##     ## summary(mutsPerClone(nca2))
##     ## The real comparison
##     expect_true( median(summary(nca)$NumClones) >
##                  median(summary(nca2)$NumClones))
##     expect_true( mean(mutsPerClone(nca)) >
##                  mean(mutsPerClone(nca2)))
    

    
## })
## cat("\n", date(), "\n")



date()
test_that("mutPropGrowth diffs with s> 0, McFL", {
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n mcf1: the seed is", pseed, "\n")
    ft <- 3 
    pops <- 50
    lni <- 100
    no <- 1e3 ## 5e1 
    ni <- c(4, rep(0, lni)) ## 5
    names(ni) <- c("a", paste0("n", seq.int(lni)))
    ni <- sample(ni) ## scramble
    fe <- allFitnessEffects(noIntGenes = ni)
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n mcf1a: the seed is", pseed, "\n")
    nca <- oncoSimulPop(pops, fe, finalTime = ft,
                        mutationPropGrowth = TRUE,
                        initSize = no,
                        sampleEvery = 0.03,
                        keepEvery = 1,
                        initMutant = "a", model = "McFL",
                        onlyCancer = FALSE, seed = NULL, mc.cores = 2)
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n mcf1c: the seed is", pseed, "\n")
    nca2 <- oncoSimulPop(pops, fe, finalTime = ft,
                        mutationPropGrowth = FALSE,
                        initSize = no,
                        sampleEvery = 0.03,
                        keepEvery = 1,
                        initMutant = "a", model = "McFL",
                        onlyCancer = FALSE, seed = NULL, mc.cores = 2)
    ## summary(nca)[1:20, c(1, 2, 3, 8, 9)]
    ## summary(nca2)[1:20, c(1, 2, 3, 8, 9)]
    ## I once saw a weird thing
    expect_true(var(summary(nca)$NumClones) > 1e-4)
    expect_true(var(summary(nca2)$NumClones) > 1e-4)
    ## summary(summary(nca)$NumClones)
    ## summary(summary(nca2)$NumClones)
    ## summary(mutsPerClone(nca))
    ## summary(mutsPerClone(nca2))
    ## The real comparison
    expect_true( median(summary(nca)$NumClones) >
                 median(summary(nca2)$NumClones))
    expect_true( mean(mutsPerClone(nca)) >
                 mean(mutsPerClone(nca2)))
})
cat("\n", date(), "\n")




date()
test_that("mutPropGrowth diffs with s> 0, oncoSimulSample", {
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n oss1: the seed is", pseed, "\n")
    ft <- 3.5 ## 4
    pops <- 50
    lni <- 200 ## 150
    no <- 1e3 ## 5e1 
    ni <- c(2, rep(0, lni)) ## 2 ## 4 ## 5
    mu <- 5e-7 ## 1e-6
    x <- 1e-9
    names(ni) <- c("a", paste0("n", seq.int(lni)))
    ni <- sample(ni) ## scramble
    fe <- allFitnessEffects(noIntGenes = ni)
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n oss1a: the seed is", pseed, "\n")
    nca <- oncoSimulSample(pops, fe, finalTime = ft,
                        mu = mu,
                        mutationPropGrowth = TRUE,
                        initSize = no, sampleEvery = 0.02,
                        initMutant = "a", 
                        onlyCancer = FALSE, seed = NULL,
                        detectionSize = 1e9,
                        detectionDrivers = 99,
                        thresholdWhole = x)
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n oss1c: the seed is", pseed, "\n")
    nca2 <- oncoSimulSample(pops, fe, finalTime = ft,
                         mu = mu,
                        mutationPropGrowth = FALSE,
                        initSize = no, sampleEvery = 0.02,
                        initMutant = "a", 
                        onlyCancer = FALSE, seed = NULL,
                        detectionSize = 1e9,
                        detectionDrivers = 99,
                        thresholdWhole = x)
    ## nca$popSummary[1:5, c(1:3, 8:9)]
    ## summary(nca$popSummary[, "NumClones"])
    ## summary(nca$popSummary[, "TotalPopSize"])
    ## nca2$popSummary[1:5, c(1:3, 8:9)]
    ## summary(nca2$popSummary[, "NumClones"])
    ## summary(nca2$popSummary[, "TotalPopSize"])
    ## ## cc1 and cc2 should all be smaller than pops, or you are maxing
    ## ## things and not seeing patterns
    ## summary(cc1 <- colSums(nca$popSample))
    ## summary(cc2 <- colSums(nca2$popSample))
    ## which(cc1 == pops)
    ## which(cc2 == pops)
    ## Of course, these are NOT really mutationsPerClone: we collapse over
    ## whole population. Ends up being very similar to NumClones, except
    ## fr the few that go extinct.
    mutsPerClone1 <- rowSums(nca$popSample)
    mutsPerClone2 <- rowSums(nca2$popSample)
    ## summary(mutsPerClone1)
    ## summary(mutsPerClone2)
    expect_true( mean(mutsPerClone1) >
                 mean(mutsPerClone2))
    expect_true( median(nca$popSummary[, "NumClones"]) >
                 median(nca2$popSummary[, "NumClones"]))
})


date()
test_that("mutPropGrowth diffs with s> 0, oncoSimulSample, McFL", {
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n ossmcf1: the seed is", pseed, "\n")
    ft <- 40 ## 4
    pops <- 50
    lni <- 200 ## 150
    no <- 1e3 ## 5e1 
    ni <- c(2, rep(0, lni)) ## 2 ## 4 ## 5
    mu <- 5e-7 ## 1e-6
    x <- 1e-9
    names(ni) <- c("a", paste0("n", seq.int(lni)))
    ni <- sample(ni) ## scramble
    fe <- allFitnessEffects(noIntGenes = ni)
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n ossmcf1a: the seed is", pseed, "\n")
    nca <- oncoSimulSample(pops, fe, finalTime = ft,
                        mu = mu, model = "McFL",
                        mutationPropGrowth = TRUE,
                        initSize = no, sampleEvery = 0.01,
                        initMutant = "a", 
                        onlyCancer = FALSE, seed = NULL,
                        detectionSize = 1e9,
                        detectionDrivers = 99,
                        thresholdWhole = x)
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n ossmcf1c: the seed is", pseed, "\n")
    nca2 <- oncoSimulSample(pops, fe, finalTime = ft,
                         mu = mu, model = "McFL",
                        mutationPropGrowth = FALSE,
                        initSize = no, sampleEvery = 0.01,
                        initMutant = "a", 
                        onlyCancer = FALSE, seed = NULL,
                        detectionSize = 1e9,
                        detectionDrivers = 99,
                        thresholdWhole = x)
    ## nca$popSummary[1:5, c(1:3, 8:9)]
    ## summary(nca$popSummary[, "NumClones"])
    ## summary(nca$popSummary[, "TotalPopSize"])
    ## nca2$popSummary[1:5, c(1:3, 8:9)]
    ## summary(nca2$popSummary[, "NumClones"])
    ## summary(nca2$popSummary[, "TotalPopSize"])
    ## ## cc1 and cc2 should all be smaller than pops, or you are maxing
    ## ## things and not seeing patterns
    ## summary(cc1 <- colSums(nca$popSample))
    ## summary(cc2 <- colSums(nca2$popSample))
    ## which(cc1 == pops)
    ## which(cc2 == pops)
    ## Of course, these are NOT really mutationsPerClone: we collapse over
    ## whole population. Ends up being very similar to NumClones, except
    ## fr the few that go extinct.
    mutsPerClone1 <- rowSums(nca$popSample)
    mutsPerClone2 <- rowSums(nca2$popSample)
    ## summary(mutsPerClone1)
    ## summary(mutsPerClone2)
    expect_true( mean(mutsPerClone1) >
                 mean(mutsPerClone2))
    expect_true( median(nca$popSummary[, "NumClones"]) >
                 median(nca2$popSummary[, "NumClones"]))

    

    
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



