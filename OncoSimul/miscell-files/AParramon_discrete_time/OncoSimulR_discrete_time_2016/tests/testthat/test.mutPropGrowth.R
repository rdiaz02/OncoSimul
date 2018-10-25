cat(paste("\n Starting at mutPropGrowth ", date(), "\n"))

## RNGkind("L'Ecuyer-CMRG") ## for the mclapplies
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

p.value.threshold <- 0.001



date()
test_that("mutPropGrowth diffs with s> 0, McFL", {
    max.tries <- 4
    for(tries in 1:max.tries) {
        T1 <- T2 <- T3 <- T4 <- T5 <- T6 <- T7 <- T8 <- TRUE

        ## stopping on time or size is about the same
        ## in these models, but we stop on size to
        ## control for different  pop size. Now, time should be very similar, or we introduce a serious distorsion.
        
        cat("\n mcf1: a runif is", runif(1), "\n")
        ft <- 26  # 3 
        pops <- 50
        lni <- 100
        no <- 1e3 ## 5e1 
        ni <- c(4, rep(0, lni)) ## 5
        names(ni) <- c("a", paste0("n", seq.int(lni)))
        ni <- sample(ni) ## scramble
        fe <- allFitnessEffects(noIntGenes = ni)
        cat("\n mcf1a: a runif is", runif(1), "\n")
        nca <- oncoSimulPop(pops, fe, finalTime = ft, detectionProb = NA,
                            mutationPropGrowth = TRUE,
                            initSize = no,
                            sampleEvery = 0.03, detectionSize = 8e4,
                            keepEvery = 1,
                            initMutant = "a", model = "McFL",
                            onlyCancer = FALSE, seed = NULL, mc.cores = 2)
        cat("\n mcf1c: a runif is", runif(1), "\n")
        nca2 <- oncoSimulPop(pops, fe, finalTime = ft, detectionProb = NA,
                             mutationPropGrowth = FALSE,
                             initSize = no,
                             sampleEvery = 0.03, detectionSize = 8e4,
                             keepEvery = 1,
                             initMutant = "a", model = "McFL",
                             onlyCancer = FALSE, seed = NULL, mc.cores = 2)
        summary(nca)[1:20, c(1, 2, 3, 8, 9)]
        summary(nca2)[1:20, c(1, 2, 3, 8, 9)]
        ## I once saw a weird thing
        expect_true(var(summary(nca)$NumClones) > 1e-4)
        expect_true(var(summary(nca2)$NumClones) > 1e-4)
        ## Pop sizes here do not really differ, as we start with the
        ## initMutant with increased s
        summary(summary(nca)[, 2])
        summary(summary(nca2)[, 2])
        ## summary(summary(nca)$NumClones)
        ## summary(summary(nca2)$NumClones)
        ## summary(mutsPerClone(nca))
        ## summary(mutsPerClone(nca2))
        ## The real comparison
        T1 <- ( wilcox.test(summary(nca)$NumClones ,
                     summary(nca2)$NumClones, alternative = "greater")$p.value < p.value.threshold)
        T2 <- ( t.test(mutsPerClone(nca) ,
                mutsPerClone(nca2), alternative = "greater")$p.value < p.value.threshold)

        
        if(T1 && T2 && T3 && T4 && T5 && T6 && T7 && T8) break;
    }
    cat(paste("\n done tries", tries, "\n"))
    expect_true(T1 && T2 && T3 && T4 && T5 && T6 && T7 && T8)
})
cat("\n", date(), "\n")



date()
test_that("mutPropGrowth diffs with s> 0, McFL, stop on time", {
    max.tries <- 4
    for(tries in 1:max.tries) {
        T1 <- T2 <- T3 <- T4 <- T5 <- T6 <- T7 <- T8 <- TRUE

        ## stopping on time or size is about the same in these models, but
        ## we stop on time now. Note that sizes are about the same as we
        ## have been in the plateau for long
        
        
        cat("\n mcf1_ontime: a runif is", runif(1), "\n")
        ft <- 6 
        pops <- 50
        lni <- 100
        no <- 1e3 ## 5e1 
        ni <- c(4, rep(0, lni)) ## 5
        ds <- 1e9
        names(ni) <- c("a", paste0("n", seq.int(lni)))
        ni <- sample(ni) ## scramble
        fe <- allFitnessEffects(noIntGenes = ni)
        cat("\n mcf1a: a runif is", runif(1), "\n")
        nca <- oncoSimulPop(pops, fe, finalTime = ft, detectionProb = NA,
                            mutationPropGrowth = TRUE,
                            initSize = no,
                            sampleEvery = 0.03, detectionSize = ds,
                            keepEvery = 1,
                            initMutant = "a", model = "McFL",
                            onlyCancer = FALSE, seed = NULL, mc.cores = 2)
        cat("\n mcf1c: a runif is", runif(1), "\n")
        nca2 <- oncoSimulPop(pops, fe, finalTime = ft, detectionProb = NA,
                             mutationPropGrowth = FALSE,
                             initSize = no,
                             sampleEvery = 0.03, detectionSize = ds,
                             keepEvery = 1,
                             initMutant = "a", model = "McFL",
                             onlyCancer = FALSE, seed = NULL, mc.cores = 2)
        summary(nca)[1:20, c(1, 2, 3, 8, 9)]
        summary(nca2)[1:20, c(1, 2, 3, 8, 9)]
        ## I once saw a weird thing
        expect_true(var(summary(nca)$NumClones) > 1e-4)
        expect_true(var(summary(nca2)$NumClones) > 1e-4)
        ## Pop sizes here do not really differ, as we start with the
        ## initMutant with increased s
        summary(summary(nca)[, 2])
        summary(summary(nca2)[, 2])
        ## summary(summary(nca)$NumClones)
        ## summary(summary(nca2)$NumClones)
        ## summary(mutsPerClone(nca))
        ## summary(mutsPerClone(nca2))
        ## The real comparison
        T1 <- ( wilcox.test(summary(nca)$NumClones ,
                     summary(nca2)$NumClones, alternative = "greater")$p.value < p.value.threshold)
        T2 <- ( t.test(mutsPerClone(nca) ,
                mutsPerClone(nca2), alternative = "greater")$p.value < p.value.threshold)

        
        if(T1 && T2 && T3 && T4 && T5 && T6 && T7 && T8) break;
    }
    cat(paste("\n done tries", tries, "\n"))
    expect_true(T1 && T2 && T3 && T4 && T5 && T6 && T7 && T8)
})
cat("\n", date(), "\n")



date()
test_that("mutPropGrowth diffs with s> 0, oncoSimulSample", {
    max.tries <- 4
    for(tries in 1:max.tries) {
        T1 <- T2 <- T3 <- T4 <- T5 <- T6 <- T7 <- T8 <- TRUE
        cat("\n oss1: a runif is", runif(1), "\n")
        
        ft <- 133 ## we stop on size way earlier 
        pops <- 50
        lni <- 200 ## 150
        no <- 1e3 ## 5e1 
        ni <- c(2, rep(0, lni)) ## 2 ## 4 ## 5
        mu <- 5e-7 ## 1e-6
        x <- 1e-9
        names(ni) <- c("a", paste0("n", seq.int(lni)))
        ni <- sample(ni) ## scramble
        fe <- allFitnessEffects(noIntGenes = ni)
        cat("\n oss1a: a runif is", runif(1), "\n")
        nca <- oncoSimulSample(pops, fe, finalTime = ft, detectionProb = NA,
                               mu = mu,
                               mutationPropGrowth = TRUE,
                               initSize = no, sampleEvery = 0.02,
                               initMutant = "a", 
                               onlyCancer = FALSE, seed = NULL,
                               detectionSize = 1e6,
                               detectionDrivers = 99,
                               thresholdWhole = x)
        cat("\n oss1c: a runif is", runif(1), "\n")
        nca2 <- oncoSimulSample(pops, fe, finalTime = ft, detectionProb = NA,
                                mu = mu,
                                mutationPropGrowth = FALSE,
                                initSize = no, sampleEvery = 0.02,
                                initMutant = "a", 
                                onlyCancer = FALSE, seed = NULL,
                                detectionSize = 1e6,
                                detectionDrivers = 99,
                                thresholdWhole = x)
        nca$popSummary[1:5, c(1:3, 8:9)]
        nca2$popSummary[1:5, c(1:3, 8:9)]
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
        T1 <- ( t.test(mutsPerClone1 ,
                     mutsPerClone2, alternative = "greater")$p.value < p.value.threshold)
        T2 <- ( wilcox.test(nca$popSummary[, "NumClones"] ,
                     nca2$popSummary[, "NumClones"], alternative = "greater")$p.value < p.value.threshold)
        ## Pop sizes here do not really differ, as we start with the
        ## initMutant with increased s
        summary(nca$popSummary[, "TotalPopSize"])
        summary(nca2$popSummary[, "TotalPopSize"])
        
        if(T1 && T2 && T3 && T4 && T5 && T6 && T7 && T8) break;
    }
    cat(paste("\n done tries", tries, "\n"))
    expect_true(T1 && T2 && T3 && T4 && T5 && T6 && T7 && T8)
})
date()





date()
test_that("mutPropGrowth diffs with s> 0, oncoSimulSample", {
    max.tries <- 4
    for(tries in 1:max.tries) {
        T1 <- T2 <- T3 <- T4 <- T5 <- T6 <- T7 <- T8 <- TRUE
        cat("\n oss1_ontime: a runif is", runif(1), "\n")
        
        ft <- 4 ## we stop on time 
        pops <- 50
        lni <- 200 ## 150
        no <- 1e3 ## 5e1 
        ni <- c(2, rep(0, lni)) ## 2 ## 4 ## 5
        mu <- 5e-7 ## 1e-6
        x <- 1e-9
        names(ni) <- c("a", paste0("n", seq.int(lni)))
        ni <- sample(ni) ## scramble
        fe <- allFitnessEffects(noIntGenes = ni)
        cat("\n oss1a: a runif is", runif(1), "\n")
        nca <- oncoSimulSample(pops, fe, finalTime = ft, detectionProb = NA,
                               mu = mu,
                               mutationPropGrowth = TRUE,
                               initSize = no, sampleEvery = 0.02,
                               initMutant = "a", 
                               onlyCancer = FALSE, seed = NULL,
                               detectionSize = 1e9,
                               detectionDrivers = 99,
                               thresholdWhole = x)
        cat("\n oss1c: a runif is", runif(1), "\n")
        nca2 <- oncoSimulSample(pops, fe, finalTime = ft, detectionProb = NA,
                                mu = mu,
                                mutationPropGrowth = FALSE,
                                initSize = no, sampleEvery = 0.02,
                                initMutant = "a", 
                                onlyCancer = FALSE, seed = NULL,
                                detectionSize = 1e9,
                                detectionDrivers = 99,
                                thresholdWhole = x)
        nca$popSummary[1:5, c(1:3, 8:9)]
        nca2$popSummary[1:5, c(1:3, 8:9)]
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
        T1 <- ( t.test(mutsPerClone1 ,
                     mutsPerClone2, alternative = "greater")$p.value < p.value.threshold)
        T2 <- ( wilcox.test(nca$popSummary[, "NumClones"] ,
                     nca2$popSummary[, "NumClones"], alternative = "greater")$p.value < p.value.threshold)
        ## Pop sizes here do not really differ, as we start with the
        ## initMutant with increased s
        summary(nca$popSummary[, "TotalPopSize"])
        summary(nca2$popSummary[, "TotalPopSize"])
        
        if(T1 && T2 && T3 && T4 && T5 && T6 && T7 && T8) break;
    }
    cat(paste("\n done tries", tries, "\n"))
    expect_true(T1 && T2 && T3 && T4 && T5 && T6 && T7 && T8)
})
date()



date()
test_that("mutPropGrowth diffs with s> 0, oncoSimulSample, McFL", {
    max.tries <- 4
    for(tries in 1:max.tries) {
        
        T1 <- T2 <- T3 <- T4 <- T5 <- T6 <- T7 <- T8 <- TRUE
        cat("\n ossmcf1: a runif is", runif(1), "\n")
        ft <- 100  ## we stop on size
        pops <- 50
        lni <- 200 ## 150
        no <- 1e3 ## 5e1 
        ni <- c(3, rep(0, lni)) ## 2 ## 4 ## 5
        mu <- 5e-7 ## 1e-6
        x <- 1e-9
        names(ni) <- c("a", paste0("n", seq.int(lni)))
        ni <- sample(ni) ## scramble
        fe <- allFitnessEffects(noIntGenes = ni)
        cat("\n ossmcf1a: a runif is", runif(1), "\n")
        nca <- oncoSimulSample(pops, fe, finalTime = ft, detectionProb = NA,
                               mu = mu, model = "McFL",
                               mutationPropGrowth = TRUE,
                               initSize = no, sampleEvery = 0.01,
                               initMutant = "a", 
                               onlyCancer = FALSE, seed = NULL,
                               detectionSize = 3e4,
                               detectionDrivers = 99,
                               thresholdWhole = x)
        cat("\n ossmcf1c: a runif is", runif(1), "\n")
        nca2 <- oncoSimulSample(pops, fe, finalTime = ft, detectionProb = NA,
                                mu = mu, model = "McFL",
                                mutationPropGrowth = FALSE,
                                initSize = no, sampleEvery = 0.01,
                                initMutant = "a", 
                                onlyCancer = FALSE, seed = NULL,
                                detectionSize = 3e4,
                                detectionDrivers = 99,
                                thresholdWhole = x)
        nca$popSummary[1:5, c(1:3, 8:9)]
        ## summary(nca$popSummary[, "NumClones"])
        ## summary(nca$popSummary[, "TotalPopSize"])
        nca2$popSummary[1:5, c(1:3, 8:9)]
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
        T1 <- ( t.test(mutsPerClone1 ,
                     mutsPerClone2, alternative = "greater")$p.value < p.value.threshold)
        T2 <- ( wilcox.test(nca$popSummary[, "NumClones"] ,
                     nca2$popSummary[, "NumClones"], alternative = "greater")$p.value < p.value.threshold)
        ## Pop sizes here do not really differ, as we start with the
        ## initMutant with increased s. And this is McFL, so bounded from above.
        summary(nca$popSummary[, "TotalPopSize"])
        summary(nca2$popSummary[, "TotalPopSize"])
        
        if(T1 && T2 && T3 && T4 && T5 && T6 && T7 && T8) break;
    }
    cat(paste("\n done tries", tries, "\n"))
    expect_true(T1 && T2 && T3 && T4 && T5 && T6 && T7 && T8)
})
date()




date()
test_that("mutPropGrowth diffs with s> 0, oncoSimulSample, McFL, stop on time", {
    max.tries <- 4
    for(tries in 1:max.tries) {
        
        T1 <- T2 <- T3 <- T4 <- T5 <- T6 <- T7 <- T8 <- TRUE
        cat("\n ossmcf1: a runif is", runif(1), "\n")
        ft <- 10 ## we stop on time but sizes are about same. We are in the plateau
        pops <- 50
        lni <- 200 ## 150
        no <- 1e3 ## 5e1 
        ni <- c(3, rep(0, lni)) ## 2 ## 4 ## 5
        mu <- 5e-7 ## 1e-6
        x <- 1e-9
        names(ni) <- c("a", paste0("n", seq.int(lni)))
        ni <- sample(ni) ## scramble
        fe <- allFitnessEffects(noIntGenes = ni)
        cat("\n ossmcf1a: a runif is", runif(1), "\n")
        nca <- oncoSimulSample(pops, fe, finalTime = ft, detectionProb = NA,
                               mu = mu, model = "McFL",
                               mutationPropGrowth = TRUE,
                               initSize = no, sampleEvery = 0.01,
                               initMutant = "a", 
                               onlyCancer = FALSE, seed = NULL,
                               detectionSize = 3e7,
                               detectionDrivers = 99,
                               thresholdWhole = x)
        cat("\n ossmcf1c: a runif is", runif(1), "\n")
        nca2 <- oncoSimulSample(pops, fe, finalTime = ft, detectionProb = NA,
                                mu = mu, model = "McFL",
                                mutationPropGrowth = FALSE,
                                initSize = no, sampleEvery = 0.01,
                                initMutant = "a", 
                                onlyCancer = FALSE, seed = NULL,
                                detectionSize = 3e7,
                                detectionDrivers = 99,
                                thresholdWhole = x)
        nca$popSummary[1:5, c(1:3, 8:9)]
        ## summary(nca$popSummary[, "NumClones"])
        ## summary(nca$popSummary[, "TotalPopSize"])
        nca2$popSummary[1:5, c(1:3, 8:9)]
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
        T1 <- ( t.test(mutsPerClone1 ,
                     mutsPerClone2, alternative = "greater")$p.value < p.value.threshold)
        T2 <- ( wilcox.test(nca$popSummary[, "NumClones"] ,
                     nca2$popSummary[, "NumClones"], alternative = "greater")$p.value < p.value.threshold)
        ## Pop sizes here do not really differ, as we start with the
        ## initMutant with increased s. And this is McFL, so bounded from above.
        summary(nca$popSummary[, "TotalPopSize"])
        summary(nca2$popSummary[, "TotalPopSize"])
        
        if(T1 && T2 && T3 && T4 && T5 && T6 && T7 && T8) break;
    }
    cat(paste("\n done tries", tries, "\n"))
    expect_true(T1 && T2 && T3 && T4 && T5 && T6 && T7 && T8)
})
date()






test_that("We crash as we should", {
    so <- 0.2
    feo <- allFitnessEffects(orderEffects = c("a > b" = so,
                                              "b > a" = -5),
                             noIntGenes = rep(0, 50))
    ## with testthat 0.11.0.9000 we should be able
    ## to use use_catch to catch the C++ exception directly
    o1 <- oncoSimulIndiv(feo,
               mu = 1e-7, detectionProb = NA,
               initSize = 1e4,
               K = 1e4,
               model = "Exp",
               detectionDrivers = 99,
               finalTime = 0.01,
               detectionSize = 1e12,
               sampleEvery = 0.001,
               keepEvery = 30,
               initMutant = "b > a",
               mutationPropGrowth = TRUE, 
               onlyCancer = FALSE,
               verbosity = 0)
    expect_true(grepl("pE not finite", o1$other$ExceptionMessage,
                      fixed = TRUE))
})

test_that("...and we don't when we shouldn't", {
    so <- 0.2
    feo <- allFitnessEffects(orderEffects = c("a > b" = so,
                                              "b > a" = -5),
                             noIntGenes = rep(0, 50))
    ## with testthat 0.11.0.9000 we should be able
    ## to use use_catch to catch the C++ exception directly
    o1 <- oncoSimulIndiv(feo,
               mu = 1e-7,
               initSize = 1e4,
               K = 1e4,
               model = "Exp",
               detectionDrivers = 99,
               finalTime = 0.01,
               detectionSize = 1e12,
               sampleEvery = 0.001,
               keepEvery = 30,
               initMutant = "b > a",
               mutationPropGrowth = FALSE, 
               onlyCancer = FALSE,
               verbosity = 0)
    expect_true(length(grepl("pE not finite", o1$other$ExceptionMessage,
                      fixed = TRUE)) == 0)
})




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



