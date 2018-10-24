cat(paste("\n Starting at mutPropGrowth long", date(), "\n"))

## Note that many of the tests below where we do a test and have a
## comparison like "whatever$p.value > p.fail " are, of course, expected
## to fail with prob. ~ p.fail even if things are perfectly OK.


## Why this does not really reflect what we want, and why number of clones
## is better that capture the idea of "more mutations". NumClones reflects
## the creation of a new clone, something that happens whenever there is a
## mutation (that does not land you on a pre-existing clone).

######################################################################
######################################################################

## ## The functions below measure number of mutated positions by summing
## ## number of alleles and dividing by pop.size. But only at final time.
## ## And large pops of clones with few muts remain and swamp.

## popS <- function(out) unlist(lapply(out, function(x) x$TotalPopSize))

## muts <- function(out) {
##     popSize <- popS(out)
##     gc <- rowSums(OncoSimulR:::geneCounts(out))
##     Muts.per.indiv <- na.omit(gc/popSize)
##     return(list(PopSize = popSize,
##                 Muts = gc,
##                 Muts.per.indiv = Muts.per.indiv,
##                 Muts.per.indiv.no.0 = Muts.per.indiv[gc > 0]))
## }


## This is better, but still you need time to allow accumulation of clones
## with many mutations, and those with few remain

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


######################################################################
######################################################################




## RNGkind("L'Ecuyer-CMRG") ## for the mclapplies
## If crashes I want to see where: thus output seed.

p.value.threshold <- 0.001

## this tests takes 10 seconds. Moved to long. So this was in the standard
## ones.
date()
test_that("mutPropGrowth diffs with s> 0", {
    max.tries <- 4
    for(tries in 1:max.tries) {
TTT <- NULL
    
    cat("\n mgp1: a runif is", runif(1), "\n")
    ft <- 100 ## 4  we are stopping on size
    pops <- 50
    lni <- 150 ## 100
    no <- 1e3 ## 5e1 
    ni <- c(2, rep(0, lni)) ## 2 ## 4 ## 5
    mu <- 1e-5 ## 1e-6
    names(ni) <- c("a", paste0("n", seq.int(lni)))
    ni <- sample(ni) ## scramble 
    fe <- allFitnessEffects(noIntGenes = ni)
    cat("\n mpg1a: a runif is", runif(1), "\n")
    nca <- oncoSimulPop(pops, fe, finalTime = ft,
                        mu = mu, detectionSize =  3e4,
                        mutationPropGrowth = TRUE,
                        initSize = no, sampleEvery = 0.1,
                        initMutant = "a", keepEvery = 1,
                        onlyCancer = FALSE, detectionProb = NA, seed = NULL, mc.cores = 2)
    cat("\n mpg1c: a runif is", runif(1), "\n")
    nca2 <- oncoSimulPop(pops, fe, finalTime = ft,
                         mu = mu, detectionSize =  3e4,
                        mutationPropGrowth = FALSE,
                        initSize = no, sampleEvery = 0.1,
                        initMutant = "a", keepEvery = 1,
                        onlyCancer = FALSE, detectionProb = NA, seed = NULL, mc.cores = 2)
    summary(nca)[1:20, c(1, 2, 3, 8, 9)]
    summary(nca2)[1:20, c(1, 2, 3, 8, 9)]
    ## I once saw a weird thing
    ## summary(summary(nca)$NumClones)
    ## summary(summary(nca2)$NumClones)
    ## summary(mutsPerClone(nca))
    ## summary(mutsPerClone(nca2))
    ## The real comparison
    TTT <- c(TTT,  wilcox.test(summary(nca)$NumClones,
                             summary(nca2)$NumClones,
                             alternative = "greater")$p.value < p.value.threshold)
    TTT <- c(TTT,  wilcox.test(mutsPerClone(nca),
                             mutsPerClone(nca2),
                             alternative = "greater")$p.value < p.value.threshold)

if( all(TTT) ) break;
    }
    cat(paste("\n done tries", tries, "\n"))
    expect_true(all(TTT))

    })
cat("\n", date(), "\n")






date()
test_that("mutPropGrowth diffs with s> 0, stopping on time", {
    max.tries <- 4
    for(tries in 1:max.tries) {
TTT <- NULL
    
    cat("\n mgp1,ontime : a runif is", runif(1), "\n")
    ft <- 2 ## 4  we are stopping on time; sizes are similar
    pops <- 50
    lni <- 150 ## 100
    no <- 1e3 ## 5e1 
    ni <- c(2, rep(0, lni)) ## 2 ## 4 ## 5
    mu <- 1e-5 ## 1e-6
    names(ni) <- c("a", paste0("n", seq.int(lni)))
    ni <- sample(ni) ## scramble 
    fe <- allFitnessEffects(noIntGenes = ni)
    cat("\n mpg1a: a runif is", runif(1), "\n")
    nca <- oncoSimulPop(pops, fe, finalTime = ft,
                        mu = mu, detectionSize =  3e9,
                        mutationPropGrowth = TRUE,
                        initSize = no, sampleEvery = 0.1,
                        initMutant = "a", keepEvery = 1,
                        onlyCancer = FALSE, detectionProb = NA, seed = NULL, mc.cores = 2)
    cat("\n mpg1c: a runif is", runif(1), "\n")
    nca2 <- oncoSimulPop(pops, fe, finalTime = ft,
                         mu = mu, detectionSize =  3e9,
                        mutationPropGrowth = FALSE,
                        initSize = no, sampleEvery = 0.1,
                        initMutant = "a", keepEvery = 1,
                        onlyCancer = FALSE, detectionProb = NA, seed = NULL, mc.cores = 2)
    summary(nca)[1:20, c(1, 2, 3, 8, 9)]
    summary(nca2)[1:20, c(1, 2, 3, 8, 9)]
    ## I once saw a weird thing
    ## summary(summary(nca)$NumClones)
    ## summary(summary(nca2)$NumClones)
    ## summary(mutsPerClone(nca))
    ## summary(mutsPerClone(nca2))
    ## The real comparison
    TTT <- c(TTT,  wilcox.test(summary(nca)$NumClones,
                             summary(nca2)$NumClones,
                             alternative = "greater")$p.value < p.value.threshold)
    TTT <- c(TTT,  wilcox.test(mutsPerClone(nca),
                             mutsPerClone(nca2),
                             alternative = "greater")$p.value < p.value.threshold)

if( all(TTT) ) break;
    }
    cat(paste("\n done tries", tries, "\n"))
    expect_true(all(TTT))

    })
cat("\n", date(), "\n")






### Now, test mutPropGrowth does not lead to differences when there are no differences in growth.
## Same settings as other tests with similar names, but now no fitness, so
## no growth. So no effect of mutationPropGrowth. Note that if we carry
## out the tests below when there are differences as in previous tests, we
## p <<< 1e-6.

## Note that these tests will fail, as they are based on p-values, with
## ... well, a frequency given approx. by the p-value. Thus, we do not use
## them in the standard tests. But I run them as part of the long tests;
## if they fail, run again to make sure it is just that an unlikely event
## happened.


date()
test_that("mutPropGrowth no diffs with s = 0", {
      max.tries <- 4
    for(tries in 1:max.tries) {
TTT <- NULL
   
    
    cat("\n mgp1ND: a runif is", runif(1), "\n")
    ft <- 4 ## There is no growth expected, so stopping on time is
            ## fine. Actually, the only  reasonable thing
    pops <- 100
    lni <- 150 ## 100
    no <- 1e3 ## 5e1 
    ni <- c(0, rep(0, lni)) ## 2 ## 4 ## 5
    mu <- 1e-5 ## 1e-6
    names(ni) <- c("a", paste0("n", seq.int(lni)))
    ni <- sample(ni) ## scramble 
    fe <- allFitnessEffects(noIntGenes = ni)
    cat("\n mpg1NDa: a runif is", runif(1), "\n")
    nca <- oncoSimulPop(pops, fe, finalTime = ft,
                        mu = mu,
                        mutationPropGrowth = TRUE,
                        initSize = no, sampleEvery = 0.1,
                        initMutant = "a", keepEvery = 1,
                        onlyCancer = FALSE, detectionProb = NA, seed = NULL, mc.cores = 2)
    cat("\n mpg1NDc: a runif is", runif(1), "\n")
    nca2 <- oncoSimulPop(pops, fe, finalTime = ft,
                         mu = mu,
                        mutationPropGrowth = FALSE,
                        initSize = no, sampleEvery = 0.1,
                        initMutant = "a", keepEvery = 1,
                        onlyCancer = FALSE, detectionProb = NA, seed = NULL, mc.cores = 2)
    summary(nca)[1:20, c(1, 2, 3, 8, 9)]
    summary(nca2)[1:20, c(1, 2, 3, 8, 9)]
    p.fail <- 1e-3
    TTT <- c(TTT, t.test(sqrt(summary(nca)$NumClones),
                       sqrt(summary(nca2)$NumClones))$p.value > p.fail)
    TTT <- c(TTT, t.test(sqrt(mutsPerClone(nca)),
                       sqrt(mutsPerClone(nca2)))$p.value > p.fail)
    ## I once saw a weird thing
    ## expect_true(var(summary(nca)$NumClones) > 1e-4)
    ## expect_true(var(summary(nca2)$NumClones) > 1e-4)
    ## summary(summary(nca)$NumClones)
    ## summary(summary(nca2)$NumClones)
    ## summary(mutsPerClone(nca))
    ## summary(mutsPerClone(nca2))

    if( all(TTT) ) break;
    }
    cat(paste("\n done tries", tries, "\n"))
    expect_true(all(TTT))
})
cat("\n", date(), "\n")















date()
test_that("mutPropGrowth no diffs with s = 0, McFL", {
    
     max.tries <- 4
    for(tries in 1:max.tries) {
TTT <- NULL
    
    cat("\n mcfND1: a runif is", runif(1), "\n")
    ft <- 3 
    pops <- 200
    lni <- 100
    no <- 1e3 ## 5e1 
    ni <- c(0, rep(0, lni)) ## 5
    names(ni) <- c("a", paste0("n", seq.int(lni)))
    ni <- sample(ni) ## scramble
    fe <- allFitnessEffects(noIntGenes = ni)
    
    
    cat("\n mcfND1a: a runif is", runif(1), "\n")
    nca <- oncoSimulPop(pops, fe, finalTime = ft,
                        mutationPropGrowth = TRUE,
                        initSize = no,
                        initMutant = "a", model = "McFL",
                        onlyCancer = FALSE, detectionProb = NA, seed = NULL, mc.cores = 2)
    
    
    cat("\n mcfND1c: a runif is", runif(1), "\n")
    nca2 <- oncoSimulPop(pops, fe, finalTime = ft,
                        mutationPropGrowth = FALSE,
                        initSize = no,
                        initMutant = "a", model = "McFL",
                        onlyCancer = FALSE, detectionProb = NA, seed = NULL, mc.cores = 2)
    ## summary(nca)[1:20, c(1, 2, 3, 8, 9)]
    ## summary(nca2)[1:20, c(1, 2, 3, 8, 9)]
    ## summary(summary(nca)$NumClones)
    ## summary(summary(nca2)$NumClones)
    ## summary(mutsPerClone(nca))
    ## summary(mutsPerClone(nca2))
    p.fail <- 1e-3
    TTT <- c(TTT, t.test(sqrt(summary(nca)$NumClones),
                       sqrt(summary(nca2)$NumClones))$p.value > p.fail)
    TTT <- c(TTT, t.test(sqrt(mutsPerClone(nca)),
                       sqrt(mutsPerClone(nca2)))$p.value > p.fail)
    if( all(TTT) ) break;
    }
    cat(paste("\n done tries", tries, "\n"))
    expect_true(all(TTT))
})
cat("\n", date(), "\n")




date()
test_that("mutPropGrowth no diffs with s = 0, oncoSimulSample", {
    
     max.tries <- 4
    for(tries in 1:max.tries) {
TTT <- NULL
    
    cat("\n oss1: a runif is", runif(1), "\n")
    ft <- 3.5 
    pops <- 250
    lni <- 200
    no <- 1e3 
    ni <- c(0, rep(0, lni)) 
    mu <- 5e-7 
    x <- 1e-9
    names(ni) <- c("a", paste0("n", seq.int(lni)))
    ni <- sample(ni) ## scramble
    fe <- allFitnessEffects(noIntGenes = ni)
    
    
    cat("\n oss1a: a runif is", runif(1), "\n")
    nca <- oncoSimulSample(pops, fe, finalTime = ft,
                        mu = mu,
                        mutationPropGrowth = TRUE,
                        initSize = no, sampleEvery = 0.1,
                        initMutant = "a", 
                        onlyCancer = FALSE, detectionProb = NA, seed = NULL,
                        detectionSize = 1e9,
                        detectionDrivers = 99,
                        thresholdWhole = x)
    
    
    cat("\n oss1c: a runif is", runif(1), "\n")
    nca2 <- oncoSimulSample(pops, fe, finalTime = ft,
                         mu = mu,
                        mutationPropGrowth = FALSE,
                        initSize = no, sampleEvery = 0.1,
                        initMutant = "a", 
                        onlyCancer = FALSE, detectionProb = NA, seed = NULL,
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
    p.fail <- 1e-3
    TTT <- c(TTT, t.test(sqrt( nca$popSummary[, "NumClones"]  ),
                       sqrt( nca2$popSummary[, "NumClones"]  ))$p.value > p.fail)
    TTT <- c(TTT, t.test(sqrt(mutsPerClone1),
                       sqrt(mutsPerClone2))$p.value > p.fail)
    if( all(TTT) ) break;
    }
    cat(paste("\n done tries", tries, "\n"))
    expect_true(all(TTT))
})


date()
test_that("mutPropGrowth no diffs with s = 0, oncoSimulSample, McFL", {
    
      max.tries <- 4
    for(tries in 1:max.tries) {
TTT <- NULL
   
    cat("\n ossmcfND1: a runif is", runif(1), "\n")
    ft <- 40 ## 4
    pops <- 250
    lni <- 200 ## 150
    no <- 1e3 ## 5e1 
    ni <- c(0, rep(0, lni)) ## 2 ## 4 ## 5
    mu <- 5e-7 ## 1e-6
    x <- 1e-9
    names(ni) <- c("a", paste0("n", seq.int(lni)))
    ni <- sample(ni) ## scramble
    fe <- allFitnessEffects(noIntGenes = ni)
    
    
    cat("\n ossmcfND1a: a runif is", runif(1), "\n")
    nca <- oncoSimulSample(pops, fe, finalTime = ft,
                        mu = mu, model = "McFL",
                        mutationPropGrowth = TRUE,
                        initSize = no, sampleEvery = 0.01,
                        initMutant = "a", 
                        onlyCancer = FALSE, detectionProb = NA, seed = NULL,
                        detectionSize = 1e9,
                        detectionDrivers = 99,
                        thresholdWhole = x)
    
    
    cat("\n ossmcfND1c: a runif is", runif(1), "\n")
    nca2 <- oncoSimulSample(pops, fe, finalTime = ft,
                         mu = mu, model = "McFL",
                        mutationPropGrowth = FALSE,
                        initSize = no, sampleEvery = 0.01,
                        initMutant = "a", 
                        onlyCancer = FALSE, detectionProb = NA, seed = NULL,
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
    p.fail <- 1e-3
    TTT <- c(TTT, t.test(sqrt( nca$popSummary[, "NumClones"]  ),
                       sqrt( nca2$popSummary[, "NumClones"]  ))$p.value > p.fail)
    TTT <- c(TTT, t.test(sqrt(mutsPerClone1),
                       sqrt(mutsPerClone2))$p.value > p.fail)
    if( all(TTT) ) break;
    }
    cat(paste("\n done tries", tries, "\n"))
    expect_true(all(TTT))
})
cat("\n", date(), "\n")


## The tests below can occasionally fail (but that probability decreases
## as we increase number of pops), as they should.
## Some of these tests take some time. For now, show times.

## I got obfuscated in many tests below. They are left as they have some
## interesting ideas for playing, but are not good for testing because
## they take forever.


## Beware of exiting because max number of subjects reached, and that
## could happen sooner for the faster growing, so less time to accumulate
## mutations. Likewise, differences between nca and nca2 depend on large
## enough differences in mutation, and thus a relatively large multiplier
## factor, so a large s.





## Idea is that the mutation will be to the module that increases growth
## rate. And in one case that will lead to increase in mutation rate. But
## we have increased variability in outcome because when that mutant hits
## is, well, stochastic.

cat("\n", date(), "\n")
test_that("oncoSimulSample Without initmutant and modules", {
    
     max.tries <- 4
     for(tries in 1:max.tries) {
         
         TTT <- NULL
         cat("\n osS: a runif is", runif(1), "\n")
    pops <- 70
    lni <- 1 ## no fitness effects genes
    fni <- 50 ## fitness effects genes
    no <- 1e4 
    ft <- 5 ## we stop on time; sizes are roughly similar
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
    cat("\n osSa: a runif is", runif(1), "\n")
    b1 <- oncoSimulSample(pops,
                          f3,
                          mu = mu,
                          mutationPropGrowth = FALSE,
                          finalTime =ft,
                          initSize = no,
                          onlyCancer = FALSE, detectionProb = NA,
                          sampleEvery = 0.01,
                          detectionSize = 3e9,
                          detectionDrivers = 99,
                          seed =NULL,
                          thresholdWhole = x)
    cat("\n osSb: a runif is", runif(1), "\n")
    b2 <- oncoSimulSample(pops,
                         f3,
                         mu = mu,
                         mutationPropGrowth = TRUE,
                         finalTime =ft,
                         initSize = no,
                         onlyCancer = FALSE, detectionProb = NA,
                         sampleEvery = 0.01,
                          detectionSize = 3e9,
                          detectionDrivers = 99,
                          seed =NULL,
                         thresholdWhole = x)
         b2$popSummary[1:5, c(1:3, 8:9)]
         b1$popSummary[1:5, c(1:3, 8:9)]
         summary(b2$popSummary[, "NumClones"])
         summary(b1$popSummary[, "NumClones"])
         summary(b2$popSummary[, "TotalPopSize"])
         summary(b1$popSummary[, "TotalPopSize"])
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
    TTT <- c(TTT,  wilcox.test(mutsPerClone2,
                             mutsPerClone1,
                             alternative = "greater")$p.value < p.value.threshold)
    TTT <- c(TTT,  wilcox.test(b2$popSummary[, "NumClones"],
                             b1$popSummary[, "NumClones"],
                             alternative = "greater")$p.value < p.value.threshold)

if( all(TTT) ) break;
    }
    cat(paste("\n done tries", tries, "\n"))
    expect_true(all(TTT))
})
cat("\n", date(), "\n")


cat("\n", date(), "\n")
test_that("oncoSimulSample Without initmutant and modules, McFL", {
      max.tries <- 4
    for(tries in 1:max.tries) {

        TTT <- NULL
    cat("\n osSMcFL: a runif is", runif(1), "\n")
    pops <- 90
    lni <- 1 ## no fitness effects genes
    fni <- 50 ## fitness effects genes
    no <- 1e4 ## note we use only 10 in the other example below
    ft <- 10 ## stopping on time and sizes are roughly similar 
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
    cat("\n osSMcFLa: a runif is", runif(1), "\n")
    b1 <- oncoSimulSample(pops,
                          f3,
                          mu = mu,
                          mutationPropGrowth = FALSE,
                          finalTime =ft,
                          initSize = no,
                          onlyCancer = FALSE, detectionProb = NA,
                          sampleEvery = 0.01,
                          detectionSize = 5e9,
                          detectionDrivers = 99,
                          seed =NULL,
                          model = "McFL",
                          thresholdWhole = x)
    cat("\n osSMcFLb: a runif is", runif(1), "\n")
    b2 <- oncoSimulSample(pops,
                         f3,
                         mu = mu,
                         mutationPropGrowth = TRUE,
                         finalTime =ft,
                         initSize = no,
                         onlyCancer = FALSE, detectionProb = NA,
                         sampleEvery = 0.01,
                          detectionSize = 5e9,
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
    TTT <- c(TTT,  wilcox.test(mutsPerClone2,
                             mutsPerClone1,
                             alternative = "greater")$p.value < p.value.threshold)
    TTT <- c(TTT,  wilcox.test(b2$popSummary[, "NumClones"],
                             b1$popSummary[, "NumClones"],
                             alternative = "greater")$p.value < p.value.threshold)
        
if( all(TTT) ) break;
    }
    cat(paste("\n done tries", tries, "\n"))
    expect_true(all(TTT))

})
cat("\n", date(), "\n")



## Same as above, but we stop in detectionSize, so no differences are
## attributable just to differences in population size.


cat("\n", date(), "\n")
test_that("oncoSimulSample Without initmutant and modules, fixed size", {
    max.tries <- 4
    for(tries in 1:max.tries) {

        TTT <- NULL
    cat("\n osSFPS: a runif is", runif(1), "\n")
    pops <- 90
    lni <- 1 ## no fitness effects genes
    fni <- 50 ## fitness effects genes
    no <- 1e4 
    ft <- 900 ## stopping on size
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
    cat("\n osSFPSa: a runif is", runif(1), "\n")
    b1 <- oncoSimulSample(pops,
                          f3,
                          mu = mu,
                          mutationPropGrowth = FALSE,
                          finalTime =ft,
                          initSize = no,
                          onlyCancer = FALSE, detectionProb = NA,
                          sampleEvery = 0.01,
                          detectionSize = 6e4,
                          detectionDrivers = 99,
                          seed =NULL,
                          thresholdWhole = x)
    cat("\n osSFPSb: a runif is", runif(1), "\n")
    b2 <- oncoSimulSample(pops,
                         f3,
                         mu = mu,
                         mutationPropGrowth = TRUE,
                         finalTime =ft,
                         initSize = no,
                         onlyCancer = FALSE, detectionProb = NA,
                         sampleEvery = 0.01,
                          detectionSize = 6e4,
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
    TTT <- c(TTT,  wilcox.test(mutsPerClone2,
                             mutsPerClone1,
                             alternative = "greater")$p.value < p.value.threshold)
    TTT <- c(TTT,  wilcox.test(b2$popSummary[, "NumClones"],
                             b1$popSummary[, "NumClones"],
                             alternative = "greater")$p.value < p.value.threshold)

if( all(TTT) ) break;
    }
    cat(paste("\n done tries", tries, "\n"))
    expect_true(all(TTT))
})
cat("\n", date(), "\n")


cat("\n", date(), "\n")
test_that("oncoSimulSample Without initmutant and modules, McFL, fixed size", {
    
     max.tries <- 4
     for(tries in 1:max.tries) {

         
         TTT <- NULL
    cat("\n osSFPSMcFL: a runif is", runif(1), "\n")
    pops <- 70
    lni <- 1 ## no fitness effects genes
    fni <- 50 ## fitness effects genes
    no <- 1e4 ## note we use only 10 in the other example below
    ft <- 1000 ## stop on size
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
    cat("\n osSFPSMcFLa: a runif is", runif(1), "\n")
    b1 <- oncoSimulSample(pops,
                          f3,
                          mu = mu,
                          mutationPropGrowth = FALSE,
                          finalTime =ft,
                          initSize = no,
                          onlyCancer = FALSE, detectionProb = NA,
                          sampleEvery = 0.01,
                          detectionSize = 5e4,
                          detectionDrivers = 99,
                          seed =NULL,
                          model = "McFL",
                          thresholdWhole = x)
    cat("\n osSFPSMcFLb: a runif is", runif(1), "\n")
    b2 <- oncoSimulSample(pops,
                         f3,
                         mu = mu,
                         mutationPropGrowth = TRUE,
                         finalTime =ft,
                         initSize = no,
                         onlyCancer = FALSE, detectionProb = NA,
                         sampleEvery = 0.01,
                          detectionSize = 5e4,
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
    TTT <- c(TTT,  wilcox.test(mutsPerClone2,
                             mutsPerClone1,
                             alternative = "greater")$p.value < p.value.threshold)
    TTT <- c(TTT,  wilcox.test(b2$popSummary[, "NumClones"],
                             b1$popSummary[, "NumClones"],
                             alternative = "greater")$p.value < p.value.threshold)



         if( all(TTT) ) break;
    }
    cat(paste("\n done tries", tries, "\n"))
    expect_true(all(TTT))

})
cat("\n", date(), "\n")






## Psss ideas.  The idea here is that you hit the module with fitness
## effects, that sets things to grow, and then you either have or not
## mutation proportional to growth.

cat("\n", date(), "\n")
test_that("Without initmutant", {
    
     max.tries <- 4
     for(tries in 1:max.tries) {
         
         TTT <- NULL
         cat("\n s3: a runif is", runif(1), "\n")
    pops <- 50
    lni <- 1 ## no fitness effects genes
    fni <- 50 ## fitness effects genes
    no <- 1e3
ft <- 5  ## we stop on time; below we repeat stopping on size
ds <- 1e9
    s3 <- 3.0
    mu <- 5e-5 ## easier to see
    ## noInt have no fitness effects, but can accumulate mutations
    ni <- rep(0, lni)
    ## Those with fitness effects in one module, so
    ## neither fitness nor mut. rate blow up
    gn <- paste(paste0("a", 1:fni), collapse = ", ")
    f3 <- allFitnessEffects(epistasis = c("A" = s3),
                            geneToModule = c("A" = gn),
                            noIntGenes = ni)
    cat("\n s3a: a runif is", runif(1), "\n")
    s3.ng <- oncoSimulPop(pops,
                          f3,
                          mu = mu,
                          mutationPropGrowth = FALSE, sampleEvery = 0.01,
                          finalTime =ft, detectionSize = ds,
                          initSize = no, detectionDrivers = 99999,
                          onlyCancer = FALSE, detectionProb = NA,
                          seed = NULL, mc.cores = 2)
    cat("\n s3b: a runif is", runif(1), "\n")
    s3.g <- oncoSimulPop(pops,
                         f3,
                         mu = mu,
                         mutationPropGrowth = TRUE, sampleEvery = 0.01,
                         finalTime =ft, detectionSize = ds,
                         initSize = no,, detectionDrivers = 99999,
                         onlyCancer = FALSE, detectionProb = NA,
                         seed = NULL, mc.cores = 2)
    summary(s3.g)[, c(1, 2, 3, 8, 9)]
    summary(s3.ng)[, c(1, 2, 3, 8, 9)]
    TTT <- c(TTT,  wilcox.test(mutsPerClone(s3.g),
                             mutsPerClone(s3.ng),
                             alternative = "greater")$p.value < p.value.threshold)
    TTT <- c(TTT,  wilcox.test(summary(s3.g)$NumClones,
                             summary(s3.ng)$NumClones,
                             alternative = "greater")$p.value < p.value.threshold)
gc() 

if( all(TTT) ) break;
    }
    cat(paste("\n done tries", tries, "\n"))
    expect_true(all(TTT))
})
cat("\n", date(), "\n")

cat("\n", date(), "\n")
test_that("Without initmutant, 2", {
    
    ## More of the above. Use smaller s2 and smaller mutation, but then to
    ## see it reliably you need large ft and we also increase
    ## init. pop. size. Variability is huge specially because we stop on time.
    max.tries <- 4
    for(tries in 1:max.tries) {

        TTT <- NULL
        cat("\n s2: a runif is", runif(1), "\n")
        s2 <- 1.0
        ft <- 6 
        ds <- 1e9
        pops <- 200
        lni <- 1 ## no fitness effects genes
        fni <- 200 ## fitness effects genes
        no <- 1e4
        mu <- 5e-6 ## easier to see
        ## noInt have no fitness effects, but can accumulate mutations
        ni <- rep(0, lni)
        ## Those with fitness effects in one module, so
        ## neither fitness nor mut. rate blow up
        gn <- paste(paste0("a", 1:fni), collapse = ", ")
        f2 <- allFitnessEffects(epistasis = c("A" = s2),
                                geneToModule = c("A" = gn),
                                noIntGenes = ni)
        s2.ng <- oncoSimulPop(pops,
                          f2,
                          mu = mu,
                          mutationPropGrowth = FALSE, sampleEvery = 0.01,
                          finalTime =ft, detectionSize = ds,
                          initSize = no, detectionDrivers = 99999,
                          onlyCancer = FALSE, detectionProb = NA,
                          seed = NULL, mc.cores = 2)
        cat("\n s2b: a runif is", runif(1), "\n")
        s2.g <- oncoSimulPop(pops,
                         f2,
                         mu = mu,
                         mutationPropGrowth = TRUE, sampleEvery = 0.01,
                         finalTime =ft, detectionSize = ds,
                         initSize = no, detectionDrivers = 99999,
                         onlyCancer = FALSE, detectionProb = NA,
                         seed = NULL, mc.cores = 2)
        if(! (inherits(s2.ng, "oncosimulpop") &&
              inherits(s2.g, "oncosimulpop"))) {
            TTT <- FALSE
        } else {
            TTT <- TRUE
        }
        if(TTT) {
        summary(s2.g)[, c(1, 2, 3, 8, 9)]
        summary(s2.ng)[, c(1, 2, 3, 8, 9)]
        TTT <- c(TTT,  wilcox.test(mutsPerClone(s2.g),
                                   mutsPerClone(s2.ng),
                                   alternative = "greater")$p.value < p.value.threshold)
        TTT <- c(TTT,  wilcox.test(summary(s2.g)$NumClones,
                                   summary(s2.ng)$NumClones,
                                   alternative = "greater")$p.value < p.value.threshold)
        gc()
        }

        
        if( all(TTT) ) break;
    }
    cat(paste("\n done tries", tries, "\n"))
    expect_true(all(TTT))
})
cat("\n", date(), "\n")

cat("\n", date(), "\n")
test_that("McFL: Without initmutant", {
    
    max.tries <- 4
    for(tries in 1:max.tries) {


        TTT <- NULL
         cat("\n mcfls2: a runif is", runif(1), "\n")
    s2 <- 2.0
    ft <- 250 ## which is similar to size, as we are in the plateau
    pops <- 100
    lni <- 1 ## no fitness effects genes
    fni <- 50 ## fitness effects genes
    no <- 1e3
    mu <- 1e-5 ## easier to see
    ## noInt have no fitness effects, but can accumulate mutations
    ni <- rep(0, lni)
    ## Those with fitness effects in one module, so
    ## neither fitness nor mut. rate blow up
    gn <- paste(paste0("a", 1:fni), collapse = ", ")
    f2 <- allFitnessEffects(epistasis = c("A" = s2),
                            geneToModule = c("A" = gn),
                            noIntGenes = ni)
    cat("\n mcfls2a: a runif is", runif(1), "\n")
    s2.ng <- oncoSimulPop(pops,
                          f2,
                          mu = mu,
                          mutationPropGrowth = FALSE,
                          finalTime =ft,
                          initSize = no, keepEvery = 5,
                          onlyCancer = FALSE, detectionProb = NA, model = "McFL",
                          seed = NULL, mc.cores = 2)
    gc(); 
    cat("\n mcfls2b: a runif is", runif(1), "\n")
    s2.g <- oncoSimulPop(pops,
                         f2,
                         mu = mu,
                         mutationPropGrowth = TRUE,
                         finalTime =ft,
                         initSize = no, keepEvery = 5, 
                         onlyCancer = FALSE, detectionProb = NA, model = "McFL",
                         seed = NULL, mc.cores = 2)
    summary(s2.g)[, c(1, 2, 3, 8, 9)]
    summary(s2.ng)[, c(1, 2, 3, 8, 9)]
    TTT <- c(TTT,  wilcox.test(mutsPerClone(s2.g),
                             mutsPerClone(s2.ng),
                             alternative = "greater")$p.value < p.value.threshold)
    TTT <- c(TTT,  wilcox.test(summary(s2.g)$NumClones,
                             summary(s2.ng)$NumClones,
                             alternative = "greater")$p.value < p.value.threshold)


        
    ## summary(mutsPerClone(s2.g))
    ## summary(mutsPerClone(s2.ng))
    ## summary(summary(s2.g)$NumClones)
    ## summary(summary(s2.ng)$NumClones)
gc() 
if( all(TTT) ) break;
    }
    cat(paste("\n done tries", tries, "\n"))
    expect_true(all(TTT))
})
cat("\n", date(), "\n")





### As before, but fixing final population size.
cat("\n", date(), "\n")
test_that("detectionSize. Without initmutant", {
    
      max.tries <- 4
    for(tries in 1:max.tries) {

        TTT <- NULL
       cat("\n s3FPS: a runif is", runif(1), "\n")
    pops <- 200
    lni <- 1 ## no fitness effects genes
    fni <- 50 ## fitness effects genes
    no <- 1e3
    ft <- 10 ## 5
    s3 <- 3.0
    mu <- 5e-5 ## easier to see
    ## noInt have no fitness effects, but can accumulate mutations
    ni <- rep(0, lni)
    ## Those with fitness effects in one module, so
    ## neither fitness nor mut. rate blow up
    gn <- paste(paste0("a", 1:fni), collapse = ", ")
    f3 <- allFitnessEffects(epistasis = c("A" = s3),
                            geneToModule = c("A" = gn),
                            noIntGenes = ni)
    
    
    cat("\n s3FPSa: a runif is", runif(1), "\n")
    s3.ng <- oncoSimulPop(pops,
                          f3,
                          mu = mu,
                          mutationPropGrowth = FALSE,
                          detectionDrivers = 9999,
                          detectionSize = 5e5,
                          sampleEvery = 0.01,
                          keepEvery = 1,
                          finalTime =ft,
                          initSize = no,
                          onlyCancer = FALSE, detectionProb = NA,
                          seed = NULL, mc.cores = 2)
    
    
    cat("\n s3FPSb: a runif is", runif(1), "\n")
    s3.g <- oncoSimulPop(pops,
                         f3,
                         mu = mu,
                         mutationPropGrowth = TRUE,
                         detectionDrivers = 9999,
                          detectionSize = 5e5,
                          sampleEvery = 0.01,
                          keepEvery = 1,
                         finalTime =ft,
                         initSize = no,
                         onlyCancer = FALSE, detectionProb = NA,
                         seed = NULL, mc.cores = 2)
    ## summary(s3.g)[, c(1, 2, 3, 8, 9)]
    ## summary(s3.ng)[, c(1, 2, 3, 8, 9)]
    TTT <- c(TTT,  wilcox.test(mutsPerClone(s3.g),
                             mutsPerClone(s3.ng),
                             alternative = "greater")$p.value < p.value.threshold)
    TTT <- c(TTT,  wilcox.test(summary(s3.g)$NumClones,
                             summary(s3.ng)$NumClones,
                             alternative = "greater")$p.value < p.value.threshold)
    summary(summary(s3.g)[, 2])
    summary(summary(s3.ng)[, 2])
    
gc() 
if( all(TTT) ) break;
    }
    cat(paste("\n done tries", tries, "\n"))
    expect_true(all(TTT))

})
cat("\n", date(), "\n")

cat("\n", date(), "\n")
test_that("detectionSize. Without initmutant, 2", {
    
    ## More of the above. Use smaller s2 and smaller mutation, but then to
    ## see it reliably you need large final popsize
    ## Can fail sometimes because differences are small
    
     max.tries <- 4
     for(tries in 1:max.tries) {

         
         TTT <- NULL
         cat("\n s2FPS: a runif is", runif(1), "\n")
    s2 <- 1.0
    ft <- 500  ## very large, so we stop on size
    pops <- 200 
    fni <- 300 ## fitness effects genes
    no <- 1e3
    mu <- 5e-6 
    ## noInt have no fitness effects, but can accumulate mutations
    ## Those with fitness effects in one module, so
    ## neither fitness nor mut. rate blow up
    gn <- paste(paste0("a", 1:fni), collapse = ", ")
    f2 <- allFitnessEffects(epistasis = c("A" = s2),
                            geneToModule = c("A" = gn))
    
    
    cat("\n s2FPSa: a runif is", runif(1), "\n")
    s2.ng <- oncoSimulPop(pops,
                          f2,
                          mu = mu,
                          detectionDrivers = 9999,
                          detectionSize = 5e6,
                          sampleEvery = 0.01,
                          mutationPropGrowth = FALSE,
                          finalTime =ft,
                          initSize = no, keepEvery = 1,
                          onlyCancer = FALSE, detectionProb = NA,
                          seed = NULL, mc.cores = 2)
    gc(); 
    
    cat("\n s2FPSb: a runif is", runif(1), "\n")
    s2.g <- oncoSimulPop(pops,
                         f2,
                         mu = mu,
                         detectionDrivers = 9999,
                          detectionSize = 5e6,
                          sampleEvery = 0.01,
                         mutationPropGrowth = TRUE,
                         finalTime =ft,
                         initSize = no, keepEvery = 1,
                         onlyCancer = FALSE, detectionProb = NA,
                         seed = NULL, mc.cores = 2)
    ## summary(s2.g)[, c(1, 2, 3, 8, 9)]
    ## summary(s2.ng)[, c(1, 2, 3, 8, 9)]
    TTT <- c(TTT,  wilcox.test(mutsPerClone(s2.g),
                             mutsPerClone(s2.ng),
                             alternative = "greater")$p.value < p.value.threshold)
    TTT <- c(TTT,  wilcox.test(summary(s2.g)$NumClones,
                             summary(s2.ng)$NumClones,
                             alternative = "greater")$p.value < p.value.threshold)
    summary(summary(s2.g)[, 2])
    summary(summary(s2.ng)[, 2])

    gc() 
if( all(TTT) ) break;
    }
    cat(paste("\n done tries", tries, "\n"))
    expect_true(all(TTT))

})
cat("\n", date(), "\n")

cat("\n", date(), "\n")
test_that("detectionSize. McFL: Without initmutant", {
    ## with McFL limiting popSize in the simuls is not that relevant, as
    ## already limited.
    
     max.tries <- 4
    for(tries in 1:max.tries) {
TTT <- NULL
    
    cat("\n FPSmcfls2: a runif is", runif(1), "\n")
    s2 <- 2.0
    ft <- 250
    pops <- 300 ## 200
    lni <- 1 ## no fitness effects genes
    fni <- 50 ## fitness effects genes
    no <- 1e3
    mu <- 1e-5 ## easier to see
    ## noInt have no fitness effects, but can accumulate mutations
    ni <- rep(0, lni)
    ## Those with fitness effects in one module, so
    ## neither fitness nor mut. rate blow up
    gn <- paste(paste0("a", 1:fni), collapse = ", ")
    f2 <- allFitnessEffects(epistasis = c("A" = s2),
                            geneToModule = c("A" = gn),
                            noIntGenes = ni)
    
    
    cat("\n FPSmcfls2a: a runif is", runif(1), "\n")
    s2.ng <- oncoSimulPop(pops,
                          f2,
                          mu = mu,
                          mutationPropGrowth = FALSE,
                          finalTime =ft,
                          detectionDrivers = 9999,
                          detectionSize = 1e4,
                          sampleEvery = 0.01,
                          initSize = no, keepEvery = 5,
                          onlyCancer = FALSE, detectionProb = NA, model = "McFL",
                          seed = NULL, mc.cores = 2)
    gc(); 
    
    cat("\n FPSmcfls2b: a runif is", runif(1), "\n")
    s2.g <- oncoSimulPop(pops,
                         f2,
                         mu = mu,
                         mutationPropGrowth = TRUE,
                         finalTime =ft,
                         detectionDrivers = 9999,
                          detectionSize = 1e4,
                          sampleEvery = 0.01,
                         initSize = no, keepEvery = 5, 
                         onlyCancer = FALSE, detectionProb = NA, model = "McFL",
                         seed = NULL, mc.cores = 2)
    ## summary(s2.g)[, c(1, 2, 3, 8, 9)]
    ## summary(s2.ng)[, c(1, 2, 3, 8, 9)]
    TTT <- c(TTT,  wilcox.test(mutsPerClone(s2.g),
                             mutsPerClone(s2.ng),
                             alternative = "greater")$p.value < p.value.threshold)
    TTT <- c(TTT,  wilcox.test(summary(s2.g)$NumClones,
                             summary(s2.ng)$NumClones,
                             alternative = "greater")$p.value < p.value.threshold)
    summary(summary(s2.g)[, 2])
    summary(summary(s2.ng)[, 2])
    ## summary(mutsPerClone(s2.g))
    ## summary(mutsPerClone(s2.ng))
    ## summary(summary(s2.g)$NumClones)
    ## summary(summary(s2.ng)$NumClones)
gc() 
if( all(TTT) ) break;
    }
    cat(paste("\n done tries", tries, "\n"))
    expect_true(all(TTT))

})
cat("\n", date(), "\n")










cat(paste("\n Ending mutPropGrwoth-long at", date(), "\n"))






## This comparison is actually a mess.  If we stop on size, the very
## fast growing ones have been growing for a lot less, so they cannot
## have accumulated the same number of mutations.

## We are mixing phenomena here. This is a bad test.

## date()
## test_that("Ordering of number of clones with mutpropgrowth", {

    

##     ## So here, we stop on time, not size, because of the huge differences
##     ## in growth rate.
    
##      max.tries <- 4
##      for(tries in 1:max.tries) {

         
##          TTT <- NULL
##     cat("\n omp1: a runif is", runif(1), "\n")
##     pops <- 100
##     lni <- 200
##     no <- 5e3
##          ni <- c(4, 2, rep(0, lni))
##          ft <- 2
##          ds <- 5e8
##     names(ni) <- c("a", "b", paste0("n", seq.int(lni)))
##     fe <- allFitnessEffects(noIntGenes = ni)
##     cat("\n omp1a: a runif is", runif(1), "\n")
##     nca <- oncoSimulPop(pops, fe, finalTime = ft, detectionSize = ds,
##                         mutationPropGrowth = TRUE,
##                         initSize = no, sampleEvery = 0.01,
##                         initMutant = "a",
##                         onlyCancer = FALSE, detectionProb = NA, seed = NULL, mc.cores = 2)
##     cat("\n omp1b: a runif is", runif(1), "\n")
##     ncb <- oncoSimulPop(pops, fe, finalTime = ft, detectionSize = ds,
##                         mutationPropGrowth = TRUE,
##                         initSize = no,sampleEvery = 0.01,
##                         initMutant = "b",
##                         onlyCancer = FALSE, detectionProb = NA, seed = NULL, mc.cores = 2)
##     cat("\n omp1c: a runif is", runif(1), "\n")
##     nca2 <- oncoSimulPop(pops, fe, finalTime = ft, detectionSize = ds,
##                         mutationPropGrowth = FALSE,
##                         initSize = no,sampleEvery = 0.01,
##                         initMutant = "a",
##                         onlyCancer = FALSE, detectionProb = NA, seed = NULL, mc.cores = 2)
##     cat("\n omp1d: a runif is", runif(1), "\n")
##     ncb2 <- oncoSimulPop(pops, fe, finalTime = ft, detectionSize = ds,
##                         mutationPropGrowth = FALSE,
##                         initSize = no,sampleEvery = 0.01,
##                         initMutant = "b",
##                         onlyCancer = FALSE, detectionProb = NA, seed = NULL, mc.cores = 2)
##          summary(nca)[, c(1, 2, 3, 8, 9)]
##          summary(ncb)[, c(1, 2, 3, 8, 9)]
##          summary(nca2)[, c(1, 2, 3, 8, 9)]
##          summary(ncb2)[, c(1, 2, 3, 8, 9)]
##          ## The real comparison
##     TTT <- c(TTT,  wilcox.test(summary(nca)$NumClones,
##                              summary(nca2)$NumClones,
##                              alternative = "greater")$p.value < p.value.threshold)
##     TTT <- c(TTT,  wilcox.test(summary(nca)$NumClones,
##                              summary(ncb)$NumClones,
##                              alternative = "greater")$p.value < p.value.threshold)
##     TTT <- c(TTT,  wilcox.test(summary(ncb)$NumClones,
##                              summary(ncb2)$NumClones,
##                              alternative = "greater")$p.value < p.value.threshold)

         
## if( all(TTT) ) break;
##     }
##     cat(paste("\n done tries", tries, "\n"))
##     expect_true(all(TTT))
## })
## cat("\n", date(), "\n")



## cat("\n", date(), "\n")
## test_that("Ordering of number of clones with mutpropgrowth, McFL", {
    
##      max.tries <- 4
##      for(tries in 1:max.tries) {
         
##          TTT <- NULL
##          cat("\n omp2: a runif is", runif(1), "\n")
##          pops <- 10
##          lni <- 200
##          no <- 5e3
##          ni <- c(5, 2, rep(0, lni))
##          ft <- 20
##          ds <- 1e5
##          names(ni) <- c("a", "b", paste0("n", seq.int(lni)))
##          fe <- allFitnessEffects(noIntGenes = ni)
##          cat("\n omp2a: a runif is", runif(1), "\n")
##          nca <- oncoSimulPop(pops, fe, finalTime = ft, detectionSize = ds,
##                              mutationPropGrowth = TRUE,
##                         initSize = no, model = "McFL",
##                         initMutant = "a",
##                         onlyCancer = FALSE, detectionProb = NA, seed = NULL, mc.cores = 2)
##     cat("\n omp2b: a runif is", runif(1), "\n")
##     ncb <- oncoSimulPop(pops, fe, finalTime = ft, detectionSize = ds,
##                         mutationPropGrowth = TRUE,
##                         initSize = no, model = "McFL",
##                         initMutant = "b",
##                         onlyCancer = FALSE, detectionProb = NA, seed = NULL, mc.cores = 2)
##     cat("\n omp2c: a runif is", runif(1), "\n")
##     nca2 <- oncoSimulPop(pops, fe, finalTime = ft, detectionSize = ds,
##                         mutationPropGrowth = FALSE,
##                         initSize = no, model = "McFL",
##                         initMutant = "a",
##                         onlyCancer = FALSE, detectionProb = NA, seed = NULL, mc.cores = 2)
##     cat("\n omp2d: a runif is", runif(1), "\n")
##     ncb2 <- oncoSimulPop(pops, fe, finalTime = ft, detectionSize = ds,
##                         mutationPropGrowth = FALSE,
##                         initSize = no, model = "McFL",
##                         initMutant = "b",
##                         onlyCancer = FALSE, detectionProb = NA, seed = NULL, mc.cores = 2)
##          summary(nca)[, c(1, 2, 3, 8, 9)]
##          summary(ncb)[, c(1, 2, 3, 8, 9)]
##          summary(nca2)[, c(1, 2, 3, 8, 9)]
##          summary(ncb2)[, c(1, 2, 3, 8, 9)]
##     ## The real comparison
##     TTT <- c(TTT,  wilcox.test(summary(nca)$NumClones,
##                              summary(nca2)$NumClones,
##                              alternative = "greater")$p.value < p.value.threshold)
##     TTT <- c(TTT,  wilcox.test(summary(nca)$NumClones,
##                              summary(ncb)$NumClones,
##                              alternative = "greater")$p.value < p.value.threshold)
##     TTT <- c(TTT,  wilcox.test(summary(ncb)$NumClones,
##                              summary(ncb2)$NumClones,
##                              alternative = "greater")$p.value < p.value.threshold)



##          if( all(TTT) ) break;
##     }
##     cat(paste("\n done tries", tries, "\n"))
##     expect_true(all(TTT))
## })







## ## This generally works, but not always. Because: a) starting with a you
## ## can get a mutation in b. Or starting with b you can get a mutation in
## ## a. When either happens is stochastic. And we are also mixing the
## ## mutation proportional to birth with the simple increase in clones due
## ## to larger pop sizes. So this ain't a good test. When you start, say,
## ## from b and you hit a soon, you can get HUGE growth rates and that
## ## induces huge variability, slow times, etc.

## cat("\n", date(), "\n")
## test_that("Ordering of number of clones and mutsPerClone with mutpropgrowth, 1", {
##     
##     
##     cat("\n mpc1: a runif is", runif(1), "\n")
##     ft <- 2.5
##     pops <- 200
##     lni <- 500 ## with, say, 40 or a 100, sometimes fails the comparisons
##                ## with small differences.
##     no <- 10
##     ni <- c(5, 3, rep(0, lni))
##     names(ni) <- c("a", "b", paste0("n", seq.int(lni)))
##     fe <- allFitnessEffects(noIntGenes = ni)
##     
##     
##     cat("\n mpc1a: a runif is", runif(1), "\n")
##     nca <- oncoSimulPop(pops, fe, finalTime = ft,
##                         mutationPropGrowth = TRUE,
##                         initSize = no,
##                         initMutant = "a",
##                         onlyCancer = FALSE, detectionProb = NA, seed = NULL, mc.cores = 2)
##     gc()
##     
##     
##     cat("\n mpc1b: a runif is", runif(1), "\n")
##     ncb <- oncoSimulPop(pops, fe, finalTime = ft,
##                         mutationPropGrowth = TRUE,
##                         initSize = no,
##                         initMutant = "b",
##                         onlyCancer = FALSE, detectionProb = NA, seed = NULL, mc.cores = 2)
##     gc()
##     
##     
##     cat("\n mpc1c: a runif is", runif(1), "\n")
##     nca2 <- oncoSimulPop(pops, fe, finalTime = ft,
##                          mutationPropGrowth = FALSE,
##                          initSize = no,
##                          initMutant = "a",
##                          onlyCancer = FALSE, detectionProb = NA, seed = NULL, mc.cores = 2)
##     
##     
##     cat("\n mpc1d: a runif is", runif(1), "\n")
##     ncb2 <- oncoSimulPop(pops, fe, finalTime = ft,
##                          mutationPropGrowth = FALSE,
##                          initSize = no,
##                          initMutant = "b",
##                          onlyCancer = FALSE, detectionProb = NA, seed = NULL, mc.cores = 2)
##     gc()
##     expect_true(var(summary(nca)$NumClones) > 1e-4)
##     expect_true(var(summary(ncb)$NumClones) > 1e-4)
##     expect_true(var(summary(nca2)$NumClones) > 1e-4)
##     expect_true(var(summary(ncb2)$NumClones) > 1e-4)
##     ## The real comparison
##     expect_true( median(summary(nca)$NumClones) >
##                  median(summary(ncb)$NumClones))
##     expect_true( median(summary(ncb)$NumClones) >
##                  median(summary(ncb2)$NumClones))
##     expect_true( mean(mutsPerClone(nca)) >
##                  mean(mutsPerClone(ncb)))
##     expect_true( mean(mutsPerClone(ncb)) >
##                  mean(mutsPerClone(ncb2)))
##     ## These can fail in this case, since small diffs. as small mutlipliers
##     expect_true( mean(mutsPerClone(nca)) >
##                  mean(mutsPerClone(nca2)))
##     expect_true( median(summary(nca)$NumClones) >
##                  median(summary(nca2)$NumClones))
##     ## In this cases, we would expect differences in total population size
##     ## between a and a2, but minor or non detectable between b and b2. In
##     ## a we expect them because, since those affect a lot the mutation
##     ## rate, we expect them to get b faster, and thus grow faster noticeably.
## gc() 
## })
## cat("\n", date(), "\n")




## cat("\n", date(), "\n")
## test_that("Ordering of number of clones and mutsPerClone with mutpropgrowth, 3", {
##     
##     
##     cat("\n mpc3: a runif is", runif(1), "\n")
##     ## The s coefficient is small, and so small differences. Here, much large
##     ## mu
##     ft <- 12 ## going beyond 13 or so, gets it to bail because of reaching max
##     ## pop in nca
##     pops <- 200
##     lni <- 50
##     no <- 10
##     ni <- c(1.0, 0.8, rep(0, lni))
##     mu <- 1e-5
##     names(ni) <- c("a", "b", paste0("n", seq.int(lni)))
##     fe <- allFitnessEffects(noIntGenes = ni)
##     
##     
##     cat("\n mpc3a: a runif is", runif(1), "\n")
##     nca <- oncoSimulPop(pops, fe, finalTime = ft,
##                         mutationPropGrowth = TRUE,
##                         mu = mu,
##                         initSize = no, keepEvery = 1,
##                         initMutant = "a",
##                         onlyCancer = FALSE, detectionProb = NA, seed = NULL, mc.cores = 2)
##     gc()
##     
##     
##     cat("\n mpc3b: a runif is", runif(1), "\n")
##     ncb <- oncoSimulPop(pops, fe, finalTime = ft,
##                         mutationPropGrowth = TRUE,
##                         mu = mu,
##                         initSize = no, keepEvery = 1,
##                         initMutant = "b",
##                         onlyCancer = FALSE, detectionProb = NA, seed = NULL, mc.cores = 2)
##     gc()
##     
##     
##     cat("\n mpc3c: a runif is", runif(1), "\n")
##     nca2 <- oncoSimulPop(pops, fe, finalTime = ft,
##                          mutationPropGrowth = FALSE,
##                          mu = mu,
##                          initSize = no, keepEvery = 1,
##                          initMutant = "a",
##                          onlyCancer = FALSE, detectionProb = NA, seed = NULL, mc.cores = 2)
##     gc()
##     
##     
##     cat("\n mpc3d: a runif is", runif(1), "\n")
##     ncb2 <- oncoSimulPop(pops, fe, finalTime = ft,
##                          mutationPropGrowth = FALSE,
##                          mu = mu,                     
##                          initSize = no, keepEvery = 1,
##                          initMutant = "b",
##                          onlyCancer = FALSE, detectionProb = NA, seed = NULL, mc.cores = 2)
##     gc()
##     expect_true(var(summary(nca)$NumClones) > 1e-4)
##     expect_true(var(summary(ncb)$NumClones) > 1e-4)
##     expect_true(var(summary(nca2)$NumClones) > 1e-4)
##     expect_true(var(summary(ncb2)$NumClones) > 1e-4)
##     ## The real comparison
##     expect_true( median(summary(nca)$NumClones) >
##                  median(summary(ncb)$NumClones))
##     expect_true( median(summary(ncb)$NumClones) >
##                  median(summary(ncb2)$NumClones))
##     expect_true( mean(mutsPerClone(nca)) >
##                  mean(mutsPerClone(ncb)))
##     expect_true( mean(mutsPerClone(ncb)) >
##                  mean(mutsPerClone(ncb2)))
##     ## These can fail in this case, since small diffs. as small mutlipliers
##     expect_true( mean(mutsPerClone(nca)) >
##                  mean(mutsPerClone(nca2)))
##     expect_true( median(summary(nca)$NumClones) >
##                  median(summary(nca2)$NumClones))
## gc() 
## })


## cat("\n", date(), "\n")
## test_that("McFL: Ordering of number of clones and mutsPerClone with mutpropgrowth, 1", {
##     
##     
##     cat("\n mpcmcf1: a runif is", runif(1), "\n")
##     ft <- 20 ## unless large you rarely get triple, etc, mutatns
##     pops <- 200
##     lni <- 50 
##     no <- 1e3
##     ni <- c(3, 1.5, rep(0, lni))
##     names(ni) <- c("a", "b", paste0("n", seq.int(lni)))
##     fe <- allFitnessEffects(noIntGenes = ni)
##     
##     
##     cat("\n mpcm1a: a runif is", runif(1), "\n")
##     nca <- oncoSimulPop(pops, fe, finalTime = ft,
##                         mutationPropGrowth = TRUE,
##                         initSize = no, model = "McFL",
##                         initMutant = "a", keepEvery = 1,
##                         onlyCancer = FALSE, detectionProb = NA, seed = NULL, mc.cores = 2)
##     gc()
##     
##     
##     cat("\n mpcm1b: a runif is", runif(1), "\n")
##     ncb <- oncoSimulPop(pops, fe, finalTime = ft,
##                         mutationPropGrowth = TRUE,
##                         initSize = no, model = "McFL",
##                         initMutant = "b", keepEvery = 1,
##                         onlyCancer = FALSE, detectionProb = NA, seed = NULL, mc.cores = 2)
##     gc()
##     
##     
##     cat("\n mpcm1c: a runif is", runif(1), "\n")
##     nca2 <- oncoSimulPop(pops, fe, finalTime = ft,
##                          mutationPropGrowth = FALSE,
##                          initSize = no, model = "McFL",
##                          initMutant = "a", keepEvery = 1,
##                          onlyCancer = FALSE, detectionProb = NA, seed = NULL, mc.cores = 2)
##     gc()
##     
##     
##     cat("\n mpcm1d: a runif is", runif(1), "\n")
##     ncb2 <- oncoSimulPop(pops, fe, finalTime = ft,
##                          mutationPropGrowth = FALSE,
##                          initSize = no, model = "McFL",
##                          initMutant = "b", keepEvery = 1,
##                          onlyCancer = FALSE, detectionProb = NA, seed = NULL, mc.cores = 2)
##     gc()
##     expect_true(var(summary(nca)$NumClones) > 1e-4)
##     expect_true(var(summary(ncb)$NumClones) > 1e-4)
##     expect_true(var(summary(nca2)$NumClones) > 1e-4)
##     expect_true(var(summary(ncb2)$NumClones) > 1e-4)
##     ## The real comparison
##     expect_true( median(summary(nca)$NumClones) >
##                  median(summary(ncb)$NumClones))
##     expect_true( median(summary(ncb)$NumClones) >
##                  median(summary(ncb2)$NumClones))
##     expect_true( mean(mutsPerClone(nca)) >
##                  mean(mutsPerClone(ncb)))
##     expect_true( mean(mutsPerClone(ncb)) >
##                  mean(mutsPerClone(ncb2)))
##     ## These can fail in this case, since small diffs. as small mutlipliers
##     expect_true( mean(mutsPerClone(nca)) >
##                  mean(mutsPerClone(nca2)))
##     expect_true( median(summary(nca)$NumClones) >
##                  median(summary(nca2)$NumClones))
## gc() 
## })
## cat("\n", date(), "\n")

## ## A variation of the former
## cat("\n", date(), "\n")
## test_that("McFL: Ordering of number of clones and mutsPerClone with mutpropgrowth, 2", {
##     
##     
##     cat("\n mpcmcf2: a runif is", runif(1), "\n")
##     ## Increase ft
##     ft <- 350 
##     pops <- 150
##     lni <- 30
##     no <- 1e3
##     ni <- c(2, 0.8, rep(0, lni))
##     names(ni) <- c("a", "b", paste0("n", seq.int(lni)))
##     fe <- allFitnessEffects(noIntGenes = ni)
##     
##     
##     cat("\n mpcmcf2a: a runif is", runif(1), "\n")
##     nca <- oncoSimulPop(pops, fe, finalTime = ft,
##                         mutationPropGrowth = TRUE,
##                         initSize = no, model = "McFL",
##                         initMutant = "a", keepEvery = 1,
##                         onlyCancer = FALSE, detectionProb = NA, seed = NULL, mc.cores = 2)
##     gc()
##     
##     
##     cat("\n mpcmcf2b: a runif is", runif(1), "\n")
##     ncb <- oncoSimulPop(pops, fe, finalTime = ft,
##                         mutationPropGrowth = TRUE,
##                         initSize = no, model = "McFL",
##                         initMutant = "b", keepEvery = 1,
##                         onlyCancer = FALSE, detectionProb = NA, seed = NULL, mc.cores = 2)
##     gc()
##     
##     
##     cat("\n mpcmcf2c: a runif is", runif(1), "\n")
##     nca2 <- oncoSimulPop(pops, fe, finalTime = ft,
##                          mutationPropGrowth = FALSE,
##                          initSize = no, model = "McFL",
##                          initMutant = "a", keepEvery = 1,
##                          onlyCancer = FALSE, detectionProb = NA, seed = NULL, mc.cores = 2)
##     gc()
##     
##     
##     cat("\n mpcmcf2d: a runif is", runif(1), "\n")
##     ncb2 <- oncoSimulPop(pops, fe, finalTime = ft,
##                          mutationPropGrowth = FALSE,
##                          initSize = no, model = "McFL",
##                          initMutant = "b", keepEvery = 1,
##                          onlyCancer = FALSE, detectionProb = NA, seed = NULL, mc.cores = 2)
##     gc()
##     ## summary(nca)[, c(1, 2, 3, 8, 9)]
##     ## summary(nca2)[, c(1, 2, 3, 8, 9)]
##     ## summary(ncb)[, c(1, 2, 3, 8, 9)]
##     ## summary(ncb2)[, c(1, 2, 3, 8, 9)]
##     expect_true(var(summary(nca)$NumClones) > 1e-4)
##     expect_true(var(summary(ncb)$NumClones) > 1e-4)
##     expect_true(var(summary(nca2)$NumClones) > 1e-4)
##     expect_true(var(summary(ncb2)$NumClones) > 1e-4)
##     ## The real comparison
##     expect_true( median(summary(nca)$NumClones) >
##                  median(summary(ncb)$NumClones))
##     expect_true( median(summary(ncb)$NumClones) >
##                  median(summary(ncb2)$NumClones))
##     expect_true( median(summary(nca)$NumClones) >
##                  median(summary(nca2)$NumClones))
##     expect_true( mean(mutsPerClone(nca)) >
##                  mean(mutsPerClone(ncb)))
##     expect_true( mean(mutsPerClone(nca)) >
##                  mean(mutsPerClone(nca2)))
##     expect_true( mean(mutsPerClone(ncb)) >
##                  mean(mutsPerClone(ncb2)))
##     ## These can fail in this case, since small diffs. as small mutlipliers
## gc() 
## })
## date()


## ## When diffs are very tiny in s, as we increase ft, it is easy to get a
## ## mutant in the second gene with non-zero s and thus the growth rates
## ## will be the same and there will be no differences. This is exacerbated
## ## with McFL as the less fit clone disappears

## ## We take this further below. When we init from a1 below, we start from a
## ## clone with smaller growth rate. But there are many b to jump
## ## to. However, if you start from b, you can only increase birth rate
## ## mutating exactly one a. Thus, over moderate finalTimes, starting a a
## ## will lead to more clones, etc

## ## But, as above, this can fail if you hit the unlikely gene. 

## cat("\n", date(), "\n")
## test_that("McFL: Ordering of number of clones and mutsPerClone with mutpropgrowth, and different mmodules",{
##     
##     
##     cat("\n mpcmcf3: a runif is", runif(1), "\n")
##     ft <- 10 
##     pops <- 200
##     mu <- 1e-5
##     lni <- 10
##     fni <- 50
##     no <- 1e4
##     s1 <- 0.9
##     s2 <- 1
##     gn <- paste(paste0("b", 1:fni), collapse = ", ")
##     ni <- rep(0, lni)
##     names(ni) <- paste0("n", seq.int(lni))
##     f1 <- allFitnessEffects(epistasis = c("A" = s1,
##                                           "B" = s2),
##                             geneToModule = c("A" = "a1",
##                                              "B" = gn),
##                             noIntGenes = ni)
##     
##     
##     cat("\n mpcmcf3a: a runif is", runif(1), "\n")
##     nca <- oncoSimulPop(pops, f1, finalTime = ft,
##                         mutationPropGrowth = TRUE,
##                         mu = mu, keepEvery = 1,
##                         initSize = no, model = "McFL",
##                         initMutant = "a1", 
##                         onlyCancer = FALSE, detectionProb = NA, seed = NULL, mc.cores = 2)
##     
##     gc()
##     
##     cat("\n mpcmcf3b: a runif is", runif(1), "\n")
##     ncb <- oncoSimulPop(pops, f1, finalTime = ft,
##                         mutationPropGrowth = TRUE,
##                         mu = mu, keepEvery = 1,
##                         initSize = no, model = "McFL",
##                         initMutant = "b1",
##                         onlyCancer = FALSE, detectionProb = NA, seed = NULL, mc.cores = 2)
##     ## OK, but it also happens below just because init a1 eventually grows
##     ## faster, so larger pop, so more mutants, etc
##     ## We could try just counting the number of mutation events, but we would
##     ## need to standardize by total time of life over all cells.
##     ## As growth rates are faster (those that start with a, here, and then
##     ## get a b), even when mutationPropGrowth is false, the population
##     ## gets larger. And thus, it is easier to get a clone with three
##     ## mutations, etc.
##     gc()
##     
##     
##     cat("\n mpcmcf3c: a runif is", runif(1), "\n")
##     nca2 <- oncoSimulPop(pops, f1, finalTime = ft,
##                          mutationPropGrowth = FALSE,
##                          mu = mu, keepEvery = 1,
##                          initSize = no, model = "McFL",
##                          initMutant = "a1",
##                          onlyCancer = FALSE, detectionProb = NA, seed = NULL, mc.cores = 2)
##     gc(); 
##     
##     cat("\n mpcmcf3d: a runif is", runif(1), "\n")
##     ncb2 <- oncoSimulPop(pops, f1, finalTime = ft,
##                          mutationPropGrowth = FALSE,
##                          mu = mu, keepEvery = 1,                    
##                          initSize = no, model = "McFL",
##                          initMutant = "b1",
##                          onlyCancer = FALSE, detectionProb = NA, seed = NULL, mc.cores = 2)
##     expect_true(var(summary(nca)$NumClones) > 1e-4)
##     expect_true(var(summary(ncb)$NumClones) > 1e-4)
##     expect_true(var(summary(nca2)$NumClones) > 1e-4)
##     expect_true(var(summary(ncb2)$NumClones) > 1e-4)
##     ## The real comparison
##     expect_true( median(summary(nca)$NumClones) >
##                  median(summary(ncb)$NumClones))
##     expect_true( median(summary(ncb)$NumClones) >
##                  median(summary(ncb2)$NumClones))
##     expect_true( mean(mutsPerClone(nca)) >
##                  mean(mutsPerClone(ncb)))
##     expect_true( mean(mutsPerClone(ncb)) >
##                  mean(mutsPerClone(ncb2)))
##     ## These can fail in this case, since small diffs. as small mutlipliers
##     expect_true( mean(mutsPerClone(nca)) >
##                  mean(mutsPerClone(nca2)))
##     expect_true( median(summary(nca)$NumClones) >
##                  median(summary(nca2)$NumClones))
## gc() 
## })





## cat("\n", date(), "\n")
## test_that("Ordering of number of clones and mutsPerClone with initMutant and modules, oncoSimulSample", {
##     
##     
##     cat("\n ossmpc1: a runif is", runif(1), "\n")
##     ft <- 2  
##     pops <- 400
##     lni <- 500 ## with, say, 40 or a 100, sometimes fails the comparisons
##                ## with small differences.
##     no <- 10
##     x <- 1e-40
##     ni <- c(5, 3, rep(0, lni))
##     names(ni) <- c("a", "b", paste0("n", seq.int(lni)))
##     fe <- allFitnessEffects(noIntGenes = ni)
##     
##     
##     cat("\n ossmpc1a: a runif is", runif(1), "\n")
##     nca <- oncoSimulSample(pops, fe, finalTime = ft,
##                         mutationPropGrowth = TRUE,
##                         initSize = no, 
##                         initMutant = "a",
##                         onlyCancer = FALSE, detectionProb = NA,  sampleEvery = 0.01,
##                           detectionSize = 1e9,
##                           detectionDrivers = 99,
##                           seed =NULL, max.wall.time = 3000,
##                           thresholdWhole = x)
##     gc(); 
##     
##     cat("\n ossmpc1b: a runif is", runif(1), "\n")
##     ncb <- oncoSimulSample(pops, fe, finalTime = ft,
##                         mutationPropGrowth = TRUE,
##                         initSize = no, 
##                         initMutant = "b",
##                         onlyCancer = FALSE, detectionProb = NA, sampleEvery = 0.01,
##                           detectionSize = 1e9,
##                           detectionDrivers = 99,
##                           seed =NULL, max.wall.time = 3000,
##                           thresholdWhole = x)
##     gc(); 
##     
##     cat("\n ossmpc1c: a runif is", runif(1), "\n")
##     nca2 <- oncoSimulSample(pops, fe, finalTime = ft,
##                          mutationPropGrowth = FALSE,
##                          initSize = no, 
##                          initMutant = "a",
##                          onlyCancer = FALSE, detectionProb = NA, sampleEvery = 0.01,
##                           detectionSize = 1e9,
##                           detectionDrivers = 99,
##                           seed =NULL, max.wall.time = 3000,
##                           thresholdWhole = x)
##     gc(); 
##     
##     cat("\n ossmpc1d: a runif is", runif(1), "\n")
##     ncb2 <- oncoSimulSample(pops, fe, finalTime = ft,
##                          mutationPropGrowth = FALSE,
##                          initSize = no, 
##                          initMutant = "b",
##                          onlyCancer = FALSE, detectionProb = NA, sampleEvery = 0.01,
##                           detectionSize = 1e9,
##                           detectionDrivers = 99,
##                           seed =NULL, max.wall.time = 3000,
##                          thresholdWhole = x)
##     ## nca$popSummary[1:5, c(1:3, 8:9)]
##     ## ncb$popSummary[1:5, c(1:3, 8:9)]
##     ## nca2$popSummary[1:5, c(1:3, 8:9)]
##     ## ncb2$popSummary[1:5, c(1:3, 8:9)]
##     ## summary(nca$popSummary[, "NumClones"])
##     ## summary(ncb$popSummary[, "NumClones"])
##     ## summary(nca2$popSummary[, "NumClones"])
##     ## summary(ncb2$popSummary[, "NumClones"])
##     ## cat("\n mutsperclone\n")
##     ## summary(mutsPerCloneNCA <- rowSums(nca$popSample))
##     ## summary(mutsPerCloneNCB <- rowSums(ncb$popSample))
##     ## summary(mutsPerCloneNCA2 <- rowSums(nca2$popSample))
##     ## summary(mutsPerCloneNCB2 <- rowSums(ncb2$popSample))
##     ## (cc1 <- colSums(nca$popSample))
##     ## Because of how we do the "mutsPerClone", mutsPerClone almost the
##     ## same ans number of clones (except for the few cases when a clone
##     ## has gone extinct)
##     if(!is.array(nca$popSample)) {
##         ## occasionally, I get funny things
##         warning("nca$popSample not an array")
##         cat(class(nca$popSample))
##         cat(dim(nca$popSample))
##         cat(str(nca$popSample))
##     }
##     mutsPerCloneNCA <- rowSums(nca$popSample)
##     mutsPerCloneNCB <- rowSums(ncb$popSample)
##     mutsPerCloneNCA2 <- rowSums(nca2$popSample)
##     mutsPerCloneNCB2 <- rowSums(ncb2$popSample)
##     expect_true( median(mutsPerCloneNCA) >
##                  median(mutsPerCloneNCA2))
##     expect_true( median(mutsPerCloneNCB) >
##                  median(mutsPerCloneNCB2))
##     expect_true( median(mutsPerCloneNCA) >
##                  median(mutsPerCloneNCB))
##     expect_true( median(nca$popSummary[, "NumClones"]) >
##                  median(nca2$popSummary[, "NumClones"]))
##     expect_true( median(ncb$popSummary[, "NumClones"]) >
##                  median(ncb2$popSummary[, "NumClones"]))
##     expect_true( median(nca$popSummary[, "NumClones"]) >
##                  median(ncb$popSummary[, "NumClones"]))
## gc() 
## })
## cat("\n", date(), "\n")


## cat("\n", date(), "\n")
## test_that("Ordering of number of clones and mutsPerClone with initMutant and modules, oncoSimulSample, McFL", {
##     
##     
##     cat("\n ossmpc1McFL: a runif is", runif(1), "\n")
##     ft <- 2  
##     pops <- 200 ## 200
##     lni <- 500 ## with, say, 40 or a 100, sometimes fails the comparisons
##                ## with small differences.
##     no <- 100
##     x <- 1e-40
##     ni <- c(5, 3, rep(0, lni))
##     names(ni) <- c("a", "b", paste0("n", seq.int(lni)))
##     fe <- allFitnessEffects(noIntGenes = ni)
##     
##     
##     cat("\n ossmpc1McFLa: a runif is", runif(1), "\n")
##     nca <- oncoSimulSample(pops, fe, finalTime = ft,
##                         mutationPropGrowth = TRUE,
##                         initSize = no, 
##                         initMutant = "a", model = "McFL",
##                         onlyCancer = FALSE, detectionProb = NA,  sampleEvery = 0.01,
##                           detectionSize = 1e9,
##                           detectionDrivers = 99,
##                           seed =NULL, max.wall.time = 3000,
##                           thresholdWhole = x)
##     gc(); 
##     
##     cat("\n ossmpc1McFLb: a runif is", runif(1), "\n")
##     ncb <- oncoSimulSample(pops, fe, finalTime = ft,
##                         mutationPropGrowth = TRUE,
##                         initSize = no, 
##                         initMutant = "b", model = "McFL",
##                         onlyCancer = FALSE, detectionProb = NA, sampleEvery = 0.01,
##                           detectionSize = 1e9,
##                           detectionDrivers = 99,
##                           seed =NULL, max.wall.time = 3000,
##                           thresholdWhole = x)
##     gc(); 
##     
##     cat("\n ossmpc1McFLc: a runif is", runif(1), "\n")
##     nca2 <- oncoSimulSample(pops, fe, finalTime = ft,
##                          mutationPropGrowth = FALSE,
##                          initSize = no, 
##                          initMutant = "a", model = "McFL",
##                          onlyCancer = FALSE, detectionProb = NA, sampleEvery = 0.01,
##                           detectionSize = 1e9,
##                           detectionDrivers = 99,
##                           seed =NULL, max.wall.time = 3000,
##                           thresholdWhole = x)
##     gc(); 
##     
##     cat("\n ossmpc1McFLd: a runif is", runif(1), "\n")
##     ncb2 <- oncoSimulSample(pops, fe, finalTime = ft,
##                          mutationPropGrowth = FALSE,
##                          initSize = no, 
##                          initMutant = "b", model = "McFL",
##                          onlyCancer = FALSE, detectionProb = NA, sampleEvery = 0.01,
##                           detectionSize = 1e9,
##                           detectionDrivers = 99,
##                           seed =NULL, max.wall.time = 3000,
##                          thresholdWhole = x)
##     ## nca$popSummary[1:5, c(1:3, 8:9)]
##     ## ncb$popSummary[1:5, c(1:3, 8:9)]
##     ## nca2$popSummary[1:5, c(1:3, 8:9)]
##     ## ncb2$popSummary[1:5, c(1:3, 8:9)]
##     ## summary(nca$popSummary[, "NumClones"])
##     ## summary(ncb$popSummary[, "NumClones"])
##     ## summary(nca2$popSummary[, "NumClones"])
##     ## summary(ncb2$popSummary[, "NumClones"])
##     ## cat("\n mutsperclone\n")
##     ## summary(mutsPerCloneNCA <- rowSums(nca$popSample))
##     ## summary(mutsPerCloneNCB <- rowSums(ncb$popSample))
##     ## summary(mutsPerCloneNCA2 <- rowSums(nca2$popSample))
##     ## summary(mutsPerCloneNCB2 <- rowSums(ncb2$popSample))
##     ## (cc1 <- colSums(nca$popSample))
##     ## Because of how we do the "mutsPerClone", mutsPerClone almost the
##     ## same ans number of clones (except for the few cases when a clone
##     ## has gone extinct)
##     if(!is.array(nca$popSample)) {
##         ## occasionally, I get funny things
##         warning("nca$popSample not an array")
##         cat(class(nca$popSample))
##         cat(dim(nca$popSample))
##         cat(str(nca$popSample))
##     }
##     if(!is.array(ncb$popSample)) {
##         ## occasionally, I get funny things
##         warning("ncb$popSample not an array")
##         cat(class(ncb$popSample))
##         cat(dim(ncb$popSample))
##         cat(str(ncb$popSample))
##     }
##     mutsPerCloneNCA <- rowSums(nca$popSample)
##     mutsPerCloneNCB <- rowSums(ncb$popSample)
##     mutsPerCloneNCA2 <- rowSums(nca2$popSample)
##     mutsPerCloneNCB2 <- rowSums(ncb2$popSample)
##     ## median, because occasionally we can get something far out.
##     expect_true( median(mutsPerCloneNCA) >
##                  median(mutsPerCloneNCA2))
##     expect_true( median(mutsPerCloneNCB) >
##                  median(mutsPerCloneNCB2))
##     expect_true( median(mutsPerCloneNCA) >
##                  median(mutsPerCloneNCB))
##     expect_true( median(nca$popSummary[, "NumClones"]) >
##                  median(nca2$popSummary[, "NumClones"]))
##     expect_true( median(ncb$popSummary[, "NumClones"]) >
##                  median(ncb2$popSummary[, "NumClones"]))
##     expect_true( median(nca$popSummary[, "NumClones"]) >
##                  median(ncb$popSummary[, "NumClones"]))
## gc() 
## })
## cat("\n", date(), "\n")




## Commented out because it is about the same as 1 and 3, but here with
## tinier diffs which are harder to detect and can fail unless we use a
## huge number of pops. This is an overkill.

## cat("\n", date(), "\n")
## test_that("Ordering of number of clones and mutsPerClone with mutpropgrowth, 2", {
##     
##     
##     cat("\n mpc2: a runif is", runif(1), "\n")
##     ## The s coefficient is small, and so small differences between nca and
##     ## nca2.
##     ft <- 15 ## going beyond 16 or so, gets it to bail because of reaching max
##     ## pop
##     pops <- 400
##     lni <- 300
##     no <- 10
##     ni <- c(1, 0.5, rep(0, lni))
##     names(ni) <- c("a", "b", paste0("n", seq.int(lni)))
##     fe <- allFitnessEffects(noIntGenes = ni)
##     
##     
##     cat("\n mpc2a: a runif is", runif(1), "\n")
##     nca <- oncoSimulPop(pops, fe, finalTime = ft,
##                         mutationPropGrowth = TRUE,
##                         initSize = no, keepEvery = 1,
##                         initMutant = "a",
##                         onlyCancer = FALSE, detectionProb = NA, seed = NULL, mc.cores = 2)
##     gc()
##     
##     
##     cat("\n mpc2b: a runif is", runif(1), "\n")
##     ncb <- oncoSimulPop(pops, fe, finalTime = ft,
##                         mutationPropGrowth = TRUE,
##                         initSize = no, keepEvery = 1,
##                         initMutant = "b",
##                         onlyCancer = FALSE, detectionProb = NA, seed = NULL, mc.cores = 2)
##     gc()
##     
##     
##     cat("\n mpc2c: a runif is", runif(1), "\n")
##     nca2 <- oncoSimulPop(pops, fe, finalTime = ft,
##                          mutationPropGrowth = FALSE,
##                          initSize = no, keepEvery = 1,
##                          initMutant = "a",
##                          onlyCancer = FALSE, detectionProb = NA, seed = NULL, mc.cores = 2)
##     gc()
##     
##     
##     cat("\n mpc2d: a runif is", runif(1), "\n")
##     ncb2 <- oncoSimulPop(pops, fe, finalTime = ft,
##                          mutationPropGrowth = FALSE,
##                          initSize = no, keepEvery = 1,
##                          initMutant = "b",
##                          onlyCancer = FALSE, detectionProb = NA, seed = NULL, mc.cores = 2)
##     gc()
##     ## summary(nca)[, c(1, 2, 3, 8, 9)]
##     ## summary(nca2)[, c(1, 2, 3, 8, 9)]
##     ## summary(ncb)[, c(1, 2, 3, 8, 9)]
##     ## summary(ncb2)[, c(1, 2, 3, 8, 9)]
##     expect_true(var(summary(nca)$NumClones) > 1e-4)
##     expect_true(var(summary(ncb)$NumClones) > 1e-4)
##     expect_true(var(summary(nca2)$NumClones) > 1e-4)
##     expect_true(var(summary(ncb2)$NumClones) > 1e-4)
##     ## The real comparison
##     expect_true( median(summary(nca)$NumClones) >
##                  median(summary(ncb)$NumClones))
##     expect_true( median(summary(ncb)$NumClones) >
##                  median(summary(ncb2)$NumClones))
##     expect_true( median(mutsPerClone(nca)) >
##                  median(mutsPerClone(ncb)))
##     ## next can fail, as differences are small
##     expect_true( median(mutsPerClone(ncb)) >
##                  median(mutsPerClone(ncb2)))
##     ## These can fail in this case, since small diffs. as small mutlipliers
##     expect_true( mean(mutsPerClone(nca)) >
##                  mean(mutsPerClone(nca2)))
##     expect_true( median(summary(nca)$NumClones) >
##                  median(summary(nca2)$NumClones))
##     gc()
## })
