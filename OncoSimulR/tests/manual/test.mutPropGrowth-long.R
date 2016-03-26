cat(paste("\n Starting at mutPropGrowth long", date(), "\n"))

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




RNGkind("L'Ecuyer-CMRG") ## for the mclapplies
## If crashes I want to see where: thus output seed.



## this tests takes 10 seconds. Moved to long. So this was in the standard
## ones.
date()
test_that("mutPropGrowth diffs with s> 0", {
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n mgp1: the seed is", pseed, "\n")
    ft <- 4 ## 2.7
    pops <- 100
    lni <- 150 ## 100
    no <- 1e3 ## 5e1 
    ni <- c(2, rep(0, lni)) ## 2 ## 4 ## 5
    mu <- 1e-5 ## 1e-6
    names(ni) <- c("a", paste0("n", seq.int(lni)))
    ni <- sample(ni) ## scramble 
    fe <- allFitnessEffects(noIntGenes = ni)
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n mpg1a: the seed is", pseed, "\n")
    nca <- oncoSimulPop(pops, fe, finalTime = ft,
                        mu = mu,
                        mutationPropGrowth = TRUE,
                        initSize = no, sampleEvery = 0.1,
                        initMutant = "a", keepEvery = 1,
                        onlyCancer = FALSE, seed = NULL, mc.cores = 2)
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n mpg1c: the seed is", pseed, "\n")
    nca2 <- oncoSimulPop(pops, fe, finalTime = ft,
                         mu = mu,
                        mutationPropGrowth = FALSE,
                        initSize = no, sampleEvery = 0.1,
                        initMutant = "a", keepEvery = 1,
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
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n mgp1ND: the seed is", pseed, "\n")
    ft <- 4 ## 2.7
    pops <- 200
    lni <- 150 ## 100
    no <- 1e3 ## 5e1 
    ni <- c(0, rep(0, lni)) ## 2 ## 4 ## 5
    mu <- 1e-5 ## 1e-6
    names(ni) <- c("a", paste0("n", seq.int(lni)))
    ni <- sample(ni) ## scramble 
    fe <- allFitnessEffects(noIntGenes = ni)
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n mpg1NDa: the seed is", pseed, "\n")
    nca <- oncoSimulPop(pops, fe, finalTime = ft,
                        mu = mu,
                        mutationPropGrowth = TRUE,
                        initSize = no, sampleEvery = 0.1,
                        initMutant = "a", keepEvery = 1,
                        onlyCancer = FALSE, seed = NULL, mc.cores = 2)
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n mpg1NDc: the seed is", pseed, "\n")
    nca2 <- oncoSimulPop(pops, fe, finalTime = ft,
                         mu = mu,
                        mutationPropGrowth = FALSE,
                        initSize = no, sampleEvery = 0.1,
                        initMutant = "a", keepEvery = 1,
                        onlyCancer = FALSE, seed = NULL, mc.cores = 2)
    ## summary(nca)[1:20, c(1, 2, 3, 8, 9)]
    ## summary(nca2)[1:20, c(1, 2, 3, 8, 9)]
    p.fail <- 1e-3
    expect_true(t.test(sqrt(summary(nca)$NumClones),
                       sqrt(summary(nca2)$NumClones))$p.value > p.fail)
    expect_true(t.test(sqrt(mutsPerClone(nca)),
                       sqrt(mutsPerClone(nca2)))$p.value > p.fail)
    ## I once saw a weird thing
    ## expect_true(var(summary(nca)$NumClones) > 1e-4)
    ## expect_true(var(summary(nca2)$NumClones) > 1e-4)
    ## summary(summary(nca)$NumClones)
    ## summary(summary(nca2)$NumClones)
    ## summary(mutsPerClone(nca))
    ## summary(mutsPerClone(nca2))
})
cat("\n", date(), "\n")



date()
test_that("mutPropGrowth no diffs with s = 0, McFL", {
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n mcfND1: the seed is", pseed, "\n")
    ft <- 3 
    pops <- 200
    lni <- 100
    no <- 1e3 ## 5e1 
    ni <- c(0, rep(0, lni)) ## 5
    names(ni) <- c("a", paste0("n", seq.int(lni)))
    ni <- sample(ni) ## scramble
    fe <- allFitnessEffects(noIntGenes = ni)
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n mcfND1a: the seed is", pseed, "\n")
    nca <- oncoSimulPop(pops, fe, finalTime = ft,
                        mutationPropGrowth = TRUE,
                        initSize = no,
                        initMutant = "a", model = "McFL",
                        onlyCancer = FALSE, seed = NULL, mc.cores = 2)
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n mcfND1c: the seed is", pseed, "\n")
    nca2 <- oncoSimulPop(pops, fe, finalTime = ft,
                        mutationPropGrowth = FALSE,
                        initSize = no,
                        initMutant = "a", model = "McFL",
                        onlyCancer = FALSE, seed = NULL, mc.cores = 2)
    ## summary(nca)[1:20, c(1, 2, 3, 8, 9)]
    ## summary(nca2)[1:20, c(1, 2, 3, 8, 9)]
    ## summary(summary(nca)$NumClones)
    ## summary(summary(nca2)$NumClones)
    ## summary(mutsPerClone(nca))
    ## summary(mutsPerClone(nca2))
    p.fail <- 1e-3
    expect_true(t.test(sqrt(summary(nca)$NumClones),
                       sqrt(summary(nca2)$NumClones))$p.value > p.fail)
    expect_true(t.test(sqrt(mutsPerClone(nca)),
                       sqrt(mutsPerClone(nca2)))$p.value > p.fail)
})
cat("\n", date(), "\n")




date()
test_that("mutPropGrowth no diffs with s = 0, oncoSimulSample", {
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n oss1: the seed is", pseed, "\n")
    ft <- 3.5 ## 4
    pops <- 200
    lni <- 200 ## 150
    no <- 1e3 ## 5e1 
    ni <- c(0, rep(0, lni)) ## 2 ## 4 ## 5
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
                        initSize = no, sampleEvery = 0.1,
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
                        initSize = no, sampleEvery = 0.1,
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
    p.fail <- 1e-3
    expect_true(t.test(sqrt( nca$popSummary[, "NumClones"]  ),
                       sqrt( nca2$popSummary[, "NumClones"]  ))$p.value > p.fail)
    expect_true(t.test(sqrt(mutsPerClone1),
                       sqrt(mutsPerClone2))$p.value > p.fail)
})


date()
test_that("mutPropGrowth no diffs with s = 0, oncoSimulSample, McFL", {
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n ossmcfND1: the seed is", pseed, "\n")
    ft <- 40 ## 4
    pops <- 200
    lni <- 200 ## 150
    no <- 1e3 ## 5e1 
    ni <- c(0, rep(0, lni)) ## 2 ## 4 ## 5
    mu <- 5e-7 ## 1e-6
    x <- 1e-9
    names(ni) <- c("a", paste0("n", seq.int(lni)))
    ni <- sample(ni) ## scramble
    fe <- allFitnessEffects(noIntGenes = ni)
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n ossmcfND1a: the seed is", pseed, "\n")
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
    cat("\n ossmcfND1c: the seed is", pseed, "\n")
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
    p.fail <- 1e-3
    expect_true(t.test(sqrt( nca$popSummary[, "NumClones"]  ),
                       sqrt( nca2$popSummary[, "NumClones"]  ))$p.value > p.fail)
    expect_true(t.test(sqrt(mutsPerClone1),
                       sqrt(mutsPerClone2))$p.value > p.fail)  
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


## This generally works, but not always. Because: a) starting with a you
## can get a mutation in b. Or starting with b you can get a mutation in
## a. When either happens is stochastic. And we are also mixing the
## mutation proportional to birth with the simple increase in clones due
## to larger pop sizes. So this ain't a good test but here finalTime is
## very short so unlikely the "hit a and hit b".


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

## Psss ideas.  The idea here is that you hit the module with fitness
## effects, that sets things to grow, and then you either have or not
## mutation proportional to growth.

cat("\n", date(), "\n")
test_that("Without initmutant", {
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n s3: the seed is", pseed, "\n")
    pops <- 200
    lni <- 1 ## no fitness effects genes
    fni <- 50 ## fitness effects genes
    no <- 1e3
    ft <- 5
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
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n s3a: the seed is", pseed, "\n")
    s3.ng <- oncoSimulPop(pops,
                          f3,
                          mu = mu,
                          mutationPropGrowth = FALSE,
                          finalTime =ft,
                          initSize = no,
                          onlyCancer = FALSE,
                          seed = NULL, mc.cores = 2)
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n s3b: the seed is", pseed, "\n")
    s3.g <- oncoSimulPop(pops,
                         f3,
                         mu = mu,
                         mutationPropGrowth = TRUE,
                         finalTime =ft,
                         initSize = no,
                         onlyCancer = FALSE,
                         seed = NULL, mc.cores = 2)
    ## summary(s3.g)[, c(1, 2, 3, 8, 9)]
    ## summary(s3.ng)[, c(1, 2, 3, 8, 9)]
    expect_true( mean(mutsPerClone(s3.g)) >
                 mean(mutsPerClone(s3.ng)))
    expect_true( median(summary(s3.g)$NumClones) >
                 median(summary(s3.ng)$NumClones))
gc() 
})
cat("\n", date(), "\n")

cat("\n", date(), "\n")
test_that("Without initmutant, 2", {
    ## More of the above. Use smaller s2 and smaller mutation, but then to
    ## see it reliably you need large ft and we also increase
    ## init. pop. size.
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n s2: the seed is", pseed, "\n")
    s2 <- 1.0
    ft <- 15
    pops <- 200
    lni <- 1 ## no fitness effects genes
    fni <- 50 ## fitness effects genes
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
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n s2a: the seed is", pseed, "\n")
    s2.ng <- oncoSimulPop(pops,
                          f2,
                          mu = mu,
                          mutationPropGrowth = FALSE,
                          finalTime =ft,
                          initSize = no, keepEvery = 1,
                          onlyCancer = FALSE,
                          seed = NULL, mc.cores = 2)
    gc(); pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n s2b: the seed is", pseed, "\n")
    s2.g <- oncoSimulPop(pops,
                         f2,
                         mu = mu,
                         mutationPropGrowth = TRUE,
                         finalTime =ft,
                         initSize = no, keepEvery = 1,
                         onlyCancer = FALSE,
                         seed = NULL, mc.cores = 2)
    ## summary(s2.g)[, c(1, 2, 3, 8, 9)]
    ## summary(s2.ng)[, c(1, 2, 3, 8, 9)]
    expect_true( mean(mutsPerClone(s2.g)) >
                 mean(mutsPerClone(s2.ng)))
    expect_true( median(summary(s2.g)$NumClones) >
                 median(summary(s2.ng)$NumClones))
gc() 
})
cat("\n", date(), "\n")

cat("\n", date(), "\n")
test_that("McFL: Without initmutant", {
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n mcfls2: the seed is", pseed, "\n")
    s2 <- 2.0
    ft <- 250
    pops <- 200 ## 200
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
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n mcfls2a: the seed is", pseed, "\n")
    s2.ng <- oncoSimulPop(pops,
                          f2,
                          mu = mu,
                          mutationPropGrowth = FALSE,
                          finalTime =ft,
                          initSize = no, keepEvery = 5,
                          onlyCancer = FALSE, model = "McFL",
                          seed = NULL, mc.cores = 2)
    gc(); pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n mcfls2b: the seed is", pseed, "\n")
    s2.g <- oncoSimulPop(pops,
                         f2,
                         mu = mu,
                         mutationPropGrowth = TRUE,
                         finalTime =ft,
                         initSize = no, keepEvery = 5, 
                         onlyCancer = FALSE, model = "McFL",
                         seed = NULL, mc.cores = 2)
    ## summary(s2.g)[, c(1, 2, 3, 8, 9)]
    ## summary(s2.ng)[, c(1, 2, 3, 8, 9)]
    expect_true( mean(mutsPerClone(s2.g)) >
                 mean(mutsPerClone(s2.ng)))
    expect_true( median(summary(s2.g)$NumClones) >
                 median(summary(s2.ng)$NumClones))
    ## summary(mutsPerClone(s2.g))
    ## summary(mutsPerClone(s2.ng))
    ## summary(summary(s2.g)$NumClones)
    ## summary(summary(s2.ng)$NumClones)
gc() 
})
cat("\n", date(), "\n")


cat(paste("\n Ending mutPropGrwoth-long at", date(), "\n"))









## ## This generally works, but not always. Because: a) starting with a you
## ## can get a mutation in b. Or starting with b you can get a mutation in
## ## a. When either happens is stochastic. And we are also mixing the
## ## mutation proportional to birth with the simple increase in clones due
## ## to larger pop sizes. So this ain't a good test. When you start, say,
## ## from b and you hit a soon, you can get HUGE growth rates and that
## ## induces huge variability, slow times, etc.

## cat("\n", date(), "\n")
## test_that("Ordering of number of clones and mutsPerClone with mutpropgrowth, 1", {
##     pseed <- sample(9999999, 1)
##     set.seed(pseed)
##     cat("\n mpc1: the seed is", pseed, "\n")
##     ft <- 2.5
##     pops <- 200
##     lni <- 500 ## with, say, 40 or a 100, sometimes fails the comparisons
##                ## with small differences.
##     no <- 10
##     ni <- c(5, 3, rep(0, lni))
##     names(ni) <- c("a", "b", paste0("n", seq.int(lni)))
##     fe <- allFitnessEffects(noIntGenes = ni)
##     pseed <- sample(9999999, 1)
##     set.seed(pseed)
##     cat("\n mpc1a: the seed is", pseed, "\n")
##     nca <- oncoSimulPop(pops, fe, finalTime = ft,
##                         mutationPropGrowth = TRUE,
##                         initSize = no,
##                         initMutant = "a",
##                         onlyCancer = FALSE, seed = NULL, mc.cores = 2)
##     gc()
##     pseed <- sample(9999999, 1)
##     set.seed(pseed)
##     cat("\n mpc1b: the seed is", pseed, "\n")
##     ncb <- oncoSimulPop(pops, fe, finalTime = ft,
##                         mutationPropGrowth = TRUE,
##                         initSize = no,
##                         initMutant = "b",
##                         onlyCancer = FALSE, seed = NULL, mc.cores = 2)
##     gc()
##     pseed <- sample(9999999, 1)
##     set.seed(pseed)
##     cat("\n mpc1c: the seed is", pseed, "\n")
##     nca2 <- oncoSimulPop(pops, fe, finalTime = ft,
##                          mutationPropGrowth = FALSE,
##                          initSize = no,
##                          initMutant = "a",
##                          onlyCancer = FALSE, seed = NULL, mc.cores = 2)
##     pseed <- sample(9999999, 1)
##     set.seed(pseed)
##     cat("\n mpc1d: the seed is", pseed, "\n")
##     ncb2 <- oncoSimulPop(pops, fe, finalTime = ft,
##                          mutationPropGrowth = FALSE,
##                          initSize = no,
##                          initMutant = "b",
##                          onlyCancer = FALSE, seed = NULL, mc.cores = 2)
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
##     pseed <- sample(9999999, 1)
##     set.seed(pseed)
##     cat("\n mpc3: the seed is", pseed, "\n")
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
##     pseed <- sample(9999999, 1)
##     set.seed(pseed)
##     cat("\n mpc3a: the seed is", pseed, "\n")
##     nca <- oncoSimulPop(pops, fe, finalTime = ft,
##                         mutationPropGrowth = TRUE,
##                         mu = mu,
##                         initSize = no, keepEvery = 1,
##                         initMutant = "a",
##                         onlyCancer = FALSE, seed = NULL, mc.cores = 2)
##     gc()
##     pseed <- sample(9999999, 1)
##     set.seed(pseed)
##     cat("\n mpc3b: the seed is", pseed, "\n")
##     ncb <- oncoSimulPop(pops, fe, finalTime = ft,
##                         mutationPropGrowth = TRUE,
##                         mu = mu,
##                         initSize = no, keepEvery = 1,
##                         initMutant = "b",
##                         onlyCancer = FALSE, seed = NULL, mc.cores = 2)
##     gc()
##     pseed <- sample(9999999, 1)
##     set.seed(pseed)
##     cat("\n mpc3c: the seed is", pseed, "\n")
##     nca2 <- oncoSimulPop(pops, fe, finalTime = ft,
##                          mutationPropGrowth = FALSE,
##                          mu = mu,
##                          initSize = no, keepEvery = 1,
##                          initMutant = "a",
##                          onlyCancer = FALSE, seed = NULL, mc.cores = 2)
##     gc()
##     pseed <- sample(9999999, 1)
##     set.seed(pseed)
##     cat("\n mpc3d: the seed is", pseed, "\n")
##     ncb2 <- oncoSimulPop(pops, fe, finalTime = ft,
##                          mutationPropGrowth = FALSE,
##                          mu = mu,                     
##                          initSize = no, keepEvery = 1,
##                          initMutant = "b",
##                          onlyCancer = FALSE, seed = NULL, mc.cores = 2)
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
##     pseed <- sample(9999999, 1)
##     set.seed(pseed)
##     cat("\n mpcmcf1: the seed is", pseed, "\n")
##     ft <- 20 ## unless large you rarely get triple, etc, mutatns
##     pops <- 200
##     lni <- 50 
##     no <- 1e3
##     ni <- c(3, 1.5, rep(0, lni))
##     names(ni) <- c("a", "b", paste0("n", seq.int(lni)))
##     fe <- allFitnessEffects(noIntGenes = ni)
##     pseed <- sample(9999999, 1)
##     set.seed(pseed)
##     cat("\n mpcm1a: the seed is", pseed, "\n")
##     nca <- oncoSimulPop(pops, fe, finalTime = ft,
##                         mutationPropGrowth = TRUE,
##                         initSize = no, model = "McFL",
##                         initMutant = "a", keepEvery = 1,
##                         onlyCancer = FALSE, seed = NULL, mc.cores = 2)
##     gc()
##     pseed <- sample(9999999, 1)
##     set.seed(pseed)
##     cat("\n mpcm1b: the seed is", pseed, "\n")
##     ncb <- oncoSimulPop(pops, fe, finalTime = ft,
##                         mutationPropGrowth = TRUE,
##                         initSize = no, model = "McFL",
##                         initMutant = "b", keepEvery = 1,
##                         onlyCancer = FALSE, seed = NULL, mc.cores = 2)
##     gc()
##     pseed <- sample(9999999, 1)
##     set.seed(pseed)
##     cat("\n mpcm1c: the seed is", pseed, "\n")
##     nca2 <- oncoSimulPop(pops, fe, finalTime = ft,
##                          mutationPropGrowth = FALSE,
##                          initSize = no, model = "McFL",
##                          initMutant = "a", keepEvery = 1,
##                          onlyCancer = FALSE, seed = NULL, mc.cores = 2)
##     gc()
##     pseed <- sample(9999999, 1)
##     set.seed(pseed)
##     cat("\n mpcm1d: the seed is", pseed, "\n")
##     ncb2 <- oncoSimulPop(pops, fe, finalTime = ft,
##                          mutationPropGrowth = FALSE,
##                          initSize = no, model = "McFL",
##                          initMutant = "b", keepEvery = 1,
##                          onlyCancer = FALSE, seed = NULL, mc.cores = 2)
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
##     pseed <- sample(9999999, 1)
##     set.seed(pseed)
##     cat("\n mpcmcf2: the seed is", pseed, "\n")
##     ## Increase ft
##     ft <- 350 
##     pops <- 150
##     lni <- 30
##     no <- 1e3
##     ni <- c(2, 0.8, rep(0, lni))
##     names(ni) <- c("a", "b", paste0("n", seq.int(lni)))
##     fe <- allFitnessEffects(noIntGenes = ni)
##     pseed <- sample(9999999, 1)
##     set.seed(pseed)
##     cat("\n mpcmcf2a: the seed is", pseed, "\n")
##     nca <- oncoSimulPop(pops, fe, finalTime = ft,
##                         mutationPropGrowth = TRUE,
##                         initSize = no, model = "McFL",
##                         initMutant = "a", keepEvery = 1,
##                         onlyCancer = FALSE, seed = NULL, mc.cores = 2)
##     gc()
##     pseed <- sample(9999999, 1)
##     set.seed(pseed)
##     cat("\n mpcmcf2b: the seed is", pseed, "\n")
##     ncb <- oncoSimulPop(pops, fe, finalTime = ft,
##                         mutationPropGrowth = TRUE,
##                         initSize = no, model = "McFL",
##                         initMutant = "b", keepEvery = 1,
##                         onlyCancer = FALSE, seed = NULL, mc.cores = 2)
##     gc()
##     pseed <- sample(9999999, 1)
##     set.seed(pseed)
##     cat("\n mpcmcf2c: the seed is", pseed, "\n")
##     nca2 <- oncoSimulPop(pops, fe, finalTime = ft,
##                          mutationPropGrowth = FALSE,
##                          initSize = no, model = "McFL",
##                          initMutant = "a", keepEvery = 1,
##                          onlyCancer = FALSE, seed = NULL, mc.cores = 2)
##     gc()
##     pseed <- sample(9999999, 1)
##     set.seed(pseed)
##     cat("\n mpcmcf2d: the seed is", pseed, "\n")
##     ncb2 <- oncoSimulPop(pops, fe, finalTime = ft,
##                          mutationPropGrowth = FALSE,
##                          initSize = no, model = "McFL",
##                          initMutant = "b", keepEvery = 1,
##                          onlyCancer = FALSE, seed = NULL, mc.cores = 2)
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
##     pseed <- sample(9999999, 1)
##     set.seed(pseed)
##     cat("\n mpcmcf3: the seed is", pseed, "\n")
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
##     pseed <- sample(9999999, 1)
##     set.seed(pseed)
##     cat("\n mpcmcf3a: the seed is", pseed, "\n")
##     nca <- oncoSimulPop(pops, f1, finalTime = ft,
##                         mutationPropGrowth = TRUE,
##                         mu = mu, keepEvery = 1,
##                         initSize = no, model = "McFL",
##                         initMutant = "a1", 
##                         onlyCancer = FALSE, seed = NULL, mc.cores = 2)
##     pseed <- sample(9999999, 1)
##     gc()
##     set.seed(pseed)
##     cat("\n mpcmcf3b: the seed is", pseed, "\n")
##     ncb <- oncoSimulPop(pops, f1, finalTime = ft,
##                         mutationPropGrowth = TRUE,
##                         mu = mu, keepEvery = 1,
##                         initSize = no, model = "McFL",
##                         initMutant = "b1",
##                         onlyCancer = FALSE, seed = NULL, mc.cores = 2)
##     ## OK, but it also happens below just because init a1 eventually grows
##     ## faster, so larger pop, so more mutants, etc
##     ## We could try just counting the number of mutation events, but we would
##     ## need to standardize by total time of life over all cells.
##     ## As growth rates are faster (those that start with a, here, and then
##     ## get a b), even when mutationPropGrowth is false, the population
##     ## gets larger. And thus, it is easier to get a clone with three
##     ## mutations, etc.
##     gc()
##     pseed <- sample(9999999, 1)
##     set.seed(pseed)
##     cat("\n mpcmcf3c: the seed is", pseed, "\n")
##     nca2 <- oncoSimulPop(pops, f1, finalTime = ft,
##                          mutationPropGrowth = FALSE,
##                          mu = mu, keepEvery = 1,
##                          initSize = no, model = "McFL",
##                          initMutant = "a1",
##                          onlyCancer = FALSE, seed = NULL, mc.cores = 2)
##     gc(); pseed <- sample(9999999, 1)
##     set.seed(pseed)
##     cat("\n mpcmcf3d: the seed is", pseed, "\n")
##     ncb2 <- oncoSimulPop(pops, f1, finalTime = ft,
##                          mutationPropGrowth = FALSE,
##                          mu = mu, keepEvery = 1,                    
##                          initSize = no, model = "McFL",
##                          initMutant = "b1",
##                          onlyCancer = FALSE, seed = NULL, mc.cores = 2)
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
##     pseed <- sample(9999999, 1)
##     set.seed(pseed)
##     cat("\n ossmpc1: the seed is", pseed, "\n")
##     ft <- 2  
##     pops <- 400
##     lni <- 500 ## with, say, 40 or a 100, sometimes fails the comparisons
##                ## with small differences.
##     no <- 10
##     x <- 1e-40
##     ni <- c(5, 3, rep(0, lni))
##     names(ni) <- c("a", "b", paste0("n", seq.int(lni)))
##     fe <- allFitnessEffects(noIntGenes = ni)
##     pseed <- sample(9999999, 1)
##     set.seed(pseed)
##     cat("\n ossmpc1a: the seed is", pseed, "\n")
##     nca <- oncoSimulSample(pops, fe, finalTime = ft,
##                         mutationPropGrowth = TRUE,
##                         initSize = no, 
##                         initMutant = "a",
##                         onlyCancer = FALSE,  sampleEvery = 0.01,
##                           detectionSize = 1e9,
##                           detectionDrivers = 99,
##                           seed =NULL, max.wall.time = 3000,
##                           thresholdWhole = x)
##     gc(); pseed <- sample(9999999, 1)
##     set.seed(pseed)
##     cat("\n ossmpc1b: the seed is", pseed, "\n")
##     ncb <- oncoSimulSample(pops, fe, finalTime = ft,
##                         mutationPropGrowth = TRUE,
##                         initSize = no, 
##                         initMutant = "b",
##                         onlyCancer = FALSE, sampleEvery = 0.01,
##                           detectionSize = 1e9,
##                           detectionDrivers = 99,
##                           seed =NULL, max.wall.time = 3000,
##                           thresholdWhole = x)
##     gc(); pseed <- sample(9999999, 1)
##     set.seed(pseed)
##     cat("\n ossmpc1c: the seed is", pseed, "\n")
##     nca2 <- oncoSimulSample(pops, fe, finalTime = ft,
##                          mutationPropGrowth = FALSE,
##                          initSize = no, 
##                          initMutant = "a",
##                          onlyCancer = FALSE, sampleEvery = 0.01,
##                           detectionSize = 1e9,
##                           detectionDrivers = 99,
##                           seed =NULL, max.wall.time = 3000,
##                           thresholdWhole = x)
##     gc(); pseed <- sample(9999999, 1)
##     set.seed(pseed)
##     cat("\n ossmpc1d: the seed is", pseed, "\n")
##     ncb2 <- oncoSimulSample(pops, fe, finalTime = ft,
##                          mutationPropGrowth = FALSE,
##                          initSize = no, 
##                          initMutant = "b",
##                          onlyCancer = FALSE, sampleEvery = 0.01,
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
##     pseed <- sample(9999999, 1)
##     set.seed(pseed)
##     cat("\n ossmpc1McFL: the seed is", pseed, "\n")
##     ft <- 2  
##     pops <- 200 ## 200
##     lni <- 500 ## with, say, 40 or a 100, sometimes fails the comparisons
##                ## with small differences.
##     no <- 100
##     x <- 1e-40
##     ni <- c(5, 3, rep(0, lni))
##     names(ni) <- c("a", "b", paste0("n", seq.int(lni)))
##     fe <- allFitnessEffects(noIntGenes = ni)
##     pseed <- sample(9999999, 1)
##     set.seed(pseed)
##     cat("\n ossmpc1McFLa: the seed is", pseed, "\n")
##     nca <- oncoSimulSample(pops, fe, finalTime = ft,
##                         mutationPropGrowth = TRUE,
##                         initSize = no, 
##                         initMutant = "a", model = "McFL",
##                         onlyCancer = FALSE,  sampleEvery = 0.01,
##                           detectionSize = 1e9,
##                           detectionDrivers = 99,
##                           seed =NULL, max.wall.time = 3000,
##                           thresholdWhole = x)
##     gc(); pseed <- sample(9999999, 1)
##     set.seed(pseed)
##     cat("\n ossmpc1McFLb: the seed is", pseed, "\n")
##     ncb <- oncoSimulSample(pops, fe, finalTime = ft,
##                         mutationPropGrowth = TRUE,
##                         initSize = no, 
##                         initMutant = "b", model = "McFL",
##                         onlyCancer = FALSE, sampleEvery = 0.01,
##                           detectionSize = 1e9,
##                           detectionDrivers = 99,
##                           seed =NULL, max.wall.time = 3000,
##                           thresholdWhole = x)
##     gc(); pseed <- sample(9999999, 1)
##     set.seed(pseed)
##     cat("\n ossmpc1McFLc: the seed is", pseed, "\n")
##     nca2 <- oncoSimulSample(pops, fe, finalTime = ft,
##                          mutationPropGrowth = FALSE,
##                          initSize = no, 
##                          initMutant = "a", model = "McFL",
##                          onlyCancer = FALSE, sampleEvery = 0.01,
##                           detectionSize = 1e9,
##                           detectionDrivers = 99,
##                           seed =NULL, max.wall.time = 3000,
##                           thresholdWhole = x)
##     gc(); pseed <- sample(9999999, 1)
##     set.seed(pseed)
##     cat("\n ossmpc1McFLd: the seed is", pseed, "\n")
##     ncb2 <- oncoSimulSample(pops, fe, finalTime = ft,
##                          mutationPropGrowth = FALSE,
##                          initSize = no, 
##                          initMutant = "b", model = "McFL",
##                          onlyCancer = FALSE, sampleEvery = 0.01,
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
##     pseed <- sample(9999999, 1)
##     set.seed(pseed)
##     cat("\n mpc2: the seed is", pseed, "\n")
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
##     pseed <- sample(9999999, 1)
##     set.seed(pseed)
##     cat("\n mpc2a: the seed is", pseed, "\n")
##     nca <- oncoSimulPop(pops, fe, finalTime = ft,
##                         mutationPropGrowth = TRUE,
##                         initSize = no, keepEvery = 1,
##                         initMutant = "a",
##                         onlyCancer = FALSE, seed = NULL, mc.cores = 2)
##     gc()
##     pseed <- sample(9999999, 1)
##     set.seed(pseed)
##     cat("\n mpc2b: the seed is", pseed, "\n")
##     ncb <- oncoSimulPop(pops, fe, finalTime = ft,
##                         mutationPropGrowth = TRUE,
##                         initSize = no, keepEvery = 1,
##                         initMutant = "b",
##                         onlyCancer = FALSE, seed = NULL, mc.cores = 2)
##     gc()
##     pseed <- sample(9999999, 1)
##     set.seed(pseed)
##     cat("\n mpc2c: the seed is", pseed, "\n")
##     nca2 <- oncoSimulPop(pops, fe, finalTime = ft,
##                          mutationPropGrowth = FALSE,
##                          initSize = no, keepEvery = 1,
##                          initMutant = "a",
##                          onlyCancer = FALSE, seed = NULL, mc.cores = 2)
##     gc()
##     pseed <- sample(9999999, 1)
##     set.seed(pseed)
##     cat("\n mpc2d: the seed is", pseed, "\n")
##     ncb2 <- oncoSimulPop(pops, fe, finalTime = ft,
##                          mutationPropGrowth = FALSE,
##                          initSize = no, keepEvery = 1,
##                          initMutant = "b",
##                          onlyCancer = FALSE, seed = NULL, mc.cores = 2)
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
