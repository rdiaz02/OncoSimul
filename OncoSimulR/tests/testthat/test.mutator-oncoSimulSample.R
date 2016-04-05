## Repeat tests in test.mutator, using oncoSimulSample.
## This is a concession to extreme paranoia.


cat(paste("\n Starting test.mutator-oncoSimulSample.R test at", date()))
RNGkind("L'Ecuyer-CMRG") ## for the mclapplies

## require(car) ## for linearHypothesis, below. In the namespace

enom <- function(name, mu, ni = no, pp = pops) {
    ## expected without a given name for init
    ii <- which(names(mu) == name)
    out <- ni * pp * mu[-ii]
    out[order(names(out))]
}

pnom <- function(name, mu, ni = no, pp = pops) {
    ee <- enom(name, mu, ni, pp)
    ee/sum(ee)
}

snomSampl <- function(name, out) {
    ## observed without the init
    cs <- colSums(out$popSample)
    ii <- which(names(cs) == name)
    cs <- cs[-ii]
    cs[order(names(cs))]
}

smSampl <- function(name, out) {
    ## totals for a given gene
    cs <- colSums(out$popSample)
    ii <- which(names(cs) == name)
    cs[ii]
}

totalindSampl <- function(out) {
    ## total num indivs
  sum(out$popSummary$TotalPopSize)  
}

medianNClonesOSS <- function(x) {
    median(x$popSummary[, "NumClones"])
}



## ugly hack. Of course, not really mutations per clone. But the closest iwth oncoSimulSample and sampling whole pop.
mutsPerCloneOSS <- function(out) {
    rowSums(out$popSample)
}

## These next two tests are probably the strongest (provided we accept
## using the initMutant) as we compare observed with expected and the
## estimated effect of mutator

## FIXME
##FIXME: Think this test: is it failing more than it should? Implement a two or three loop approach?

date()
test_that("Mutator increases by given factor with per-gene-mut rates: major axis and chi-sq test", {

    ## Two cases: mutator and no mutator, with variable mutation rates.
    ## rates such that rates of no mutator = rates of mutator * mutator.

    ## Why not compare mutlitplication factor keeping mutation rates
    ## constant? Because specially with mutator and large diffs in mut
    ## rates, with oncoSimulSample you undersample variation with
    ## wholePop, etc.

    ## Setings similar to oss11 in per-gene-mutation-rates but with the mutator
    
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n AEu8: the seed is", pseed, "\n")
    pops <- 2000
    ft <- 5e-3
    lni <- 7
    no <- 5e5
    ni <- c(0, 0, 0, rep(0, lni))
    ## scramble around names
    names(ni) <- c("hereisoneagene",
                   "oreoisasabgene",
                   "nnhsisthecgene",
                   replicate(lni,
                             paste(sample(letters, 12), collapse = "")))
    ni <- ni[order(names(ni))]
    fe <- allFitnessEffects(noIntGenes = ni)
    mutator1 <- rep(1, lni + 3)
    pg1 <- seq(from = 1e-9, to = 1e-6, length.out = lni + 3) ## max should not be
                                                  ## huge here as mutator
                                                  ## is 34. Can get beyond
                                                  ## 1
    names(mutator1) <- sample(names(ni))
    names(pg1) <- sample(names(ni))
    mutator1["oreoisasabgene"] <- 100
    m1 <- allMutatorEffects(noIntGenes = mutator1)
    ## pg1["hereisoneagene"] <- 1e-4 ## if this gets huge, then you are
    ##                               ## undersampling and the chi-square will
    ##                               ## fail. But then, we probably are
    ##                               ## running into numerical issues: 3
    ##                               ## orders of magnitude differences.
    m1.pg1.b <- oncoSimulSample(pops,
                           fe,
                           mu = pg1,
                           muEF = m1,
                           finalTime = ft,
                           mutationPropGrowth = FALSE,
                           initSize = no,
                           initMutant ="oreoisasabgene",
                           sampleEvery = 0.01, thresholdWhole = 1e-20,
                           detectionSize = 1e9,
                           detectionDrivers = 9999,
                           onlyCancer = FALSE, seed = NULL)
    ## m1.pg1.b$popSummary[, c(1:3, 8:9)]
    summary(m1.pg1.b$popSummary[, "NumClones"])
    ## Recall that init-mutant tests check always present of initMutant
    ## against a thresholWhole of 1. Here it is slightly different.
    expect_true(smSampl("oreoisasabgene", m1.pg1.b) == pops)
    ## catch a pattern that would make the previous trivially true
    expect_false(sum(m1.pg1.b$popSample) == pops * (lni + 3))
    ## next two, to compare with oss1a
    sort(enom("oreoisasabgene", pg1, no, pops))
    sort(snomSampl("oreoisasabgene", m1.pg1.b))
    ## Compare with the expected for this scenario
    p.fail <- 1e-3
    expect_true(chisq.test(snomSampl("oreoisasabgene", m1.pg1.b),
                           p = pnom("oreoisasabgene", pg1, no, pops))$p.value > p.fail)

    
    pg2 <- seq(from = 1e-7, to = 1e-4, length.out = lni + 3)
    names(pg2) <- names(pg1)
    m1.pg2.b <- oncoSimulSample(pops,
                           fe,
                           mu = pg2,
                           finalTime = ft,
                           mutationPropGrowth = FALSE,
                           initSize = no,
                           initMutant ="oreoisasabgene",
                           sampleEvery = 0.01, thresholdWhole = 1e-20,
                           detectionSize = 1e9,
                           detectionDrivers = 9999,
                           onlyCancer = FALSE, seed = NULL)
    ## m1.pg2.b$popSummary[, c(1:3, 8:9)]
    summary(m1.pg2.b$popSummary[, "NumClones"])
    ## Recall that init-mutant tests check always present of initMutant
    ## against a thresholWhole of 1. Here it is slightly different.
    expect_true(smSampl("oreoisasabgene", m1.pg2.b) == pops)
    ## catch a pattern that would make the previous trivially true
    expect_false(sum(m1.pg2.b$popSample) == pops * (lni + 3))
    ## next two, to compare with oss1a
    sort(enom("oreoisasabgene", pg2, no, pops))
    sort(snomSampl("oreoisasabgene", m1.pg2.b))
    p.fail <- 1e-3
    expect_true(chisq.test(snomSampl("oreoisasabgene", m1.pg2.b),
                           p = pnom("oreoisasabgene", pg2, no, pops))$p.value > p.fail)


    ## Compare the mutator with the no mutator
    expect_true(chisq.test(snomSampl("oreoisasabgene", m1.pg1.b),
                           snomSampl("oreoisasabgene", m1.pg2.b))$p.value > p.fail)
    y <- sqrt(snomSampl("oreoisasabgene", m1.pg1.b))
    x <- sqrt(snomSampl("oreoisasabgene", m1.pg2.b))
    mma <- smatr::ma(y ~ x, slope.test = 1, elev.test = 0) ## From smatr package, for major axis
    ## intercept not different from 0
    expect_true(mma$elevtest[[1]]$p > p.fail)
    expect_true(mma$slopetest[[1]]$p > p.fail)
    
    ## We could use a lm and do a simultaneous test on both slope and
    ## intercept as. But this is really asking for major axis regression

    ## lm1 <- lm(snomSampl("oreoisasabgene", m1.pg1.b) ~
    ##               snomSampl("oreoisasabgene", m1.pg2.b))
    ## ## test intercept is 0, slope is 1. Not technically fully correct, as
    ## ## X variable has noise. We should do major axis or similar and these
    ## ## are counts.
    ## expect_true(linearHypothesis(lm1, diag(2), c(0, 1))[["Pr(>F)"]][2] >
    ##             p.fail)

})
date()





date()
test_that("McFL: Mutator increases by given factor with per-gene-mut rates: major axis and chi-sq test", {

    ## Two cases: mutator and no mutator, with variable mutation rates.
    ## rates such that rates of no mutator = rates of mutator * mutator.

    ## Why not compare mutlitplication factor keeping mutation rates
    ## constant? Because specially with mutator and large diffs in mut
    ## rates, with oncoSimulSample you undersample variation with
    ## wholePop, etc.

    ## Setings similar to oss11 in per-gene-mutation-rates but with the mutator
    
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n MCFL: AEu8: the seed is", pseed, "\n")
    pops <- 2000
    ft <- 5e-3
    lni <- 7
    no <- 5e5
    ni <- c(0, 0, 0, rep(0, lni))
    ## scramble around names
    names(ni) <- c("hereisoneagene",
                   "oreoisasabgene",
                   "nnhsisthecgene",
                   replicate(lni,
                             paste(sample(letters, 12), collapse = "")))
    ni <- ni[order(names(ni))]
    fe <- allFitnessEffects(noIntGenes = ni)
    mutator1 <- rep(1, lni + 3)
    pg1 <- seq(from = 1e-9, to = 1e-6, length.out = lni + 3) ## max should not be
                                                  ## huge here as mutator
                                                  ## is 34. Can get beyond
                                                  ## 1
    names(mutator1) <- sample(names(ni))
    names(pg1) <- sample(names(ni))
    mutator1["oreoisasabgene"] <- 100
    m1 <- allMutatorEffects(noIntGenes = mutator1)
    ## pg1["hereisoneagene"] <- 1e-4 ## if this gets huge, then you are
    ##                               ## undersampling and the chi-square will
    ##                               ## fail. But then, we probably are
    ##                               ## running into numerical issues: 3
    ##                               ## orders of magnitude differences.
    m1.pg1.b <- oncoSimulSample(pops,
                           fe,
                           mu = pg1,
                           muEF = m1,
                           finalTime = ft,
                           mutationPropGrowth = FALSE,
                           initSize = no,
                           initMutant ="oreoisasabgene",
                           model = "McFL",
                           sampleEvery = 0.01, thresholdWhole = 1e-20,
                           detectionSize = 1e9,
                           detectionDrivers = 9999,
                           onlyCancer = FALSE, seed = NULL)
    ## m1.pg1.b$popSummary[, c(1:3, 8:9)]
    summary(m1.pg1.b$popSummary[, "NumClones"])
    ## Recall that init-mutant tests check always present of initMutant
    ## against a thresholWhole of 1. Here it is slightly different.
    expect_true(smSampl("oreoisasabgene", m1.pg1.b) == pops)
    ## catch a pattern that would make the previous trivially true
    expect_false(sum(m1.pg1.b$popSample) == pops * (lni + 3))
    ## next two, to compare with oss1a
    sort(enom("oreoisasabgene", pg1, no, pops))
    sort(snomSampl("oreoisasabgene", m1.pg1.b))
    ## Compare with the expected for this scenario
    p.fail <- 1e-3
    expect_true(chisq.test(snomSampl("oreoisasabgene", m1.pg1.b),
                           p = pnom("oreoisasabgene", pg1, no, pops))$p.value > p.fail)

    pg2 <- seq(from = 1e-7, to = 1e-4, length.out = lni + 3)
    names(pg2) <- names(pg1)
    m1.pg2.b <- oncoSimulSample(pops,
                           fe,
                           mu = pg2,
                           finalTime = ft,
                           mutationPropGrowth = FALSE,
                           initSize = no,
                           initMutant ="oreoisasabgene",
                           model = "McFL",
                           sampleEvery = 0.01, thresholdWhole = 1e-20,
                           detectionSize = 1e9,
                           detectionDrivers = 9999,
                           onlyCancer = FALSE, seed = NULL)
    ## m1.pg2.b$popSummary[, c(1:3, 8:9)]
    summary(m1.pg2.b$popSummary[, "NumClones"])
    ## Recall that init-mutant tests check always present of initMutant
    ## against a thresholWhole of 1. Here it is slightly different.
    expect_true(smSampl("oreoisasabgene", m1.pg2.b) == pops)
    ## catch a pattern that would make the previous trivially true
    expect_false(sum(m1.pg2.b$popSample) == pops * (lni + 3))
    ## next two, to compare with oss1a
    sort(enom("oreoisasabgene", pg2, no, pops))
    sort(snomSampl("oreoisasabgene", m1.pg2.b))
    p.fail <- 1e-3
    expect_true(chisq.test(snomSampl("oreoisasabgene", m1.pg2.b),
                           p = pnom("oreoisasabgene", pg2, no, pops))$p.value > p.fail)


    ## Compare mutator with no mutator
    expect_true(chisq.test(snomSampl("oreoisasabgene", m1.pg1.b),
                           snomSampl("oreoisasabgene", m1.pg2.b))$p.value > p.fail)
    y <- sqrt(snomSampl("oreoisasabgene", m1.pg1.b))
    x <- sqrt(snomSampl("oreoisasabgene", m1.pg2.b))
    mma <- smatr::ma(y ~ x, slope.test = 1, elev.test = 0) ## From smatr package, for major axis
    ## intercept not different from 0
    expect_true(mma$elevtest[[1]]$p > p.fail)
    expect_true(mma$slopetest[[1]]$p > p.fail)
    
    ## We could use a lm and do a simultaneous test on both slope and
    ## intercept as. But this is really asking for major axis regression

    ## lm1 <- lm(snomSampl("oreoisasabgene", m1.pg1.b) ~
    ##               snomSampl("oreoisasabgene", m1.pg2.b))
    ## ## test intercept is 0, slope is 1. Not technically fully correct, as
    ## ## X variable has noise. We should do major axis or similar and these
    ## ## are counts.
    ## expect_true(linearHypothesis(lm1, diag(2), c(0, 1))[["Pr(>F)"]][2] >
    ##             p.fail)

})
date()




## Slow (~ 3 seconds) but tests modules of mutator nicely.
date() ## Beware: this uses a lot of RAM without the gc()
test_that("Mutator modules differences", {
    
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n mmd1: the seed is", pseed, "\n")
    reps <- 10
    no <- 5e3
    ft <- 100
    mu <- 1e-5
    lni <- 50
    m1 <- 1
    m2 <- 25
    m3 <- 50
    ni <- rep(0, lni)
    gn <- paste0("a", 1:lni)
    names(ni) <- gn
    gn <- paste(gn, collapse = ", ")
    mut1 <- allMutatorEffects(epistasis = c("A" = m1),
                              geneToModule = c("A" = gn))
    mut2 <- allMutatorEffects(epistasis = c("A" = m2),
                              geneToModule = c("A" = gn))
    mut3 <- allMutatorEffects(epistasis = c("A" = m3),
                              geneToModule = c("A" = gn))
    f1 <- allFitnessEffects(noIntGenes = ni)
    b1 <- oncoSimulSample(reps,
                       f1,
                       mu = mu,
                       muEF = mut1,
                       onlyCancer = FALSE,
                       initSize = no,
                       finalTime = ft,
                        sampleEvery = 0.01, thresholdWhole = 1e-20,
                           detectionSize = 1e9,
                           detectionDrivers = 9999,
                       seed = NULL
                       )
    gc()
    b2 <- oncoSimulSample(reps,
                       f1,
                       mu = mu,
                       muEF = mut2,
                       onlyCancer = FALSE,
                       initSize = no,
                       finalTime = ft,
                        sampleEvery = 0.01, thresholdWhole = 1e-20,
                           detectionSize = 1e9,
                           detectionDrivers = 9999,
                       seed = NULL
                       )
    gc()
    b3 <- oncoSimulSample(reps,
                       f1,
                       mu = mu,
                       muEF = mut3,
                       onlyCancer = FALSE,
                       initSize = no,
                       finalTime = ft,
                        sampleEvery = 0.01, thresholdWhole = 1e-20,
                           detectionSize = 1e9,
                           detectionDrivers = 9999,
                       seed = NULL
                       )
    gc()
    b3$popSummary[, c(1:3, 8:9)]
    b2$popSummary[, c(1:3, 8:9)]
    b1$popSummary[, c(1:3, 8:9)]
    ## mean(rowSums(b3$popSample))
    ## mean(rowSums(b2$popSample))
    ## mean(rowSums(b1$popSample))
    ## This is, of course, affected by sampling only at end: we do not see
    ## the many intermediate events.
    p.fail <- 0.05
    expect_true( t.test( b3$popSummary[, "NumClones"], 
                 b2$popSummary[, "NumClones"], alternative = "greater")$p.value < p.fail)
    expect_true( t.test( b2$popSummary[, "NumClones"], 
                 b1$popSummary[, "NumClones"], alternative = "greater")$p.value < p.fail)
    expect_true( t.test( rowSums(b3$popSample) ,
                 rowSums(b2$popSample), alternative = "greater")$p.value < p.fail)
    expect_true( t.test( rowSums(b2$popSample) ,
                 rowSums(b1$popSample), alternative = "greater")$p.value < p.fail)
    
})
date()




date() ## Beware: this uses a lot of RAM without the gc()
test_that("McFL: Mutator modules differences", {
    
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n MCFLmmd1: the seed is", pseed, "\n")
    reps <- 10
    no <- 5e3
    ft <- 100
    mu <- 1e-5
    lni <- 50
    m1 <- 1
    m2 <- 25
    m3 <- 50
    ni <- rep(0, lni)
    gn <- paste0("a", 1:lni)
    names(ni) <- gn
    gn <- paste(gn, collapse = ", ")
    mut1 <- allMutatorEffects(epistasis = c("A" = m1),
                              geneToModule = c("A" = gn))
    mut2 <- allMutatorEffects(epistasis = c("A" = m2),
                              geneToModule = c("A" = gn))
    mut3 <- allMutatorEffects(epistasis = c("A" = m3),
                              geneToModule = c("A" = gn))
    f1 <- allFitnessEffects(noIntGenes = ni)
    b1 <- oncoSimulSample(reps,
                       f1,
                       mu = mu,
                       muEF = mut1,
                       onlyCancer = FALSE,
                       initSize = no,
                       finalTime = ft,
                        sampleEvery = 0.01, thresholdWhole = 1e-20,
                           detectionSize = 1e9,
                           detectionDrivers = 9999,
                        model = "McFL",
                       seed = NULL
                       )
    gc()
    b2 <- oncoSimulSample(reps,
                       f1,
                       mu = mu,
                       muEF = mut2,
                       onlyCancer = FALSE,
                       initSize = no,
                       finalTime = ft,
                        sampleEvery = 0.01, thresholdWhole = 1e-20,
                           detectionSize = 1e9,
                           detectionDrivers = 9999,
                        model = "McFL",
                       seed = NULL
                       )
    gc()
    b3 <- oncoSimulSample(reps,
                       f1,
                       mu = mu,
                       muEF = mut3,
                       onlyCancer = FALSE,
                       initSize = no,
                       finalTime = ft,
                        sampleEvery = 0.01, thresholdWhole = 1e-20,
                           detectionSize = 1e9,
                           detectionDrivers = 9999,
                        model = "McFL",
                       seed = NULL
                       )
    gc()
    b3$popSummary[, c(1:3, 8:9)]
    b2$popSummary[, c(1:3, 8:9)]
    b1$popSummary[, c(1:3, 8:9)]
    ## mean(rowSums(b3$popSample))
    ## mean(rowSums(b2$popSample))
    ## mean(rowSums(b1$popSample))
    ## This is, of course, affected by sampling only at end: we do not see
    ## the many intermediate events.
    p.fail <- 0.05
    expect_true( t.test( b3$popSummary[, "NumClones"], 
                 b2$popSummary[, "NumClones"], alternative = "greater")$p.value < p.fail)
    expect_true( t.test( b2$popSummary[, "NumClones"], 
                 b1$popSummary[, "NumClones"], alternative = "greater")$p.value < p.fail)
    expect_true( t.test( rowSums(b3$popSample) ,
                 rowSums(b2$popSample), alternative = "greater")$p.value < p.fail)
    expect_true( t.test( rowSums(b2$popSample) ,
                 rowSums(b1$popSample), alternative = "greater")$p.value < p.fail)
    
})
date()

    

date() 
test_that("Mutator, several modules differences", {
    pseed <- sample(99999999, 1)
    set.seed(pseed)
    cat("\n mmd1: the seed is", pseed, "\n")
    reps <- 100
    no <- 5e3
    ft <- 50 ## you need it large enough to get enough hits
    mu <- 1e-5
    ln <- 50 
    m1 <- 5 ## if this is too large, easy to get it to blow.
    ni <- rep(0, 3 * ln)
    gna <- paste0("a", 1:ln)
    gnb <- paste0("b", 1:ln)
    gnc <- paste0("c", 1:ln)
    names(ni) <- c(gna, gnb, gnc)
    gn1 <- paste(c(gna, gnb, gnc), collapse = ", ")
    gna <- paste(gna, collapse = ", ")
    gnb <- paste(gnb, collapse = ", ")
    gnc <- paste(gnc, collapse = ", ")
    mut1 <- allMutatorEffects(epistasis = c("A" = m1),
                              geneToModule = c("A" = gn1))
    mut2 <- allMutatorEffects(epistasis = c("A" = m1,
                                            "B" = m1,
                                            "C" = m1),
                              geneToModule = c("A" = gna,
                                               "B" = gnb,
                                               "C" = gnc))
    f1 <- allFitnessEffects(noIntGenes = ni)
    b1 <- oncoSimulSample(reps,
                       f1,
                       mu = mu,
                       muEF = mut1,
                       onlyCancer = FALSE,
                       initSize = no,
                       finalTime = ft,
                        sampleEvery = 0.01, thresholdWhole = 1e-20,
                           detectionSize = 1e9,
                           detectionDrivers = 9999,
                       seed = NULL
                       )
    gc()
    b2 <- oncoSimulSample(reps,
                       f1,
                       mu = mu,
                       muEF = mut2,
                       onlyCancer = FALSE,
                       initSize = no,
                       finalTime = ft,
                        sampleEvery = 0.01, thresholdWhole = 1e-20,
                           detectionSize = 1e9,
                           detectionDrivers = 9999,
                       seed = NULL
                       )
    gc()
    ## b2$popSummary[, c(1:3, 8:9)]
    ## b1$popSummary[, c(1:3, 8:9)]
    ## mean(rowSums(b2$popSample))
    ## mean(rowSums(b1$popSample))
    ## This is, of course, affected by sampling only at end: we do not see
    ## the many intermediate events.
    ## Variances for NumClones are hugely unequal, even after log transform.;
    ## might want Wilcoxon? Similar for rowSums of popSample
    p.fail <- 0.05
    expect_true( wilcox.test( b2$popSummary[, "NumClones"], 
                 b1$popSummary[, "NumClones"], alternative = "greater")$p.value < p.fail)
    expect_true( wilcox.test( rowSums(b2$popSample) ,
                 rowSums(b1$popSample), alternative = "greater")$p.value < p.fail)
})
date()


## Remember that numClones is underestimated, possibly severly, by
## oncoSimulSample compared to oncoSimulPop, since we only look at the
## clones that exist at the end.

date() 
test_that("Mutator, several modules differences, McFL", {
    
    pseed <- sample(99999999, 1)
    pseed <- 91339980
    set.seed(pseed)
    cat("\n mmd1: the seed is", pseed, "\n")
    reps <- 60
    no <- 5e3
    ft <- 50 ## you need it large enough to get enough hits
    mu <- 1e-5
    ln <- 50 
    m1 <- 5 ## if this is too large, easy to get it to blow.
    ni <- rep(0, 3 * ln)
    gna <- paste0("a", 1:ln)
    gnb <- paste0("b", 1:ln)
    gnc <- paste0("c", 1:ln)
    names(ni) <- c(gna, gnb, gnc)
    gn1 <- paste(c(gna, gnb, gnc), collapse = ", ")
    gna <- paste(gna, collapse = ", ")
    gnb <- paste(gnb, collapse = ", ")
    gnc <- paste(gnc, collapse = ", ")
    mut1 <- allMutatorEffects(epistasis = c("A" = m1),
                              geneToModule = c("A" = gn1))
    mut2 <- allMutatorEffects(epistasis = c("A" = m1,
                                            "B" = m1,
                                            "C" = m1),
                              geneToModule = c("A" = gna,
                                               "B" = gnb,
                                               "C" = gnc))
    f1 <- allFitnessEffects(noIntGenes = ni)
    b1 <- oncoSimulSample(reps,
                       f1,
                       mu = mu,
                       muEF = mut1,
                       onlyCancer = FALSE,
                       initSize = no,
                       finalTime = ft,
                        sampleEvery = 0.01, thresholdWhole = 1e-20,
                           detectionSize = 1e9,
                           detectionDrivers = 9999,
                       seed = NULL, model = "McFL"
                       )
    gc()
    b2 <- oncoSimulSample(reps,
                       f1,
                       mu = mu,
                       muEF = mut2,
                       onlyCancer = FALSE,
                       initSize = no,
                       finalTime = ft,
                        sampleEvery = 0.01, thresholdWhole = 1e-20,
                           detectionSize = 1e9,
                           detectionDrivers = 9999,
                       seed = NULL, model = "McFL"
                       )
    gc()
    b2$popSummary[, c(1:3, 8:9)]
    b1$popSummary[, c(1:3, 8:9)]
    mean(rowSums(b2$popSample))
    mean(rowSums(b1$popSample))
    ## This is, of course, affected by sampling only at end: we do not see
    ## the many intermediate events.
    p.fail <- 0.05
    expect_true( t.test( b2$popSummary[, "NumClones"], 
                 b1$popSummary[, "NumClones"], alternative = "greater")$p.value < p.fail)
    expect_true( t.test( rowSums(b2$popSample) ,
                 rowSums(b1$popSample), alternative = "greater")$p.value < p.fail)

})
date()






## ## FIXME: fix this with detectionSize and note that rowSums cannot work as
## ## done now, as all have all muts.
## test_that("Relative ordering of number of clones with mutator effects", {
    
##     ## Can occasionally blow up with pE.f: pE not finite.
##     pseed <-sample(9999999, 1)
##     set.seed(pseed)
##     cat("\n x1: the seed is", pseed, "\n")
##     pops <- 20
##     fe <- allFitnessEffects(noIntGenes = c("a" = 0.12,
##                                            "b" = 0.14,
##                                            "c" = 0.16,
##                                            "d" = 0.11))
##     fm6 <- allMutatorEffects(noIntGenes = c("a" = 5,
##                                             "b" = 10,
##                                             "c" = 12,
##                                             "d" = 14))
##     nc1 <- oncoSimulSample(pops, fe, muEF = fm6, finalTime =250,
##                         mutationPropGrowth = FALSE,
##                         sampleEvery = 0.01, thresholdWhole = 1e-20,
##                         initSize = 1e6,
##                         onlyCancer = FALSE)
##     fm8 <- allMutatorEffects(noIntGenes = c("a" = 1,
##                                             "b" = 1,
##                                             "c" = 1,
##                                             "d" = 1))
##     nc2 <- oncoSimulSample(pops, fe, muEF = fm8, finalTime =250,
##                         mutationPropGrowth = FALSE,
##                         sampleEvery = 0.01, thresholdWhole = 1e-20,
##                         initSize = 1e6,
##                         onlyCancer = FALSE)
##     fm7 <- allMutatorEffects(noIntGenes = c("a" = 1e-6,
##                                             "b" = 1e-6,
##                                             "c" = 1e-6,
##                                             "d" = 1e-6))
##     nc3 <- oncoSimulSample(pops, fe, muEF = fm7, finalTime =250,
##                         mutationPropGrowth = FALSE,
##                         sampleEvery = 0.01, thresholdWhole = 1e-20,
##                         initSize = 1e6,
##                         onlyCancer = FALSE)
##     expect_true(median(nc1$popSummary[, "NumClones"]) > median(nc2$popSummary[, "NumClones"]))
##     expect_true(median(nc2$popSummary[, "NumClones"]) > median(nc3$popSummary[, "NumClones"]))
##     expect_true(mean(rowSums(nc1$popSample)) > mean(rowSums(nc2$popSample)))
##     expect_true(mean(rowSums(nc2$popSample)) > mean(rowSums(nc3$popSample)))
##     nc1$popSummary[, c(1:3, 8:9)]
##     nc2$popSummary[, c(1:3, 8:9)]
##     nc3$popSummary[, c(1:3, 8:9)]
    
## })
## date()
## ## FIXME: same as above
## test_that("Relative ordering of number of clones with init mutant of mutator effects", {

##     ## Here we stop on finalTime, not popSize
##     ## Can occasionally blow up with pE.f: pE not finite.
##     pseed <-sample(9999999, 1)
##     set.seed(pseed)
##     cat("\n x2bc: the seed is", pseed, "\n")
##     pops <- 10
##     ni <- rep(0.01, 50)
##     names(ni) <- c("a", "b", "c", "d", paste0("n", 1:46))
##     fe <- allFitnessEffects(noIntGenes = ni)
##     fm6 <- allMutatorEffects(noIntGenes = c("a" = .05,
##                                             "b" = 1,
##                                             "c" = 10,
##                                             "d" = 50))
##     nca <- oncoSimulSample(pops, fe, muEF = fm6, finalTime =50,
##                         mutationPropGrowth = FALSE,
##                         initSize = 1e4,
##                         initMutant = "a",
##                         sampleEvery = 0.01, thresholdWhole = 1e-20,
##                         onlyCancer = FALSE)
##     ncb <- oncoSimulSample(pops, fe, muEF = fm6, finalTime =50,
##                         mutationPropGrowth = FALSE,
##                         initSize = 1e4,
##                         initMutant = "b",
##                         sampleEvery = 0.01, thresholdWhole = 1e-20,
##                         onlyCancer = FALSE)
##     ncc <- oncoSimulSample(pops, fe, muEF = fm6, finalTime =50,
##                         mutationPropGrowth = FALSE,
##                         initSize = 1e4,
##                         initMutant = "c",
##                         sampleEvery = 0.01, thresholdWhole = 1e-20,
##                         onlyCancer = FALSE)
##     ncd <- oncoSimulSample(pops, fe, muEF = fm6, finalTime =50,
##                         mutationPropGrowth = FALSE,
##                         initSize = 1e4,
##                         initMutant = "d",
##                         sampleEvery = 0.01, thresholdWhole = 1e-20,
##                         onlyCancer = FALSE)
##     expect_true( median(nca$popSummary[, "NumClones"]) <
##                  median(ncb$popSummary[, "NumClones"]))
##     expect_true(median(ncb$popSummary[, "NumClones"]) <
##                 median(ncc$popSummary[, "NumClones"]) )
##     expect_true( median(ncc$popSummary[, "NumClones"]) <
##                  median(ncd$popSummary[, "NumClones"]) )
##     expect_true(mean(rowSums(nca$popSample)) < mean(rowSums(ncb$popSample)))
##     expect_true(mean(rowSums(ncb$popSample)) > mean(rowSums(ncc$popSample)))
##     expect_true(mean(rowSums(ncc$popSample)) > mean(rowSums(ncd$popSample)))
##     nca$popSummary[, c(1:3, 8:9)]
##     ncb$popSummary[, c(1:3, 8:9)]
##     ncc$popSummary[, c(1:3, 8:9)]
##     ncd$popSummary[, c(1:3, 8:9)]
       
## })



## test_that("Relative ordering of number of clones with init mutant of mutator effects and s = 0", {
##     ## Can occasionally blow up with pE.f: pE not finite.

##     pseed <-sample(9999999, 1)
##     set.seed(pseed)
##     cat("\n x2cd: the seed is", pseed, "\n")
##     pops <- 40
##     ni <- rep(0, 50)
##     names(ni) <- c("a", "b", "c", "d", paste0("n", 1:46))
##     fe <- allFitnessEffects(noIntGenes = ni)
##     fm6 <- allMutatorEffects(noIntGenes = c("a" = .05,
##                                             "b" = 1,
##                                             "c" = 10,
##                                             "d" = 50))
##     nca <- oncoSimulSample(pops, fe, muEF = fm6, finalTime =50,
##                         mutationPropGrowth = FALSE,
##                         initSize = 1e4,
##                         initMutant = "a",
##                         sampleEvery = 0.01, thresholdWhole = 1e-20,
##                         onlyCancer = FALSE)
##     ncb <- oncoSimulSample(pops, fe, muEF = fm6, finalTime =50,
##                         mutationPropGrowth = FALSE,
##                         initSize = 1e4,
##                         initMutant = "b",
##                         sampleEvery = 0.01, thresholdWhole = 1e-20,
##                         onlyCancer = FALSE)
##     ncc <- oncoSimulSample(pops, fe, muEF = fm6, finalTime =50,
##                         mutationPropGrowth = FALSE,
##                         initSize = 1e4,
##                         initMutant = "c",
##                         sampleEvery = 0.01, thresholdWhole = 1e-20,
##                         onlyCancer = FALSE)
##     ncd <- oncoSimulSample(pops, fe, muEF = fm6, finalTime =50,
##                         mutationPropGrowth = FALSE,
##                         initSize = 1e4,
##                         initMutant = "d",
##                         sampleEvery = 0.01, thresholdWhole = 1e-20,
##                         onlyCancer = FALSE)
##     ## I once saw a weird thing
##     expect_true(var(nca$popSummary[, c(1:3, 8:9)]$NumClones) > 1e-4)
##     expect_true(var(ncb$popSummary[, c(1:3, 8:9)]$NumClones) > 1e-4)
##     expect_true(var(ncc$popSummary[, c(1:3, 8:9)]$NumClones) > 1e-4)
##     expect_true(var(ncd$popSummary[, c(1:3, 8:9)]$NumClones) > 1e-4)
##     ## These are the real tests
##     expect_true( median(nca$popSummary[, "NumClones"]) <
##                  median(ncb$popSummary[, "NumClones"]))
##     expect_true(median(ncb$popSummary[, "NumClones"]) <
##                 median(ncc$popSummary[, "NumClones"]) )
##     expect_true( median(ncc$popSummary[, "NumClones"]) <
##                  median(ncd$popSummary[, "NumClones"]) )
##     expect_true(mean(rowSums(nca$popSample)) < mean(rowSums(ncb$popSample)))
##     expect_true(mean(rowSums(ncb$popSample)) < mean(rowSums(ncc$popSample)))
##     expect_true(mean(rowSums(ncc$popSample)) < mean(rowSums(ncd$popSample)))

## })

## test_that("MCFL Relative ordering of number of clones with mutator effects", {
    
##     ## Can occasionally blow up with pE.f: pE not finite.
##     pseed <-sample(9999999, 1)
##     set.seed(pseed)
##     cat("\n mcx1: the seed is", pseed, "\n")
##     pops <- 30
##     mu <- 1e-7
##     ft <- 200
##     fe <- allFitnessEffects(noIntGenes = c("a" = 0.12,
##                                            "b" = 0.14,
##                                            "c" = 0.16,
##                                            "d" = 0.11))
##     fm6 <- allMutatorEffects(noIntGenes = c("a" = 30,
##                                             "b" = 30,
##                                             "c" = 30,
##                                             "d" = 30))
##     nc1 <- oncoSimulSample(pops, mu = mu,
##                            fe, muEF = fm6, finalTime = ft,
##                         mutationPropGrowth = FALSE,
##                         sampleEvery = 0.01, thresholdWhole = 1e-20,
##                         initSize = 1e6, model = "McFL",
##                         onlyCancer = FALSE)
##     fm8 <- allMutatorEffects(noIntGenes = c("a" = 1,
##                                             "b" = 1,
##                                             "c" = 1,
##                                             "d" = 1))
##     nc2 <- oncoSimulSample(pops, mu = mu,
##                            fe, muEF = fm8, finalTime = ft,
##                         mutationPropGrowth = FALSE,
##                         sampleEvery = 0.01, thresholdWhole = 1e-20,
##                         initSize = 1e6, model = "McFL",
##                         onlyCancer = FALSE)
##     fm7 <- allMutatorEffects(noIntGenes = c("a" = 1e-6,
##                                             "b" = 1e-6,
##                                             "c" = 1e-6,
##                                             "d" = 1e-6))
##     nc3 <- oncoSimulSample(pops, mu = mu,
##                            fe, muEF = fm7, finalTime = ft,
##                         mutationPropGrowth = FALSE,
##                         sampleEvery = 0.01, thresholdWhole = 1e-20,
##                         initSize = 1e6, model = "McFL",
##                         onlyCancer = FALSE)
##     expect_true(median(nc1$popSummary[, "NumClones"]) > median(nc2$popSummary[, "NumClones"]))
##     expect_true(median(nc2$popSummary[, "NumClones"]) > median(nc3$popSummary[, "NumClones"]))
##     nc1$popSummary[, c(1:3, 8:9)]
##     nc2$popSummary[, c(1:3, 8:9)]
##     nc3$popSummary[, c(1:3, 8:9)]
##     expect_true(mean(rowSums(nc1$popSample)) > mean(rowSums(nc2$popSample)))
##     expect_true(mean(rowSums(nc2$popSample)) > mean(rowSums(nc3$popSample)))
    
## })
## date()





## test_that("MCFL Relative ordering of number of clones with init mutant of mutator effects", {

##     ## Here we stop on finalTime, not popSize
##     ## Can occasionally blow up with pE.f: pE not finite.
##     pseed <-sample(9999999, 1)
##     set.seed(pseed)
##     cat("\n mcx2bc: the seed is", pseed, "\n")
##     pops <- 10
##     ni <- rep(0.01, 50)
##     names(ni) <- c("a", "b", "c", "d", paste0("n", 1:46))
##     fe <- allFitnessEffects(noIntGenes = ni)
##     fm6 <- allMutatorEffects(noIntGenes = c("a" = .05,
##                                             "b" = 1,
##                                             "c" = 10,
##                                             "d" = 50))
##     nca <- oncoSimulSample(pops, fe, muEF = fm6, finalTime =50,
##                         mutationPropGrowth = FALSE,
##                         initSize = 1e4,
##                         initMutant = "a",
##                         sampleEvery = 0.01, thresholdWhole = 1e-20,
##                         onlyCancer = FALSE, model = "McFL")
##     ncb <- oncoSimulSample(pops, fe, muEF = fm6, finalTime =50,
##                         mutationPropGrowth = FALSE,
##                         initSize = 1e4,
##                         initMutant = "b",
##                         sampleEvery = 0.01, thresholdWhole = 1e-20,
##                         onlyCancer = FALSE, model = "McFL")
##     ncc <- oncoSimulSample(pops, fe, muEF = fm6, finalTime =50,
##                         mutationPropGrowth = FALSE,
##                         initSize = 1e4,
##                         initMutant = "c",
##                         sampleEvery = 0.01, thresholdWhole = 1e-20,
##                         onlyCancer = FALSE, model = "McFL")
##     ncd <- oncoSimulSample(pops, fe, muEF = fm6, finalTime =50,
##                         mutationPropGrowth = FALSE,
##                         initSize = 1e4,
##                         initMutant = "d",
##                         sampleEvery = 0.01, thresholdWhole = 1e-20,
##                         onlyCancer = FALSE, model = "McFL")
##     expect_true( median(nca$popSummary[, "NumClones"]) <
##                  median(ncb$popSummary[, "NumClones"]))
##     expect_true(median(ncb$popSummary[, "NumClones"]) <
##                 median(ncc$popSummary[, "NumClones"]) )
##     expect_true( median(ncc$popSummary[, "NumClones"]) <
##                  median(ncd$popSummary[, "NumClones"]) )
##     nca$popSummary[, c(1:3, 8:9)]
##     ncb$popSummary[, c(1:3, 8:9)]
##     ncc$popSummary[, c(1:3, 8:9)]
##     ncd$popSummary[, c(1:3, 8:9)]
##     expect_true(mean(rowSums(nca$popSample)) < mean(rowSums(ncb$popSample)))
##     expect_true(mean(rowSums(ncb$popSample)) < mean(rowSums(ncc$popSample)))
##     expect_true(mean(rowSums(ncc$popSample)) < mean(rowSums(ncd$popSample)))
       
## })



## test_that("MCFL Relative ordering of number of clones with init mutant of mutator effects and s = 0", {
##     ## Can occasionally blow up with pE.f: pE not finite.

##     pseed <-sample(9999999, 1)
##     set.seed(pseed)
##     cat("\n mcx2cd: the seed is", pseed, "\n")
##     pops <- 20
##     ni <- rep(0, 50)
##     names(ni) <- c("a", "b", "c", "d", paste0("n", 1:46))
##     fe <- allFitnessEffects(noIntGenes = ni)
##     fm6 <- allMutatorEffects(noIntGenes = c("a" = .05,
##                                             "b" = 1,
##                                             "c" = 10,
##                                             "d" = 50))
##     nca <- oncoSimulSample(pops, fe, muEF = fm6, finalTime =50,
##                         mutationPropGrowth = FALSE,
##                         initSize = 1e4,
##                         initMutant = "a",
##                         sampleEvery = 0.01, thresholdWhole = 1e-20,
##                         onlyCancer = FALSE, model = "McFL")
##     ncb <- oncoSimulSample(pops, fe, muEF = fm6, finalTime =50,
##                         mutationPropGrowth = FALSE,
##                         initSize = 1e4,
##                         initMutant = "b",
##                         sampleEvery = 0.01, thresholdWhole = 1e-20,
##                         onlyCancer = FALSE, model = "McFL")
##     ncc <- oncoSimulSample(pops, fe, muEF = fm6, finalTime =50,
##                         mutationPropGrowth = FALSE,
##                         initSize = 1e4,
##                         initMutant = "c",
##                         sampleEvery = 0.01, thresholdWhole = 1e-20,
##                         onlyCancer = FALSE, model = "McFL")
##     ncd <- oncoSimulSample(pops, fe, muEF = fm6, finalTime =50,
##                         mutationPropGrowth = FALSE,
##                         initSize = 1e4,
##                         initMutant = "d",
##                         sampleEvery = 0.01, thresholdWhole = 1e-20,
##                         onlyCancer = FALSE, model = "McFL")
##     ## I once saw a weird thing
##     expect_true(var(nca$popSummary[, c(1:3, 8:9)]$NumClones) > 1e-4)
##     expect_true(var(ncb$popSummary[, c(1:3, 8:9)]$NumClones) > 1e-4)
##     expect_true(var(ncc$popSummary[, c(1:3, 8:9)]$NumClones) > 1e-4)
##     expect_true(var(ncd$popSummary[, c(1:3, 8:9)]$NumClones) > 1e-4)
##     ## These are the real tests
##     expect_true( median(nca$popSummary[, "NumClones"]) <
##                  median(ncb$popSummary[, "NumClones"]))
##     expect_true(median(ncb$popSummary[, "NumClones"]) <
##                 median(ncc$popSummary[, "NumClones"]) )
##     expect_true( median(ncc$popSummary[, "NumClones"]) <
##                  median(ncd$popSummary[, "NumClones"]) )
##     expect_true(mean(rowSums(nca$popSample)) < mean(rowSums(ncb$popSample)))
##     expect_true(mean(rowSums(ncb$popSample)) < mean(rowSums(ncc$popSample)))
##     expect_true(mean(rowSums(ncc$popSample)) < mean(rowSums(ncd$popSample)))
    
## })



## test_that("Relative ordering of number of clones with mut prop growth and init and scrambled names", {

##     ## Can occasionally blow up with pE.f: pE not finite.
##     pseed <- sample(99999999, 1)
##     set.seed(pseed)
##     cat("\n x2ef: the seed is", pseed, "\n")
##     pops <- 10
##     ft <- 1
##     lni <- 200
##     no <- 5e3
##     ni <- c(5, 2, rep(0, lni))
##     ## scramble around names
##     names(ni) <- c("thisistheagene",
##                    "thisisthebgene",
##                    replicate(lni,
##                              paste(sample(letters, 12), collapse = "")))
##     ni <- ni[order(names(ni))]
##     fe <- allFitnessEffects(noIntGenes = ni)
##     fm1 <- allMutatorEffects(noIntGenes = c("thisistheagene" = 5))
##     mpg <- oncoSimulSample(pops, fe, muEF = fm1,
##                         finalTime = ft,
##                         mutationPropGrowth = TRUE,
##                         initSize = no,
##                         sampleEvery = 0.01, thresholdWhole = 1e-20,
##                         initMutant = "thisistheagene",
##                         onlyCancer = FALSE)
##     mnpg <- oncoSimulSample(pops, fe, muEF = fm1,
##                          finalTime = ft,
##                          mutationPropGrowth = FALSE,
##                          initSize = no,
##                          sampleEvery = 0.01, thresholdWhole = 1e-20,
##                          initMutant = "thisistheagene",
##                          onlyCancer = FALSE)
##     pg <- oncoSimulSample(pops, fe, 
##                        finalTime = ft,
##                        mutationPropGrowth = TRUE,
##                        initSize = no,
##                        sampleEvery = 0.01, thresholdWhole = 1e-20,
##                        initMutant = "thisistheagene",
##                        onlyCancer = FALSE)
##     npg <- oncoSimulSample(pops, fe, 
##                         finalTime = ft,
##                         mutationPropGrowth = FALSE,
##                         initSize = no,
##                         sampleEvery = 0.01, thresholdWhole = 1e-20,
##                         initMutant = "thisistheagene",
##                         onlyCancer = FALSE)
##       ## I once saw a weird thing
##     expect_true(var(mpg$popSummary[, c(1:3, 8:9)]$NumClones) > 1e-4)
##     expect_true(var(mnpg$popSummary[, c(1:3, 8:9)]$NumClones) > 1e-4)
##     expect_true(var(pg$popSummary[, c(1:3, 8:9)]$NumClones) > 1e-4)
##     expect_true(var(npg$popSummary[, c(1:3, 8:9)]$NumClones) > 1e-4)
##     ## These are the real tests
##     expect_true( median(mpg$popSummary[, "NumClones"]) >
##                  median(mnpg$popSummary[, "NumClones"]))
##     expect_true(median(mpg$popSummary[, "NumClones"]) >
##                 median(pg$popSummary[, "NumClones"]) )
##     expect_true( median(mnpg$popSummary[, "NumClones"]) >
##                  median(npg$popSummary[, "NumClones"]) )
##     expect_true( median(pg$popSummary[, "NumClones"]) >
##                  median(npg$popSummary[, "NumClones"]) )
##     expect_true(mean(rowSums(mpg$popSample)) > mean(rowSums(mnpg$popSample)))
##     expect_true(mean(rowSums(mpg$popSample)) > mean(rowSums(pg$popSample)))
##     expect_true(mean(rowSums(mnpg$popSample)) > mean(rowSums(npg$popSample)))
##     expect_true(mean(rowSums(pg$popSample)) > mean(rowSums(npg$popSample)))
    
## })


## test_that("McFL: Relative ordering of number of clones with mut prop growth and init and scrambled names", {
##     ## Can occasionally blow up with pE.f: pE not finite.
##     pseed <-sample(9999999, 1)
##     set.seed(pseed)
##     cat("\n x2gh: the seed is", pseed, "\n")
##     pops <- 10
##     ft <- 1
##     lni <- 200
##     no <- 1e3
##     ni <- c(5, 2, rep(0, lni))
##     ## scramble around names
##     names(ni) <- c("thisistheagene",
##                    "thisisthebgene",
##                    replicate(lni,
##                              paste(sample(letters, 12), collapse = "")))
##     ni <- ni[order(names(ni))]
##     fe <- allFitnessEffects(noIntGenes = ni)
##     fm1 <- allMutatorEffects(noIntGenes = c("thisistheagene" = 5))
##     mpg <- oncoSimulSample(pops, fe, muEF = fm1,
##                         finalTime = ft,
##                         mutationPropGrowth = TRUE,
##                         initSize = no, model = "McFL",
##                         initMutant = "thisistheagene",
##                         sampleEvery = 0.01, thresholdWhole = 1e-20,
##                         onlyCancer = FALSE)
##     mnpg <- oncoSimulSample(pops, fe, muEF = fm1,
##                          finalTime = ft,
##                          mutationPropGrowth = FALSE,
##                          initSize = no, model = "McFL",
##                          initMutant = "thisistheagene",
##                          sampleEvery = 0.01, thresholdWhole = 1e-20,
##                          onlyCancer = FALSE)
##     pg <- oncoSimulSample(pops, fe, 
##                        finalTime = ft,
##                        mutationPropGrowth = TRUE,
##                        initSize = no, model = "McFL",
##                        initMutant = "thisistheagene",
##                        sampleEvery = 0.01, thresholdWhole = 1e-20,
##                        onlyCancer = FALSE)
##     npg <- oncoSimulSample(pops, fe, 
##                         finalTime = ft,
##                         mutationPropGrowth = FALSE,
##                         initSize = no, model = "McFL",
##                         initMutant = "thisistheagene",
##                         sampleEvery = 0.01, thresholdWhole = 1e-20,
##                         onlyCancer = FALSE)
##       ## I once saw a weird thing
##     expect_true(var(mpg$popSummary[, c(1:3, 8:9)]$NumClones) > 1e-4)
##     expect_true(var(mnpg$popSummary[, c(1:3, 8:9)]$NumClones) > 1e-4)
##     expect_true(var(pg$popSummary[, c(1:3, 8:9)]$NumClones) > 1e-4)
##     expect_true(var(npg$popSummary[, c(1:3, 8:9)]$NumClones) > 1e-4)
##     ## These are the real tests
##     expect_true( median(mpg$popSummary[, "NumClones"]) >
##                  median(mnpg$popSummary[, "NumClones"]))
##     expect_true(median(mpg$popSummary[, "NumClones"]) >
##                 median(pg$popSummary[, "NumClones"]) )
##     expect_true( median(mnpg$popSummary[, "NumClones"]) >
##                  median(npg$popSummary[, "NumClones"]) )
##     expect_true( median(pg$popSummary[, "NumClones"]) >
##                  median(npg$popSummary[, "NumClones"]) )
##     expect_true(mean(rowSums(mpg$popSample)) > mean(rowSums(mnpg$popSample)))
##     expect_true(mean(rowSums(mpg$popSample)) > mean(rowSums(pg$popSample)))
##     expect_true(mean(rowSums(mnpg$popSample)) > mean(rowSums(npg$popSample)))
##     expect_true(mean(rowSums(pg$popSample)) > mean(rowSums(npg$popSample)))
## })


##### Comparisons against expected freqs, using a chi-square

## If any mu is very large or any lni is very large, it can fail unless
## pops is very large. And having a large mutator effect is like having a
## very large mu. We want to use very small finalTime: It is birth and
## rate that compound processes and of course we have non-independent
## sampling (overdispersion) which can make chisq a bad idea.

## Thus, we use a tiny final time so we are basically getting just
## mutation events. We need to use a large number of pops to try to avoid
## empty cells with low mutation frequencies.

## We will play with the mutator effects. Note also that here mutator is
## specified passing a vector of same size as genome. It would be faster
## to use just the name of mutator gene.


## Do next also with equal freqs.
date()
test_that("Expect freq genotypes, mutator and var mut rates", {
    
    ## We test that mutator does not affect expected frequencies of
    ## mutated genes: they are given by the mutation rate of each gene.
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n u6: the seed is", pseed, "\n")
    pops <- 1000
    ft <- 1e-5 ## small, as we cannot afford to accumulate many mutations
               ## or else, given that we have a wholePopulation sample, we
               ## get the wrong result. Not the case with single cell sampling.
    lni <- 70  ##80 
    no <- 5e5
    ni <- c(0, 0, 0, rep(0, lni))
    ## scramble around names
    names(ni) <- c("hereisoneagene",
                   "oreoisasabgene",
                   "nnhsisthecgene",
                   replicate(lni,
                             paste(sample(letters, 12), collapse = "")))
    ni <- ni[order(names(ni))]
    fe <- allFitnessEffects(noIntGenes = ni)
    ## of course, passing a mutator of 1 makes everything slow.
    mutator1 <- rep(1, lni + 3)
    ## pg1 <- rep(1e-5, lni + 3)
    pg1 <- runif(lni + 3, min = 1e-5, max = 5e-4) ## max should not be
                                                  ## huge here as mutator
                                                  ## is 34. Can get beyond
                                                  ## 1
    names(mutator1) <- sample(names(ni))
    names(pg1) <- sample(names(ni))
    mutator1["oreoisasabgene"] <- 10 ## 34    ## 53
    m1 <- allMutatorEffects(noIntGenes = mutator1)
    ## have something with much larger mutation rate
    pg1["hereisoneagene"] <- 1e-3 ## have something huge
    m1.pg1.b <- oncoSimulSample(pops,
                           fe,
                           mu = pg1,
                           muEF = m1,
                           finalTime = ft,
                           mutationPropGrowth = FALSE,
                           initSize = no,
                           initMutant ="oreoisasabgene",
                           sampleEvery = 0.01, thresholdWhole = 1e-20,
                           detectionSize = 1e9,
                           detectionDrivers = 9999,
                           onlyCancer = FALSE, seed = NULL)
    ## If numclones is much larger than 2, that signals trouble as you are
    ## smoothing differences between frequencies with oncoSimulSample,
    ## whole pop
    summary(m1.pg1.b$popSummary[, "NumClones"])
    ## Recall that init-mutant tests check always present of initMutant
    ## against a thresholWhole of 1. Here it is slightly different.
    expect_true(smSampl("oreoisasabgene", m1.pg1.b) == pops)
    ## catch a pattern that would make the previous trivially true
    expect_false(sum(m1.pg1.b$popSample) == pops * (lni + 3)) 
    pnom("oreoisasabgene", pg1, no, pops)
    snomSampl("oreoisasabgene", m1.pg1.b)
    plot(snomSampl("oreoisasabgene", m1.pg1.b)/sum(snomSampl("oreoisasabgene", m1.pg1.b)) ~ 
         pnom("oreoisasabgene", pg1, no, pops)); abline(a = 0, b = 1)
         ## yes, if very large prob for one, it is slightly underestimated
    p.fail <- 1e-3
    expect_true(chisq.test(snomSampl("oreoisasabgene", m1.pg1.b),
                           p = pnom("oreoisasabgene", pg1, no, pops))$p.value > p.fail)

})
date()



date()
test_that("Expect freq genotypes, mutator and var mut rates", {
    
    ## We test that mutator does not affect expected frequencies of
    ## mutated genes: they are given by the mutation rate of each gene.
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n sameu6: the seed is", pseed, "\n")
    pops <- 1000
    ft <- 1e-5 ## small, as we cannot afford to accumulate many mutations
               ## or else, given that we have a wholePopulation sample, we
               ## get the wrong result. Not the case with single cell sampling.
    lni <- 70  ##80 
    no <- 5e5
    ni <- c(0, 0, 0, rep(0, lni))
    ## scramble around names
    names(ni) <- c("hereisoneagene",
                   "oreoisasabgene",
                   "nnhsisthecgene",
                   replicate(lni,
                             paste(sample(letters, 12), collapse = "")))
    ni <- ni[order(names(ni))]
    fe <- allFitnessEffects(noIntGenes = ni)
    ## of course, passing a mutator of 1 makes everything slow.
    mutator1 <- rep(1, lni + 3)
    ## pg1 <- rep(1e-5, lni + 3)
    pg1 <- runif(lni + 3, min = 5e-4, max = 5e-4) ## max should not be
                                                  ## huge here as mutator
                                                  ## is 34. Can get beyond
                                                  ## 1
    names(mutator1) <- sample(names(ni))
    names(pg1) <- sample(names(ni))
    mutator1["oreoisasabgene"] <- 10 ## 34    ## 53
    m1 <- allMutatorEffects(noIntGenes = mutator1)
    ## have something with much larger mutation rate
    pg1["hereisoneagene"] <- 1e-3 ## have something huge
    m1.pg1.b <- oncoSimulSample(pops,
                           fe,
                           mu = pg1,
                           muEF = m1,
                           finalTime = ft,
                           mutationPropGrowth = FALSE,
                           initSize = no,
                           initMutant ="oreoisasabgene",
                           sampleEvery = 0.01, thresholdWhole = 1e-20,
                           detectionSize = 1e9,
                           detectionDrivers = 9999,
                           onlyCancer = FALSE, seed = NULL)
    ## If numclones is much larger than 2, that signals trouble as you are
    ## smoothing differences between frequencies with oncoSimulSample,
    ## whole pop
    summary(m1.pg1.b$popSummary[, "NumClones"])
    ## Recall that init-mutant tests check always present of initMutant
    ## against a thresholWhole of 1. Here it is slightly different.
    expect_true(smSampl("oreoisasabgene", m1.pg1.b) == pops)
    ## catch a pattern that would make the previous trivially true
    expect_false(sum(m1.pg1.b$popSample) == pops * (lni + 3)) 
    pnom("oreoisasabgene", pg1, no, pops)
    snomSampl("oreoisasabgene", m1.pg1.b)
    plot(snomSampl("oreoisasabgene", m1.pg1.b)/sum(snomSampl("oreoisasabgene", m1.pg1.b)) ~ 
         pnom("oreoisasabgene", pg1, no, pops)); abline(a = 0, b = 1)
         ## yes, if very large prob for one, it is slightly underestimated
    p.fail <- 1e-3
    expect_true(chisq.test(snomSampl("oreoisasabgene", m1.pg1.b),
                           p = pnom("oreoisasabgene", pg1, no, pops))$p.value > p.fail)

})
date()


date()
test_that("Expect freq genotypes, mutator and var mut rates", {
    
    ## Similar to above, but mutator has a single element, not the whole
    ## vector.
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n u7: the seed is", pseed, "\n")
    pops <- 2000
    ft <- 1e-7
    lni <- 80 
    no <- 5e7
    ni <- c(0, 0, 0, rep(0, lni))
    ## scramble around names
    names(ni) <- c("hereisoneagene",
                   "oreoisasabgene",
                   "nnhsisthecgene",
                   replicate(lni,
                             paste(sample(letters, 12), collapse = "")))
    ni <- ni[order(names(ni))]
    fe <- allFitnessEffects(noIntGenes = ni)
    pg1 <- runif(lni + 3, min = 1e-7, max = 1e-4) ## max should not be
                                                  ## huge here as mutator
                                                  ## is 34. Can get beyond
                                                  ## 1
    names(pg1) <- sample(names(ni))
    mutator1 <- c("oreoisasabgene" = 50) ## a single entry
    m1 <- allMutatorEffects(noIntGenes = mutator1)
    pg1["hereisoneagene"] <- 1e-3 ## to compare with a laarge one
    m1.pg1.b <- oncoSimulSample(pops,
                           fe,
                           mu = pg1,
                           muEF = m1,
                           finalTime = ft,
                           mutationPropGrowth = FALSE,
                           initSize = no,
                           initMutant ="oreoisasabgene",
                           sampleEvery = 0.01, thresholdWhole = 1e-20,
                           detectionSize = 1e9,
                           detectionDrivers = 9999,
                           onlyCancer = FALSE, seed = NULL)
    ## If numclones is much larger than 2, that signals trouble as you are
    ## smoothing differences between frequencies with oncoSimulSample,
    ## whole pop
    m1.pg1.b$popSummary[, c(1:3, 8:9)]
    summary(m1.pg1.b$popSummary[, "NumClones"])
    ## Recall that init-mutant tests check always present of initMutant
    ## against a thresholWhole of 1. Here it is slightly different.
    expect_true(smSampl("oreoisasabgene", m1.pg1.b) == pops)
    ## catch a pattern that would make the previous trivially true
    expect_false(sum(m1.pg1.b$popSample) == pops * (lni + 3)) 
    pnom("oreoisasabgene", pg1, no, pops)
    snomSampl("oreoisasabgene", m1.pg1.b)
    plot(snomSampl("oreoisasabgene", m1.pg1.b)/sum(snomSampl("oreoisasabgene", m1.pg1.b)) ~ 
         pnom("oreoisasabgene", pg1, no, pops)); abline(a = 0, b = 1)
         ## yes, if very large prob for one, it is slightly underestimated
    p.fail <- 1e-3
    expect_true(chisq.test(snomSampl("oreoisasabgene", m1.pg1.b),
                           p = pnom("oreoisasabgene", pg1, no, pops))$p.value > p.fail)


})
date()


## ## FIXME: check this test again
## date()
## test_that("Expect freq genotypes, mutator and var mut rates", {
    
##     ## increase mutator, decrease max mu

##     ## similar to oss11 in per-gene-mutation-rates but with the mutator
    
##     pseed <- sample(9999999, 1)
##     set.seed(pseed)
##     cat("\n u8: the seed is", pseed, "\n")
##     pops <- 5000
##     ft <- 5e-3
##     lni <- 7
##     no <- 5e5
##     ni <- c(0, 0, 0, rep(0, lni))
##     ## scramble around names
##     names(ni) <- c("hereisoneagene",
##                    "oreoisasabgene",
##                    "nnhsisthecgene",
##                    replicate(lni,
##                              paste(sample(letters, 12), collapse = "")))
##     ni <- ni[order(names(ni))]
##     fe <- allFitnessEffects(noIntGenes = ni)
##     mutator1 <- rep(1, lni + 3)
##     pg1 <- seq(from = 1e-9, to = 1e-6, length.out = lni + 3) ## max should not be
##                                                   ## huge here as mutator
##                                                   ## is 34. Can get beyond
##                                                   ## 1
##     names(mutator1) <- sample(names(ni))
##     names(pg1) <- sample(names(ni))
##     mutator1["oreoisasabgene"] <- 100
##     m1 <- allMutatorEffects(noIntGenes = mutator1)
##     ## pg1["hereisoneagene"] <- 1e-4 ## if this gets huge, then you are
##     ##                               ## undersampling and the chi-square will
##     ##                               ## fail. But then, we probably are
##     ##                               ## running into numerical issues: 3
##     ##                               ## orders of magnitude differences.
##     m1.pg1.b <- oncoSimulSample(pops,
##                            fe,
##                            mu = pg1,
##                            muEF = m1,
##                            finalTime = ft,
##                            mutationPropGrowth = FALSE,
##                            initSize = no,
##                            initMutant ="oreoisasabgene",
##                            sampleEvery = 0.01, thresholdWhole = 1e-20,
##                            detectionSize = 1e9,
##                            detectionDrivers = 9999,
##                            onlyCancer = FALSE, seed = NULL)
##     m1.pg1.b$popSummary[, c(1:3, 8:9)]
##     ## If numclones is much larger than 2, that signals trouble as you are
##     ## smoothing differences between frequencies with oncoSimulSample,
##     ## whole pop
##     summary(m1.pg1.b$popSummary[, "NumClones"])
##     ## Recall that init-mutant tests check always present of initMutant
##     ## against a thresholWhole of 1. Here it is slightly different.
##     expect_true(smSampl("oreoisasabgene", m1.pg1.b) == pops)
##     ## catch a pattern that would make the previous trivially true
##     expect_false(sum(m1.pg1.b$popSample) == pops * (lni + 3))
##     ## next two, to compare with oss1a
##     sort(enom("oreoisasabgene", pg1, no, pops))
##     sort(snomSampl("oreoisasabgene", m1.pg1.b))
##     ## pnom("oreoisasabgene", pg1, no, pops)
##     ## snomSampl("oreoisasabgene", m1.pg1.b)
##     plot(snomSampl("oreoisasabgene", m1.pg1.b)/sum(snomSampl("oreoisasabgene", m1.pg1.b)) ~ 
##          pnom("oreoisasabgene", pg1, no, pops)); abline(a = 0, b = 1)
##          ## yes, if very large prob for one, it is slightly underestimated
##     p.fail <- 1e-3
##     expect_true(chisq.test(snomSampl("oreoisasabgene", m1.pg1.b),
##                            p = pnom("oreoisasabgene", pg1, no, pops))$p.value > p.fail)

## })
## date()


date()
test_that("McFL, Expect freq genotypes, mutator and var mut rates", {
    
    ## We test that mutator does not affect expected frequencies of
    ## mutated genes: they are given by the mutation rate of each gene.
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n mcfu6: the seed is", pseed, "\n")
    pops <- 2000
    ft <- 1e-7
    lni <- 80 
    no <- 2e7
    ni <- c(0, 0, 0, rep(0, lni))
    ## scramble around names
    names(ni) <- c("hereisoneagene",
                   "oreoisasabgene",
                   "nnhsisthecgene",
                   replicate(lni,
                             paste(sample(letters, 12), collapse = "")))
    ni <- ni[order(names(ni))]
    fe <- allFitnessEffects(noIntGenes = ni)
    ## of course, passing a mutator of 1 makes everything slow.
    mutator1 <- rep(1, lni + 3)
    ## pg1 <- rep(1e-5, lni + 3)
    pg1 <- runif(lni + 3, min = 1e-7, max = 1e-4) ## max should not be
                                                  ## huge here as mutator
                                                  ## is 34. Can get beyond
                                                  ## 1
    names(mutator1) <- sample(names(ni))
    names(pg1) <- sample(names(ni))
    mutator1["oreoisasabgene"] <- 47
    m1 <- allMutatorEffects(noIntGenes = mutator1)
    ## have something with much larger mutation rate
    pg1["hereisoneagene"] <- 1e-3 ## 1e-3
    m1.pg1.b <- oncoSimulSample(pops,
                           fe,
                           mu = pg1,
                           muEF = m1,
                           finalTime = ft,
                           mutationPropGrowth = FALSE,
                           initSize = no,
                           model = "McFL",
                           initMutant ="oreoisasabgene",
                           sampleEvery = 0.01, thresholdWhole = 1e-20,
                           detectionSize = 1e9,
                           detectionDrivers = 9999,
                           onlyCancer = FALSE, seed = NULL)
    summary(m1.pg1.b$popSummary[, "NumClones"])
    ## m1.pg1.b$popSummary[, c(1:3, 8:9)]
    expect_true(smSampl("oreoisasabgene", m1.pg1.b) == pops)
    enom("oreoisasabgene", pg1, no, pops)
    snomSampl("oreoisasabgene", m1.pg1.b)
    plot(snomSampl("oreoisasabgene", m1.pg1.b)/sum(snomSampl("oreoisasabgene", m1.pg1.b)) ~ 
         pnom("oreoisasabgene", pg1, no, pops)); abline(a = 0, b = 1)
    p.fail <- 1e-3
    expect_true(chisq.test(snomSampl("oreoisasabgene", m1.pg1.b),
                           p = pnom("oreoisasabgene", pg1, no, pops))$p.value > p.fail)
    
    
})
date()


date()
test_that("McFL, Expect freq genotypes, mutator and var mut rates", {
    
    ## We test that mutator does not affect expected frequencies of
    ## mutated genes: they are given by the mutation rate of each gene.
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n mcfu7: the seed is", pseed, "\n")
    pops <- 2000
    ft <- 3e-7
    lni <- 80 
    no <- 2e7
    ni <- c(0, 0, 0, rep(0, lni))
    ## scramble around names
    names(ni) <- c("hereisoneagene",
                   "oreoisasabgene",
                   "nnhsisthecgene",
                   replicate(lni,
                             paste(sample(letters, 12), collapse = "")))
    ni <- ni[order(names(ni))]
    fe <- allFitnessEffects(noIntGenes = ni)
    pg1 <- seq(from = 1e-7, to = 1e-4, length.out = lni + 3) ## max should not be
                                                  ## huge here as mutator
                                                  ## is 34. Can get beyond
                                                  ## 1
    names(pg1) <- sample(names(ni))
    mutator1 <- c("oreoisasabgene" = 20) ## a single entry
    m1 <- allMutatorEffects(noIntGenes = mutator1)
    m1 <- allMutatorEffects(noIntGenes = mutator1)
    ## have something with much larger mutation rate
    pg1["hereisoneagene"] <- 1e-3 ## 1e-3
    m1.pg1.b <- oncoSimulSample(pops,
                           fe,
                           mu = pg1,
                           muEF = m1,
                           finalTime = ft,
                           mutationPropGrowth = FALSE,
                           initSize = no,
                           model = "McFL",
                           initMutant ="oreoisasabgene",
                           sampleEvery = 0.01, thresholdWhole = 1e-20,
                           detectionSize = 1e9,
                           detectionDrivers = 9999,
                           onlyCancer = FALSE, seed = NULL)
    summary(m1.pg1.b$popSummary[, "NumClones"])
    ## m1.pg1.b$popSummary[, c(1:3, 8:9)]
    expect_true(smSampl("oreoisasabgene", m1.pg1.b) == pops)
    enom("oreoisasabgene", pg1, no, pops)
    snomSampl("oreoisasabgene", m1.pg1.b)
    plot(snomSampl("oreoisasabgene", m1.pg1.b)/sum(snomSampl("oreoisasabgene", m1.pg1.b)) ~ 
         pnom("oreoisasabgene", pg1, no, pops)); abline(a = 0, b = 1)
    p.fail <- 1e-3
    expect_true(chisq.test(snomSampl("oreoisasabgene", m1.pg1.b),
                           p = pnom("oreoisasabgene", pg1, no, pops))$p.value > p.fail)
    
    
})
date()






test_that("Same mu vector, different mutator; diffs in number muts, tiny t", {

    ## Here, there is no reproduction or death. Just mutation. And no double
    ## mutants either.
    ## We test:
    ##  - mutator increases mutation rates as seen in:
    ##        - number of clones created
    ##        - number of total mutation events
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n nm0: the seed is", pseed, "\n")
    pops <- 20
    ft <- .0001
    lni <- 100
    no <- 1e7
    fi <- rep(0, lni)
    muvector <- rep(5e-6, lni)
    ## scrambling names
    names(fi) <- replicate(lni,
                           paste(sample(letters, 12), collapse = ""))
    names(muvector) <- sample(names(fi))
    ## choose something for mutator
    mutator10 <- mutator100 <- fi[5]
    mutator10[] <- 10
    mutator100[] <- 100
    fe <- allFitnessEffects(noIntGenes = fi)
    m10 <- allMutatorEffects(noIntGenes = mutator10)
    m100 <- allMutatorEffects(noIntGenes = mutator100)
    pop10 <- oncoSimulSample(pops,
                        fe,
                        mu = muvector,
                        muEF = m10,
                        finalTime = ft,
                        mutationPropGrowth = FALSE,
                        initSize = no,
                        initMutant = names(mutator10),
                        sampleEvery = 0.01, thresholdWhole = 1e-20,
                           detectionSize = 1e9,
                           detectionDrivers = 9999,
                        onlyCancer = FALSE)
    pop100 <- oncoSimulSample(pops,
                        fe,
                        mu = muvector,
                        muEF = m100,
                        finalTime = ft,
                        mutationPropGrowth = FALSE,
                        initSize = no,
                        initMutant = names(mutator10),
                        sampleEvery = 0.01, thresholdWhole = 1e-20,
                           detectionSize = 1e9,
                           detectionDrivers = 9999,
                        onlyCancer = FALSE)
    ## number of total mutations do not make sense with oncoSimulSample,
    ## since we cannot estimate them. we approximate wit sum of
    ## mutations. but that is too thick grain.  number of clones is much cleaner
    expect_true(medianNClonesOSS(pop10) < medianNClonesOSS(pop100))
    expect_true(mean(rowSums(pop10$popSample)) <
                mean(rowSums(pop100$popSample)))
    
})


test_that("Same mu vector, different mutator; diffs in number muts, larger t", {
    
    ## reproduction, death, and double and possibly triple mutants. We
    ## decrease init pop size to make this fast.
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n nm1: the seed is", pseed, "\n")
    pops <- 20
    ft <- 1
    lni <- 100
    no <- 1e5
    fi <- rep(0, lni)
    muvector <- rep(5e-6, lni)
    ## scrambling names
    names(fi) <- replicate(lni,
                           paste(sample(letters, 12), collapse = ""))
    names(muvector) <- sample(names(fi))
    ## choose something for mutator
    mutator10 <- mutator100 <- fi[5]
    mutator10[] <- 10
    mutator100[] <- 100
    fe <- allFitnessEffects(noIntGenes = fi)
    m10 <- allMutatorEffects(noIntGenes = mutator10)
    m100 <- allMutatorEffects(noIntGenes = mutator100)
    pop10 <- oncoSimulSample(pops,
                        fe,
                        mu = muvector,
                        muEF = m10,
                        finalTime = ft,
                        mutationPropGrowth = FALSE,
                        initSize = no,
                        initMutant = names(mutator10),
                        sampleEvery = 0.01, thresholdWhole = 1e-20,
                           detectionSize = 1e9,
                           detectionDrivers = 9999,
                        onlyCancer = FALSE)
    pop100 <- oncoSimulSample(pops,
                        fe,
                        mu = muvector,
                        muEF = m100,
                        finalTime = ft,
                        mutationPropGrowth = FALSE,
                        initSize = no,
                        initMutant = names(mutator10),
                        sampleEvery = 0.01, thresholdWhole = 1e-20,
                           detectionSize = 1e9,
                           detectionDrivers = 9999,
                        onlyCancer = FALSE)
        ## number of clones
    expect_true(medianNClonesOSS(pop10) < medianNClonesOSS(pop100))
    expect_true(mean(rowSums(pop10$popSample)) <
                mean(rowSums(pop100$popSample)))

    
})
date()




date()
test_that("McFL: Same mu vector, different mutator; diffs in number muts, tiny t", {

    ## Here, there is no reproduction or death. Just mutation. And no double
    ## mutants either.
    ## We test:
    ##  - mutator increases mutation rates as seen in:
    ##        - number of clones created
    ##        - number of total mutation events
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n nm2: the seed is", pseed, "\n")
    pops <- 20
    ft <- .0001
    lni <- 100
    no <- 1e7
    fi <- rep(0, lni)
    muvector <- rep(5e-6, lni)
    ## scrambling names
    names(fi) <- replicate(lni,
                           paste(sample(letters, 12), collapse = ""))
    names(muvector) <- sample(names(fi))
    ## choose something for mutator
    mutator10 <- mutator100 <- fi[5]
    mutator10[] <- 10
    mutator100[] <- 100
    fe <- allFitnessEffects(noIntGenes = fi)
    m10 <- allMutatorEffects(noIntGenes = mutator10)
    m100 <- allMutatorEffects(noIntGenes = mutator100)
    pop10 <- oncoSimulSample(pops,
                        fe,
                        mu = muvector,
                        muEF = m10,
                        model = "McFL",
                        finalTime = ft,
                        mutationPropGrowth = FALSE,
                        initSize = no,
                        initMutant = names(mutator10),
                        sampleEvery = 0.01, thresholdWhole = 1e-20,
                           detectionSize = 1e9,
                           detectionDrivers = 9999,
                        onlyCancer = FALSE)
    pop100 <- oncoSimulSample(pops,
                        fe,
                        mu = muvector,
                        muEF = m100,
                        model = "McFL",                        
                        finalTime = ft,
                        mutationPropGrowth = FALSE,
                        initSize = no,
                        initMutant = names(mutator10),
                        sampleEvery = 0.01, thresholdWhole = 1e-20,
                           detectionSize = 1e9,
                           detectionDrivers = 9999,
                        onlyCancer = FALSE)
    expect_true(medianNClonesOSS(pop10) < medianNClonesOSS(pop100))
    expect_true(mean(rowSums(pop10$popSample)) <
                mean(rowSums(pop100$popSample)))


})
date()


date()
test_that("McFL: Same mu vector, different mutator; diffs in number muts, larger t", {
    
    ## reproduction, death, and double and possibly triple mutants. We
    ## decrease init pop size to make this fast.
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n nm3: the seed is", pseed, "\n")
    pops <- 20
    ft <- 1
    lni <- 100
    no <- 1e5
    fi <- rep(0, lni)
    muvector <- rep(5e-6, lni)
    ## scrambling names
    names(fi) <- replicate(lni,
                           paste(sample(letters, 12), collapse = ""))
    names(muvector) <- sample(names(fi))
    ## choose something for mutator
    mutator10 <- mutator100 <- fi[5]
    mutator10[] <- 10
    mutator100[] <- 100
    fe <- allFitnessEffects(noIntGenes = fi)
    m10 <- allMutatorEffects(noIntGenes = mutator10)
    m100 <- allMutatorEffects(noIntGenes = mutator100)
    pop10 <- oncoSimulSample(pops,
                        fe,
                        mu = muvector,
                        muEF = m10,
                        model = "McFL",                        
                        finalTime = ft,
                        mutationPropGrowth = FALSE,
                        initSize = no,
                        initMutant = names(mutator10),
                        sampleEvery = 0.01, thresholdWhole = 1e-20,
                           detectionSize = 1e9,
                           detectionDrivers = 9999,
                        onlyCancer = FALSE)
    pop100 <- oncoSimulSample(pops,
                        fe,
                        mu = muvector,
                        muEF = m100,
                        model = "McFL",                        
                        finalTime = ft,
                        mutationPropGrowth = FALSE,
                        initSize = no,
                        initMutant = names(mutator10),
                        sampleEvery = 0.01, thresholdWhole = 1e-20,
                           detectionSize = 1e9,
                           detectionDrivers = 9999,
                        onlyCancer = FALSE)
    expect_true(medianNClonesOSS(pop10) < medianNClonesOSS(pop100))
    expect_true(mean(rowSums(pop10$popSample)) <
                mean(rowSums(pop100$popSample)))

    
})
date()





date()
test_that(" Init with different mutators", {
    
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n z2: the seed is", pseed, "\n")
    pops <- 40
    ft <- .005
    lni <- 50
    no <- 1e7
    ni <- c(0, 0, 0, rep(0, lni))
    ## scramble around names
    names(ni) <- c("hereisoneagene",
                   "oreoisasabgene",
                   "nnhsisthecgene",
                   replicate(lni,
                             paste(sample(letters, 12), collapse = "")))
    ni <- ni[order(names(ni))]
    fe <- allFitnessEffects(noIntGenes = ni)
    mutator1 <- mutator2 <- rep(1, lni + 3)
    pg1 <- rep(5e-6, lni + 3)
    ## scramble names of mutator and per-gene too
    names(mutator1) <- sample(names(ni))
    names(pg1) <- sample(names(ni))
    mutator1["hereisoneagene"] <- 100
    mutator1["oreoisasabgene"] <- 1
    mutator1["nnhsisthecgene"] <- 0.01
    m1 <- allMutatorEffects(noIntGenes = mutator1)
    m1.pg1.a <- oncoSimulSample(pops,
                           fe,
                           mu = pg1,
                           muEF = m1,
                           finalTime = ft,
                           mutationPropGrowth = FALSE,
                           initSize = no,
                           initMutant = "hereisoneagene",
                           detectionSize = 1e9,
                           detectionDrivers = 9999,
                           sampleEvery = 0.01, thresholdWhole = 1e-20,  seed = NULL,
                           onlyCancer = FALSE)
    m1.pg1.b <- oncoSimulSample(pops,
                             fe,
                             mu = pg1,
                             muEF = m1,
                             finalTime = ft,
                             mutationPropGrowth = FALSE,
                             initSize = no,
                             initMutant = "oreoisasabgene",
                           detectionSize = 1e9,
                           detectionDrivers = 9999,
                             sampleEvery = 0.01, thresholdWhole = 1e-20,  seed = NULL,
                             onlyCancer = FALSE)
    m1.pg1.c <- oncoSimulSample(pops,
                           fe,
                           mu = pg1,
                           muEF = m1,
                           finalTime = ft,
                           mutationPropGrowth = FALSE,
                           initSize = no,
                           initMutant = "nnhsisthecgene",
                           detectionSize = 1e9,
                           detectionDrivers = 9999,
                           sampleEvery = 0.01, thresholdWhole = 1e-20,  seed = NULL,
                           onlyCancer = FALSE)
    expect_true(medianNClonesOSS(m1.pg1.a) > medianNClonesOSS(m1.pg1.b))
    expect_true(medianNClonesOSS(m1.pg1.b) > medianNClonesOSS(m1.pg1.c))
    expect_true(mean(rowSums(m1.pg1.a$popSample)) >
                mean(rowSums(m1.pg1.b$popSample)))
    expect_true(mean(rowSums(m1.pg1.b$popSample)) >
                mean(rowSums(m1.pg1.c$popSample)))


})
date()



date()
test_that(" MCFL Init with different mutators", {
    
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n mcz2: the seed is", pseed, "\n")
    pops <- 40
    ft <- .005
    lni <- 50
    no <- 1e7
    ni <- c(0, 0, 0, rep(0, lni))
    ## scramble around names
    names(ni) <- c("hereisoneagene",
                   "oreoisasabgene",
                   "nnhsisthecgene",
                   replicate(lni,
                             paste(sample(letters, 12), collapse = "")))
    ni <- ni[order(names(ni))]
    fe <- allFitnessEffects(noIntGenes = ni)
    mutator1 <- mutator2 <- rep(1, lni + 3)
    pg1 <- rep(5e-6, lni + 3)
    ## scramble names of mutator and per-gene too
    names(mutator1) <- sample(names(ni))
    names(pg1) <- sample(names(ni))
    mutator1["hereisoneagene"] <- 100
    mutator1["oreoisasabgene"] <- 1
    mutator1["nnhsisthecgene"] <- 0.01
    m1 <- allMutatorEffects(noIntGenes = mutator1)
    m1.pg1.a <- oncoSimulSample(pops,
                           fe,
                           mu = pg1,
                           muEF = m1,
                           finalTime = ft,
                           mutationPropGrowth = FALSE,
                           initSize = no,
                           model = "McFL",
                           initMutant = "hereisoneagene",
                           detectionSize = 1e9,
                           detectionDrivers = 9999,
                           sampleEvery = 0.01, thresholdWhole = 1e-20,  seed = NULL,
                           onlyCancer = FALSE)
    m1.pg1.b <- oncoSimulSample(pops,
                             fe,
                             mu = pg1,
                             muEF = m1,
                             finalTime = ft,
                             mutationPropGrowth = FALSE,
                             initSize = no,
                             model = "McFL",                             
                             initMutant = "oreoisasabgene",
                                                        detectionSize = 1e9,
                           detectionDrivers = 9999,
                             sampleEvery = 0.01, thresholdWhole = 1e-20,  seed = NULL,
                             onlyCancer = FALSE)
    m1.pg1.c <- oncoSimulSample(pops,
                           fe,
                           mu = pg1,
                           muEF = m1,
                           finalTime = ft,
                           mutationPropGrowth = FALSE,
                           initSize = no,
                           model = "McFL",                           
                           initMutant = "nnhsisthecgene",
                           detectionSize = 1e9,
                           detectionDrivers = 9999,
                           sampleEvery = 0.01, thresholdWhole = 1e-20,  seed = NULL,
                           onlyCancer = FALSE)
    expect_true(medianNClonesOSS(m1.pg1.a) > medianNClonesOSS(m1.pg1.b))
    expect_true(medianNClonesOSS(m1.pg1.b) > medianNClonesOSS(m1.pg1.c))
    expect_true(mean(rowSums(m1.pg1.a$popSample)) >
                mean(rowSums(m1.pg1.b$popSample)))
    expect_true(mean(rowSums(m1.pg1.b$popSample)) >
                mean(rowSums(m1.pg1.c$popSample)))
  


})
date()


## ## FIXME: move to long, later, and increase reps. from here till end
## ## very slow, because huge number of clones. But tests several phenomena comprehensively.
## ## same with McFL below
## date()
## test_that("per-gene-mut rates and mutator", {

##     pseed <- sample(9999999, 1)
##     set.seed(pseed)
##     cat("\n oss11: the seed is", pseed, "\n")
##     ng <- 10
##     ni <- rep(0, ng)
##     m1 <- runif(ng, min = 1e-7, max = 5e-6)
##     m2 <- runif(ng, min = 1e-5, max = 1e-4)
##     names(ni) <- names(m1) <- names(m2) <- c(replicate(ng,
##                                  paste(sample(letters, 12), collapse = "")))
##     fe1 <- allFitnessEffects(noIntGenes = ni)
##     ft <- 50
##     no <- 5e5 
##     reps <- 40
##     gn <- paste(names(ni), collapse = ", ")
##     mutator1 <- allMutatorEffects(epistasis = c("MU" = 20),
##                                   geneToModule = c("MU" = gn))
##     mutator2 <- allMutatorEffects(epistasis = c("MU" = 40),
##                                   geneToModule = c("MU" = gn))
##     m1.mutator0 <- oncoSimulSample(reps,
##                            fe1,
##                            mu = m1,
##                            onlyCancer = FALSE,
##                            initSize = no,
##                            finalTime = ft,
##                            detectionSize = 1e9,
##                            detectionDrivers = 9999,
##                            sampleEvery = 0.01, thresholdWhole = 1e-20,
##                            seed = NULL
##                            )
##     m1.mutator1 <- oncoSimulSample(reps,
##                            fe1,
##                            mu = m1,
##                            muEF = mutator1,
##                            onlyCancer = FALSE,
##                            initSize = no,
##                            finalTime = ft,
##                            detectionSize = 1e9,
##                            detectionDrivers = 9999,
##                            sampleEvery = 0.01, thresholdWhole = 1e-20,
##                            seed = NULL
##                            )
##     m1.mutator2 <- oncoSimulSample(reps,
##                            fe1,
##                            mu = m1,
##                            muEF = mutator2,
##                            onlyCancer = FALSE,
##                            initSize = no,
##                            finalTime = ft,
##                            detectionSize = 1e9,
##                            detectionDrivers = 9999,
##                            sampleEvery = 0.01, thresholdWhole = 1e-20,
##                            seed = NULL
##                            )
##     runif(1)
##     m2.mutator0 <- oncoSimulSample(reps,
##                            fe1,
##                            mu = m2,
##                            onlyCancer = FALSE,
##                            initSize = no,
##                            finalTime = ft,
##                            detectionSize = 1e9,
##                            detectionDrivers = 9999,
##                            sampleEvery = 0.01, thresholdWhole = 1e-20,
##                            seed = NULL
##                            )
##     m2.mutator1 <- oncoSimulSample(reps,
##                            fe1,
##                            mu = m2,
##                            muEF = mutator1,
##                            onlyCancer = FALSE,
##                            initSize = no,
##                            finalTime = ft,
##                            detectionSize = 1e9,
##                            detectionDrivers = 9999,
##                            sampleEvery = 0.01, thresholdWhole = 1e-20,
##                            seed = NULL
##                            )
##     m2.mutator2 <- oncoSimulSample(reps,
##                            fe1,
##                            mu = m2,
##                            muEF = mutator2,
##                            onlyCancer = FALSE,
##                            initSize = no,
##                            finalTime = ft,
##                            detectionSize = 1e9,
##                            detectionDrivers = 9999,
##                            sampleEvery = 0.01, thresholdWhole = 1e-20,
##                            seed = NULL
##                        )
##     m1.mutator0$popSummary[, c(1:3, 8:9)]
##     m1.mutator1$popSummary[, c(1:3, 8:9)]
##     m1.mutator2$popSummary[, c(1:3, 8:9)]
##     m2.mutator0$popSummary[, c(1:3, 8:9)]
##     m2.mutator1$popSummary[, c(1:3, 8:9)]
##     m2.mutator2$popSummary[, c(1:3, 8:9)]
##     ## Mutator increases if larger mutator and compared to no mutator
##     ## within levels of per-gene mutation rates
##     expect_true( median(m1.mutator2$popSummary[, "NumClones"]) >
##                  median(m1.mutator1$popSummary[, "NumClones"]))
##     expect_true( median(m1.mutator1$popSummary[, "NumClones"]) >
##                  median(m1.mutator0$popSummary[, "NumClones"]))
##     expect_true( median(m2.mutator2$popSummary[, "NumClones"]) >
##                  median(m2.mutator1$popSummary[, "NumClones"]))
##     expect_true( median(m2.mutator1$popSummary[, "NumClones"]) >
##                  median(m2.mutator0$popSummary[, "NumClones"]))
##     ## can fail easily, because all are equal to pops: all genes mutated in at least one indiv.
##     ## expect_true( mean(mutsPerCloneOSS(m1.mutator2)) >
##     ##              mean(mutsPerCloneOSS(m1.mutator1)))
##     ## expect_true( mean(mutsPerCloneOSS(m1.mutator1)) >
##     ##              mean(mutsPerCloneOSS(m1.mutator0)))
##     ## expect_true( mean(mutsPerCloneOSS(m2.mutator2)) >
##     ##              mean(mutsPerCloneOSS(m2.mutator1)))
##     ## expect_true( mean(mutsPerCloneOSS(m2.mutator1)) >
##     ##              mean(mutsPerCloneOSS(m2.mutator0)))
##     ## Increases in mutation rates increase clones, etc, within levels of
##     ## mutator.
##     expect_true( median(m1.mutator0$popSummary[, "NumClones"]) <
##                  median(m2.mutator0$popSummary[, "NumClones"]))
##     expect_true( median(m1.mutator1$popSummary[, "NumClones"]) <
##                  median(m2.mutator1$popSummary[, "NumClones"]))
##     expect_true( median(m1.mutator2$popSummary[, "NumClones"]) <
##                  median(m2.mutator2$popSummary[, "NumClones"]))
##     ## expect_true( mean(mutsPerCloneOSS(m1.mutator0)) <
##     ##              mean(mutsPerCloneOSS(m2.mutator0)))
##     ## expect_true( mean(mutsPerCloneOSS(m1.mutator1)) <
##     ##              mean(mutsPerCloneOSS(m2.mutator1)))
##     ## expect_true( mean(mutsPerCloneOSS(m1.mutator2)) <
##     ##              mean(mutsPerCloneOSS(m2.mutator2)))
    
## })





## date()
## test_that("McFL: per-gene-mut rates and mutator", {

##     pseed <- sample(9999999, 1)
##     set.seed(pseed)
##     cat("\n mcfloss11: the seed is", pseed, "\n")
##     ng <- 10
##     ni <- rep(0, ng)
##     m1 <- rep(5e-6, ng) ## too much variation and hard to pick the diffs.;
##                         ## runif(ng, min = 1e-7, max = 5e-6) And if too
##                         ## tiny, you do not pick them up unless huge ft
##                         ## and then it is way too slow for m2, etc.
##     m2 <- runif(ng, min = 1e-5, max = 1e-4)
##     names(ni) <- names(m1) <- names(m2) <- c(replicate(ng,
##                                  paste(sample(letters, 12), collapse = "")))
##     fe1 <- allFitnessEffects(noIntGenes = ni)
##     ft <- 50
##     no <- 5e5 
##     reps <- 40
##     gn <- paste(names(ni), collapse = ", ")
##     mutator1 <- allMutatorEffects(epistasis = c("MU" = 20),
##                                   geneToModule = c("MU" = gn))
##     mutator2 <- allMutatorEffects(epistasis = c("MU" = 40),
##                                   geneToModule = c("MU" = gn))
##     m1.mutator0 <- oncoSimulSample(reps,
##                            fe1,
##                            mu = m1,
##                            onlyCancer = FALSE,
##                            initSize = no,
##                            finalTime = ft,
##                            detectionSize = 1e9,
##                            detectionDrivers = 9999,
##                            sampleEvery = 0.01, thresholdWhole = 1e-20,
##                            seed = NULL, model = "McFL"
##                            )
##     m1.mutator1 <- oncoSimulSample(reps,
##                            fe1,
##                            mu = m1,
##                            muEF = mutator1,
##                            onlyCancer = FALSE,
##                            initSize = no,
##                            finalTime = ft,
##                            detectionSize = 1e9,
##                            detectionDrivers = 9999,
##                            sampleEvery = 0.01, thresholdWhole = 1e-20,
##                            seed = NULL, model = "McFL"
##                            )
##     m1.mutator2 <- oncoSimulSample(reps,
##                            fe1,
##                            mu = m1,
##                            muEF = mutator2,
##                            onlyCancer = FALSE,
##                            initSize = no,
##                            finalTime = ft,
##                            detectionSize = 1e9,
##                            detectionDrivers = 9999,
##                            sampleEvery = 0.01, thresholdWhole = 1e-20,
##                            seed = NULL, model = "McFL"
##                            )
##     cat("\n starting m2\n")
##     m2.mutator0 <- oncoSimulSample(reps,
##                            fe1,
##                            mu = m2,
##                            onlyCancer = FALSE,
##                            initSize = no,
##                            finalTime = ft,
##                            detectionSize = 1e9,
##                            detectionDrivers = 9999,
##                            sampleEvery = 0.01, thresholdWhole = 1e-20,
##                            seed = NULL, model = "McFL"
##                            )
##     m2.mutator1 <- oncoSimulSample(reps,
##                            fe1,
##                            mu = m2,
##                            muEF = mutator1,
##                            onlyCancer = FALSE,
##                            initSize = no,
##                            finalTime = ft,
##                            detectionSize = 1e9,
##                            detectionDrivers = 9999,
##                            sampleEvery = 0.01, thresholdWhole = 1e-20,
##                            seed = NULL, model = "McFL"
##                            )
##     m2.mutator2 <- oncoSimulSample(reps,
##                            fe1,
##                            mu = m2,
##                            muEF = mutator2,
##                            onlyCancer = FALSE,
##                            initSize = no,
##                            finalTime = ft,
##                            detectionSize = 1e9,
##                            detectionDrivers = 9999,     
##                            sampleEvery = 0.01, thresholdWhole = 1e-20,
##                            seed = NULL, model = "McFL"
##                            )
##     m1.mutator0$popSummary[, c(1:3, 8:9)]
##     m1.mutator1$popSummary[, c(1:3, 8:9)]
##     m1.mutator2$popSummary[, c(1:3, 8:9)]
##     m2.mutator0$popSummary[, c(1:3, 8:9)]
##     m2.mutator1$popSummary[, c(1:3, 8:9)]
##     m2.mutator2$popSummary[, c(1:3, 8:9)]
##     ## Mutator increases if larger mutator and compared to no mutator
##     ## within levels of per-gene mutation rates
##     ## we could use wilcoxon or t tests actually, specially because often diffs
##     ## are not huge.
##     p.fail <- 0.05
##     expect_true( t.test( m1.mutator2$popSummary[, "NumClones"] ,
##                  m1.mutator1$popSummary[, "NumClones"], alternative = "greater")$p.value < p.fail)
##     expect_true( t.test( m1.mutator1$popSummary[, "NumClones"] ,
##                  m1.mutator0$popSummary[, "NumClones"], alternative = "greater")$p.value < p.fail)
##     expect_true( t.test( m2.mutator2$popSummary[, "NumClones"] ,
##                  m2.mutator1$popSummary[, "NumClones"], alternative = "greater")$p.value < p.fail)
##     expect_true( t.test( m2.mutator1$popSummary[, "NumClones"] ,
##                  m2.mutator0$popSummary[, "NumClones"], alternative = "greater")$p.value < p.fail)
##     ## expect_true( mean(mutsPerCloneOSS(m1.mutator2)) >
##     ##              mean(mutsPerCloneOSS(m1.mutator1)))
##     ## expect_true( mean(mutsPerCloneOSS(m1.mutator1)) >
##     ##              mean(mutsPerCloneOSS(m1.mutator0)))
##     ## expect_true( mean(mutsPerCloneOSS(m2.mutator2)) >
##     ##              mean(mutsPerCloneOSS(m2.mutator1)))
##     ## expect_true( mean(mutsPerCloneOSS(m2.mutator1)) >
##     ##              mean(mutsPerCloneOSS(m2.mutator0)))
##     ## Increases in mutation rates increase clones, etc, within levels of
##     ## mutator.
##     expect_true( t.test( m1.mutator0$popSummary[, "NumClones"] , 
##                  m2.mutator0$popSummary[, "NumClones"], alternative = "less")$p.value < p.fail)
##     expect_true( t.test( m1.mutator1$popSummary[, "NumClones"] , 
##                  m2.mutator1$popSummary[, "NumClones"], alternative = "less")$p.value < p.fail)
##     expect_true( t.test( m1.mutator2$popSummary[, "NumClones"] , 
##                  m2.mutator2$popSummary[, "NumClones"], alternative = "less")$p.value < p.fail)
##     ## expect_true( mean(mutsPerCloneOSS(m1.mutator0)) <
##     ##              mean(mutsPerCloneOSS(m2.mutator0)))
##     ## expect_true( mean(mutsPerCloneOSS(m1.mutator1)) <
##     ##              mean(mutsPerCloneOSS(m2.mutator1)))
##     ## expect_true( mean(mutsPerCloneOSS(m1.mutator2)) <
##     ##              mean(mutsPerCloneOSS(m2.mutator2)))
    
## })







cat(paste("\n Finished test.mutator-oncoSimulSample.R test at", date()))



## singleCell. Stop on 1 driver, mark all except init as drivers.??? Nope,
## as driver present, but not abundant. So stop on two. Too convoluted. If
## I want to test sampling, test sampling. Period.





## converted from test.mutator using

## sed -i 's/median(summary(\([A-Za-z0-9]*\))$NumClones)/median(\1$popSummary\[, "NumClones"\])/g' test.mutator-oncoSimulSample.R
## sed -i  's/mutsPerClone(\([A-Za-z0-9]*\))/rowSums(\1$popSample)/g' test.mutator-oncoSimulSample.R
## sed -i 's/oncoSimulPop(/oncoSimulSample(/' test.mutator-oncoSimulSample.R
## sed -i 's/, mc.cores = 2//' test.mutator-oncoSimulSample.R
## sed -i 's/mc.cores = 2,//' test.mutator-oncoSimulSample.R
## sed -i 's/mc.cores = 2)/)/' test.mutator-oncoSimulSample.R

## sed -i 's/keepEvery = [0-9],//' test.mutator-oncoSimulSample.R
## sed -i 's/, keepEvery = [0-9]//' test.mutator-oncoSimulSample.R
## sed -i 's/keepEvery = [0-9])/)/' test.mutator-oncoSimulSample.R
## sed -i 's/summary(\([A-Za-z0-9]*\))/\1$popSummary\[, c(1:3, 8:9)\]/g' test.mutator-oncoSimulSample.R
## the last is not quite ok. Leaves to sets of the [, c(1:3, 8:9)][, c(1:3, 8:9)]. Replace in emacs.
## and a few others are missed. 
