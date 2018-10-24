## Repeat tests in test.mutator, using oncoSimulSample.
## This is a concession to extreme paranoia.


cat(paste("\n Starting test.mutator-oncoSimulSample-long.R test at", date()))
cat(paste("\n         a runif ", runif(1), "\n"))
## RNGkind("L'Ecuyer-CMRG") ## for the mclapplies


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

NClonesOSS <- function(x) {
    x$popSummary[, "NumClones"]
}



## ugly hack. Of course, not really mutations per clone. But the closest iwth oncoSimulSample and sampling whole pop.
mutsPerCloneOSS <- function(out) {
    rowSums(out$popSample)
}


p.value.threshold <- 0.005


## very slow, because huge number of clones. But tests several phenomena comprehensively.
## same with McFL below

date()
test_that("per-gene-mut rates and mutator", {
    max.tries <- 4
    for(tries in 1:max.tries) {
        T1 <- T2 <- T3 <- T4 <- T5 <- T6 <- T7 <- T8 <- TRUE
        cat("\n oss11-ossl: a runif is", runif(1), "\n")
        ng <- 40
        ni <- rep(0, ng)
        m1 <- runif(ng, min = 1e-7, max = 5e-6)
        m2 <- rep(1e-5, ng) ## runif(ng, min = 1e-5, max = 1e-4): ## crazy num of clones
        names(ni) <- names(m1) <- names(m2) <- c(replicate(ng,
                                                           paste(sample(letters, 12), collapse = "")))
        fe1 <- allFitnessEffects(noIntGenes = ni)
        ft <- 25 ## 50 this is crazy and takes forever
        no <- 5e5 
        reps <- 40
        gn <- paste(names(ni), collapse = ", ")
        mutator1 <- allMutatorEffects(epistasis = c("MU" = 20),
                                      geneToModule = c("MU" = gn))
        mutator2 <- allMutatorEffects(epistasis = c("MU" = 40),
                                      geneToModule = c("MU" = gn))
        m1.mutator0 <- oncoSimulSample(reps,
                                       fe1,
                                       mu = m1,
                                       onlyCancer = FALSE, detectionProb = NA,
                                       initSize = no,
                                       finalTime = ft,
                                       detectionSize = 1e9,
                                       detectionDrivers = 9999,
                                       sampleEvery = 0.01, thresholdWhole = 1e-20,
                                       seed = NULL, max.wall.time = 2000
                                        )
        m1.mutator1 <- oncoSimulSample(reps,
                                       fe1,
                                       mu = m1,
                                       muEF = mutator1,
                                       onlyCancer = FALSE, detectionProb = NA,
                                       initSize = no,
                                       finalTime = ft,
                                       detectionSize = 1e9,
                                       detectionDrivers = 9999,
                                       sampleEvery = 0.01, thresholdWhole = 1e-20,
                                       seed = NULL, max.wall.time = 2000
                                       )
        m1.mutator2 <- oncoSimulSample(reps,
                                       fe1,
                                       mu = m1,
                                       muEF = mutator2,
                                       onlyCancer = FALSE, detectionProb = NA,
                                       initSize = no,
                                       finalTime = ft,
                                       detectionSize = 1e9,
                                       detectionDrivers = 9999,
                                       sampleEvery = 0.01, thresholdWhole = 1e-20,
                                       seed = NULL, max.wall.time = 2000
                                       )
        runif(1)
        m2.mutator0 <- oncoSimulSample(reps,
                                       fe1,
                                       mu = m2,
                                       onlyCancer = FALSE, detectionProb = NA,
                                       initSize = no,
                                       finalTime = ft,
                                       detectionSize = 1e9,
                                       detectionDrivers = 9999,
                                       sampleEvery = 0.01, thresholdWhole = 1e-20,
                                       seed = NULL, max.wall.time = 2000
                                       )
        m2.mutator1 <- oncoSimulSample(reps,
                                       fe1,
                                       mu = m2,
                                       muEF = mutator1,
                                       onlyCancer = FALSE, detectionProb = NA,
                                       initSize = no,
                                       finalTime = ft,
                                       detectionSize = 1e9,
                                       detectionDrivers = 9999,
                                       sampleEvery = 0.01, thresholdWhole = 1e-20,
                                       seed = NULL, max.wall.time = 2000
                                       )
        m2.mutator2 <- oncoSimulSample(reps,
                                       fe1,
                                       mu = m2,
                                       muEF = mutator2,
                                       onlyCancer = FALSE, detectionProb = NA,
                                       initSize = no,
                                       finalTime = ft,
                                       detectionSize = 1e9,
                                       detectionDrivers = 9999,
                                       sampleEvery = 0.01, thresholdWhole = 1e-20,
                                       seed = NULL, max.wall.time = 2000
                                       )
        if(! (
            inherits(m1.mutator0$popSummary, "data.frame") &&
            inherits(m1.mutator1$popSummary, "data.frame") &&
            inherits(m1.mutator2$popSummary, "data.frame") &&
            inherits(m2.mutator0$popSummary, "data.frame") &&
            inherits(m2.mutator1$popSummary, "data.frame") &&
            inherits(m2.mutator1$popSummary, "data.frame") ) ) {
            T8 <- FALSE
            cat("\n     not a data frame?\n")
        }
        if(T8) {
        m1.mutator0$popSummary[, c(1:3, 8:9)]
        m1.mutator1$popSummary[, c(1:3, 8:9)]
        m1.mutator2$popSummary[, c(1:3, 8:9)]
        m2.mutator0$popSummary[, c(1:3, 8:9)]
        m2.mutator1$popSummary[, c(1:3, 8:9)]
        m2.mutator2$popSummary[, c(1:3, 8:9)]
        
        ## Mutator increases if larger mutator and compared to no mutator
        ## within levels of per-gene mutation rates
        T1 <- ( wilcox.test(m1.mutator2$popSummary[, "NumClones"], m1.mutator1$popSummary[, "NumClones"],
                                 alternative = "greater")$p.value < p.value.threshold)
        T2 <- ( wilcox.test(m1.mutator1$popSummary[, "NumClones"], m1.mutator0$popSummary[, "NumClones"],
                                 alternative = "greater")$p.value < p.value.threshold)
        T3 <- ( wilcox.test(m2.mutator2$popSummary[, "NumClones"], m2.mutator1$popSummary[, "NumClones"],
                                 alternative = "greater")$p.value < p.value.threshold)
        T4 <- ( wilcox.test(m2.mutator1$popSummary[, "NumClones"], m2.mutator0$popSummary[, "NumClones"],
                                 alternative = "greater")$p.value < p.value.threshold)
        ## Increases in mutation rates increase clones, etc, within levels of
        ## mutator.
        T5 <- ( wilcox.test(m2.mutator0$popSummary[, "NumClones"], m1.mutator0$popSummary[, "NumClones"],
                                 alternative = "greater")$p.value < p.value.threshold)
        T6 <- ( wilcox.test(m2.mutator1$popSummary[, "NumClones"], m1.mutator1$popSummary[, "NumClones"],
                                 alternative = "greater")$p.value < p.value.threshold)
        T7 <- ( wilcox.test(m2.mutator2$popSummary[, "NumClones"], m1.mutator2$popSummary[, "NumClones"],
                                 alternative = "greater")$p.value < p.value.threshold)
        ## expect_true( mean(mutsPerCloneOSS(m1.mutator0)) <
        ##              mean(mutsPerCloneOSS(m2.mutator0)))
        ## expect_true( mean(mutsPerCloneOSS(m1.mutator1)) <
        ##              mean(mutsPerCloneOSS(m2.mutator1)))
        ## expect_true( mean(mutsPerCloneOSS(m1.mutator2)) <
        ##              mean(mutsPerCloneOSS(m2.mutator2)))
        }
        if(T1 && T2 && T3 && T4 && T5 && T6 && T7 && T8) break;
    }
    cat(paste("\n done tries", tries, "\n"))
    expect_true(T1 && T2 && T3 && T4 && T5 && T6 && T7 && T8)

})
date()

date()
test_that("McFL: per-gene-mut rates and mutator", {
    max.tries <- 4
    for(tries in 1:max.tries) {
        T1 <- T2 <- T3 <- T4 <- T5 <- T6 <- T7 <- T8 <- TRUE
        cat("\n mcfloss11-ossl: a runif is", runif(1), "\n")
        ng <- 40
        ni <- rep(0, ng)
        m1 <- rep(5e-6, ng) ## too much variation and hard to pick the diffs.;
        ## runif(ng, min = 1e-7, max = 5e-6) And if too
        ## tiny, you do not pick them up unless huge ft
        ## and then it is way too slow for m2, etc.
        m2 <- rep(5e-5, ng)
        names(ni) <- names(m1) <- names(m2) <- c(replicate(ng,
                                                           paste(sample(letters, 12), collapse = "")))
        fe1 <- allFitnessEffects(noIntGenes = ni)
        ft <- 20 ## 50
        no <- 5e5 
        reps <- 20 ## 40
        gn <- paste(names(ni), collapse = ", ")
        mutator1 <- allMutatorEffects(epistasis = c("MU" = 20),
                                      geneToModule = c("MU" = gn))
        mutator2 <- allMutatorEffects(epistasis = c("MU" = 40),
                                      geneToModule = c("MU" = gn))
        m1.mutator0 <- oncoSimulSample(reps,
                                       fe1,
                                       mu = m1,
                                       onlyCancer = FALSE, detectionProb = NA,
                                       initSize = no,
                                       finalTime = ft,
                                       detectionSize = 1e9,
                                       detectionDrivers = 9999,
                                       sampleEvery = 0.01, thresholdWhole = 1e-20,
                                       seed = NULL, max.wall.time = 2000, model = "McFL"
                                       )
        m1.mutator1 <- oncoSimulSample(reps,
                                       fe1,
                                       mu = m1,
                                       muEF = mutator1,
                                       onlyCancer = FALSE, detectionProb = NA,
                                       initSize = no,
                                       finalTime = ft,
                                       detectionSize = 1e9,
                                       detectionDrivers = 9999,
                                       sampleEvery = 0.01, thresholdWhole = 1e-20,
                                       seed = NULL, max.wall.time = 2000, model = "McFL"
                                       )
        m1.mutator2 <- oncoSimulSample(reps,
                                       fe1,
                                       mu = m1,
                                       muEF = mutator2,
                                       onlyCancer = FALSE, detectionProb = NA,
                                       initSize = no,
                                       finalTime = ft,
                                       detectionSize = 1e9,
                                       detectionDrivers = 9999,
                                       sampleEvery = 0.01, thresholdWhole = 1e-20,
                                       seed = NULL, max.wall.time = 2000, model = "McFL"
                                       )
        cat("\n starting m2\n")
        m2.mutator0 <- oncoSimulSample(reps,
                                       fe1,
                                       mu = m2,
                                       onlyCancer = FALSE, detectionProb = NA,
                                       initSize = no,
                                       finalTime = ft,
                                       detectionSize = 1e9,
                                       detectionDrivers = 9999,
                                       sampleEvery = 0.01, thresholdWhole = 1e-20,
                                       seed = NULL, max.wall.time = 2000, model = "McFL"
                                       )
        m2.mutator1 <- oncoSimulSample(reps,
                                       fe1,
                                       mu = m2,
                                       muEF = mutator1,
                                       onlyCancer = FALSE, detectionProb = NA,
                                       initSize = no,
                                       finalTime = ft,
                                       detectionSize = 1e9,
                                       detectionDrivers = 9999,
                                       sampleEvery = 0.01, thresholdWhole = 1e-20,
                                       seed = NULL, max.wall.time = 2000, model = "McFL"
                                       )
        m2.mutator2 <- oncoSimulSample(reps,
                                       fe1,
                                       mu = m2,
                                       muEF = mutator2,
                                       onlyCancer = FALSE, detectionProb = NA,
                                       initSize = no,
                                       finalTime = ft,
                                       detectionSize = 1e9,
                                       detectionDrivers = 9999,     
                                       sampleEvery = 0.01, thresholdWhole = 1e-20,
                                       seed = NULL, max.wall.time = 2000, model = "McFL"
                                       )
        
        if(! (
            inherits(m1.mutator0$popSummary, "data.frame") &&
            inherits(m1.mutator1$popSummary, "data.frame") &&
            inherits(m1.mutator2$popSummary, "data.frame") &&
            inherits(m2.mutator0$popSummary, "data.frame") &&
            inherits(m2.mutator1$popSummary, "data.frame") &&
            inherits(m2.mutator1$popSummary, "data.frame") ) ) {
            cat("\n     not a data frame?\n")
            T8 <- FALSE
        }
        if(T8) {
        m1.mutator0$popSummary[, c(1:3, 8:9)]
        m1.mutator1$popSummary[, c(1:3, 8:9)]
        m1.mutator2$popSummary[, c(1:3, 8:9)]
        m2.mutator0$popSummary[, c(1:3, 8:9)]
        m2.mutator1$popSummary[, c(1:3, 8:9)]
        m2.mutator2$popSummary[, c(1:3, 8:9)]
        ## Mutator increases if larger mutator and compared to no mutator
        ## within levels of per-gene mutation rates
        ## we could use wilcoxon or t tests actually, specially because often diffs
        ## are not huge.
        p.fail <- 0.005
        T1 <- ( t.test( m1.mutator2$popSummary[, "NumClones"] ,
                            m1.mutator1$popSummary[, "NumClones"], alternative = "greater")$p.value < p.fail)
        T2 <- ( t.test( m1.mutator1$popSummary[, "NumClones"] ,
                            m1.mutator0$popSummary[, "NumClones"], alternative = "greater")$p.value < p.fail)
        T3 <- ( t.test( m2.mutator2$popSummary[, "NumClones"] ,
                            m2.mutator1$popSummary[, "NumClones"], alternative = "greater")$p.value < p.fail)
        T4 <- ( t.test( m2.mutator1$popSummary[, "NumClones"] ,
                            m2.mutator0$popSummary[, "NumClones"], alternative = "greater")$p.value < p.fail)
        ## Increases in mutation rates increase clones, etc, within levels of
        ## mutator.
        T5 <- ( t.test( m1.mutator0$popSummary[, "NumClones"] , 
                            m2.mutator0$popSummary[, "NumClones"], alternative = "less")$p.value < p.fail)
        T6 <- ( t.test( m1.mutator1$popSummary[, "NumClones"] , 
                            m2.mutator1$popSummary[, "NumClones"], alternative = "less")$p.value < p.fail)
        T7 <- ( t.test( m1.mutator2$popSummary[, "NumClones"] , 
                       m2.mutator2$popSummary[, "NumClones"], alternative = "less")$p.value < p.fail)
        }
        if(T1 && T2 && T3 && T4 && T5 && T6 && T7 && T8) break;
    }
    cat(paste("\n done tries", tries, "\n"))
    expect_true(T1 && T2 && T3 && T4 && T5 && T6 && T7 && T8)

})
date()


date()
test_that("Mutator increases by given factor with per-gene-mut rates: major axis and chi-sq test", {
    max.tries <- 4
    for(tries in 1:max.tries) {
        T1 <- T2 <- T3 <- T4 <- T5 <- T6 <- T7 <- T8 <- TRUE
        ## Two cases: mutator and no mutator, with variable mutation rates.
        ## rates such that rates of no mutator = rates of mutator * mutator.
        ## Why not compare mutlitplication factor keeping mutation rates
        ## constant? Because specially with mutator and large diffs in mut
        ## rates, with oncoSimulSample you undersample variation with
        ## wholePop, etc.
        ## Setings similar to oss11 in per-gene-mutation-rates but with the mutator
        cat("\n AEu8_long-ossl: a runif is", runif(1), "\n")
        pops <- 8000
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
                                    onlyCancer = FALSE, detectionProb = NA, seed = NULL)
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
        T1 <- (chisq.test(snomSampl("oreoisasabgene", m1.pg1.b),
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
                                    onlyCancer = FALSE, detectionProb = NA, seed = NULL)
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
        T2 <- (chisq.test(snomSampl("oreoisasabgene", m1.pg2.b),
                               p = pnom("oreoisasabgene", pg2, no, pops))$p.value > p.fail)
        ## Compare the mutator with the no mutator
        T3 <- (chisq.test(snomSampl("oreoisasabgene", m1.pg1.b),
                               snomSampl("oreoisasabgene", m1.pg2.b))$p.value > p.fail)
        y <- sqrt(snomSampl("oreoisasabgene", m1.pg1.b))
        x <- sqrt(snomSampl("oreoisasabgene", m1.pg2.b))
        mma <- smatr::ma(y ~ x, slope.test = 1, elev.test = 0) ## From smatr package, for major axis
        ## intercept not different from 0
        T4 <- (mma$elevtest[[1]]$p > p.fail)
        T5 <- (mma$slopetest[[1]]$p > p.fail)
        if(T1 && T2 && T3 && T4 && T5 && T6 && T7 && T8) break;
    }
    cat(paste("\n done tries", tries, "\n"))
    expect_true(T1 && T2 && T3 && T4 && T5 && T6 && T7 && T8)

})
date()



date()
test_that("McFL: Mutator increases by given factor with per-gene-mut rates: major axis and chi-sq test", {
    max.tries <- 4
    for(tries in 1:max.tries) {
        T1 <- T2 <- T3 <- T4 <- T5 <- T6 <- T7 <- T8 <- TRUE
        ## Two cases: mutator and no mutator, with variable mutation rates.
        ## rates such that rates of no mutator = rates of mutator * mutator.
        ## Why not compare mutlitplication factor keeping mutation rates
        ## constant? Because specially with mutator and large diffs in mut
        ## rates, with oncoSimulSample you undersample variation with
        ## wholePop, etc.
        ## Setings similar to oss11 in per-gene-mutation-rates but with the mutator
        cat("\n MCFL: long_AEu8-ossl: a runif is", runif(1), "\n")
        pops <- 8000
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
                                    onlyCancer = FALSE, detectionProb = NA, seed = NULL)
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
        T1 <- (chisq.test(snomSampl("oreoisasabgene", m1.pg1.b),
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
                                    onlyCancer = FALSE, detectionProb = NA, seed = NULL)
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
        T2 <- (chisq.test(snomSampl("oreoisasabgene", m1.pg2.b),
                               p = pnom("oreoisasabgene", pg2, no, pops))$p.value > p.fail)
        ## Compare mutator with no mutator
        T3 <- (chisq.test(snomSampl("oreoisasabgene", m1.pg1.b),
                               snomSampl("oreoisasabgene", m1.pg2.b))$p.value > p.fail)
        y <- sqrt(snomSampl("oreoisasabgene", m1.pg1.b))
        x <- sqrt(snomSampl("oreoisasabgene", m1.pg2.b))
        mma <- smatr::ma(y ~ x, slope.test = 1, elev.test = 0) ## From smatr package, for major axis
        ## intercept not different from 0
        T4 <- (mma$elevtest[[1]]$p > p.fail)
        T5 <- (mma$slopetest[[1]]$p > p.fail)
        if(T1 && T2 && T3 && T4 && T5 && T6 && T7 && T8) break;
    }
    cat(paste("\n done tries", tries, "\n"))
    expect_true(T1 && T2 && T3 && T4 && T5 && T6 && T7 && T8)

})
date()

date()
test_that("Mutator, several modules differences", {
    max.tries <- 4
    for(tries in 1:max.tries) {
        T1 <- T2 <- T3 <- T4 <- T5 <- T6 <- T7 <- T8 <- TRUE
        cat("\n mmd1_2-ossl: a runif is", runif(1), "\n")
        reps <- 140
        no <- 5e3
        ft <- 80 ## you need it large enough to get enough hits
        mu <- 1e-5
        ln <- 50 
        m1 <- 7 ## if this is too large, easy to get it to blow.
        ni <- rep(0, 2 * ln)
        gna <- paste0("a", 1:ln)
        gnb <- paste0("b", 1:ln)
        names(ni) <- c(gna, gnb)
        gn1 <- paste(c(gna, gnb), collapse = ", ")
        gna <- paste(gna, collapse = ", ")
        gnb <- paste(gnb, collapse = ", ")
        mut1 <- allMutatorEffects(epistasis = c("A" = m1),
                                  geneToModule = c("A" = gn1))
        mut2 <- allMutatorEffects(epistasis = c("A" = m1,
                                                "B" = m1),
                                  geneToModule = c("A" = gna,
                                                   "B" = gnb))
        f1 <- allFitnessEffects(noIntGenes = ni)
        b1 <- oncoSimulSample(reps,
                              f1,
                              mu = mu,
                              muEF = mut1,
                              onlyCancer = FALSE, detectionProb = NA,
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
                              onlyCancer = FALSE, detectionProb = NA,
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
        T1 <- ( wilcox.test( b2$popSummary[, "NumClones"], 
                                 b1$popSummary[, "NumClones"], alternative = "greater")$p.value < p.fail)
        T2 <- ( wilcox.test( rowSums(b2$popSample) ,
                                 rowSums(b1$popSample), alternative = "greater")$p.value < p.fail)
        if(T1 && T2 && T3 && T4 && T5 && T6 && T7 && T8) break;
    }
    cat(paste("\n done tries", tries, "\n"))
    expect_true(T1 && T2 && T3 && T4 && T5 && T6 && T7 && T8)

})
date()


## Remember that numClones is underestimated, possibly severly, by
## oncoSimulSample compared to oncoSimulPop, since we only look at the
## clones that exist at the end.

date()
test_that("Mutator, several modules differences, McFL", {
    max.tries <- 4
    for(tries in 1:max.tries) {
        T1 <- T2 <- T3 <- T4 <- T5 <- T6 <- T7 <- T8 <- TRUE
        cat("\n mmd1_2_mc-ossl: a runif is", runif(1), "\n")
        reps <- 80
        no <- 5e3
        ft <- 50 ## you need it large enough to get enough hits
        mu <- 1e-5
        ln <- 50 
        m1 <- 7 ## if this is too large, easy to get it to blow.
        ni <- rep(0, 2 * ln)
        gna <- paste0("a", 1:ln)
        gnb <- paste0("b", 1:ln)
        names(ni) <- c(gna, gnb)
        gn1 <- paste(c(gna, gnb), collapse = ", ")
        gna <- paste(gna, collapse = ", ")
        gnb <- paste(gnb, collapse = ", ")
        mut1 <- allMutatorEffects(epistasis = c("A" = m1),
                                  geneToModule = c("A" = gn1))
        mut2 <- allMutatorEffects(epistasis = c("A" = m1,
                                                "B" = m1),
                                  geneToModule = c("A" = gna,
                                                   "B" = gnb))
        f1 <- allFitnessEffects(noIntGenes = ni)
        b1 <- oncoSimulSample(reps,
                              f1,
                              mu = mu,
                              muEF = mut1,
                              onlyCancer = FALSE, detectionProb = NA,
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
                              onlyCancer = FALSE, detectionProb = NA,
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
        T1 <- ( t.test( b2$popSummary[, "NumClones"], 
                            b1$popSummary[, "NumClones"], alternative = "greater")$p.value < p.fail)
        T2 <- ( t.test( rowSums(b2$popSample) ,
                            rowSums(b1$popSample), alternative = "greater")$p.value < p.fail)
        if(T1 && T2 && T3 && T4 && T5 && T6 && T7 && T8) break;
    }
    cat(paste("\n done tries", tries, "\n"))
    expect_true(T1 && T2 && T3 && T4 && T5 && T6 && T7 && T8)

})
date()


## Slow (~ 3 seconds) but tests modules of mutator nicely.

date() ## Beware: this uses a lot of RAM without the gc()
test_that("Mutator modules differences", {
    max.tries <- 4
    for(tries in 1:max.tries) {
        T1 <- T2 <- T3 <- T4 <- T5 <- T6 <- T7 <- T8 <- TRUE
        cat("\n mmd1-ossl: a runif is", runif(1), "\n")
        reps <- 40
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
                              onlyCancer = FALSE, detectionProb = NA,
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
                              onlyCancer = FALSE, detectionProb = NA,
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
                              onlyCancer = FALSE, detectionProb = NA,
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
        T1 <- ( t.test( b3$popSummary[, "NumClones"], 
                            b2$popSummary[, "NumClones"], alternative = "greater")$p.value < p.fail)
        T2 <- ( t.test( b2$popSummary[, "NumClones"], 
                            b1$popSummary[, "NumClones"], alternative = "greater")$p.value < p.fail)
        T3 <- ( t.test( rowSums(b3$popSample) ,
                            rowSums(b2$popSample), alternative = "greater")$p.value < p.fail)
        T4 <- ( t.test( rowSums(b2$popSample) ,
                            rowSums(b1$popSample), alternative = "greater")$p.value < p.fail)
        if(T1 && T2 && T3 && T4 && T5 && T6 && T7 && T8) break;
    }
    cat(paste("\n done tries", tries, "\n"))
    expect_true(T1 && T2 && T3 && T4 && T5 && T6 && T7 && T8)

})
date()

date()
test_that("Relative ordering of number of clones with mutator effects", {
    max.tries <- 4
    for(tries in 1:max.tries) {
        T1 <- T2 <- T3 <- T4 <- T5 <- T6 <- T7 <- T8 <- TRUE
        cat("\n x1-ossl: a runif is", runif(1), "\n")
        pops <- 40
        fe <- allFitnessEffects(noIntGenes = c("a" = 0.12,
                                               "b" = 0.14,
                                               "c" = 0.16,
                                               "d" = 0.11))
        fm6 <- allMutatorEffects(noIntGenes = c("a" = 5,
                                                "b" = 5,
                                                "c" = 5,
                                                "d" = 5))
        nc1 <- oncoSimulSample(pops, fe, muEF = fm6, finalTime =250,
                               mutationPropGrowth = FALSE,
                               sampleEvery = 0.01, thresholdWhole = 1e-20,
                               initSize = 1e6,
                               detectionSize = 1e8,
                               detectionDrivers = 9999, seed = NULL,
                               onlyCancer = FALSE, detectionProb = NA)
        fm8 <- allMutatorEffects(noIntGenes = c("a" = 1,
                                                "b" = 1,
                                                "c" = 1,
                                                "d" = 1))
        nc2 <- oncoSimulSample(pops, fe, muEF = fm8, finalTime =250,
                               mutationPropGrowth = FALSE,
                               sampleEvery = 0.01, thresholdWhole = 1e-20,
                               initSize = 1e6,
                               detectionSize = 1e8,
                               detectionDrivers = 9999, seed = NULL,
                               onlyCancer = FALSE, detectionProb = NA)
        fm7 <- allMutatorEffects(noIntGenes = c("a" = 1e-3,
                                                "b" = 1e-3,
                                                "c" = 1e-3,
                                                "d" = 1e-3))
        nc3 <- oncoSimulSample(pops, fe, muEF = fm7, finalTime =250,
                               mutationPropGrowth = FALSE,
                               sampleEvery = 0.01, thresholdWhole = 1e-20,
                               initSize = 1e6,
                               detectionSize = 1e8,
                               detectionDrivers = 9999, seed = NULL,
                               onlyCancer = FALSE, detectionProb = NA)
        T1 <- (wilcox.test(nc1$popSummary[, "NumClones"], nc2$popSummary[, "NumClones"],
                                alternative = "greater")$p.value < p.value.threshold)
        T2 <- (wilcox.test(nc2$popSummary[, "NumClones"], nc3$popSummary[, "NumClones"],
                                alternative = "greater")$p.value < p.value.threshold)
        ## rowSums cannot work as all have all muts.
        ## expect_true(t.test(rowSums(nc1$popSample),rowSums(nc2$popSample),
        ##                    alternative = "greater")$p.value < p.value.threshold)
        ## expect_true(t.test(rowSums(nc2$popSample),rowSums(nc3$popSample),
        ##                    alternative = "greater")$p.value < p.value.threshold)
        nc1$popSummary[, c(1:3, 8:9)]
        nc2$popSummary[, c(1:3, 8:9)]
        nc3$popSummary[, c(1:3, 8:9)]
        if(T1 && T2 && T3 && T4 && T5 && T6 && T7 && T8) break;
    }
    cat(paste("\n done tries", tries, "\n"))
    expect_true(T1 && T2 && T3 && T4 && T5 && T6 && T7 && T8)

})
date()


date()
test_that("Expect freq genotypes, mutator and var mut rates", {
    max.tries <- 4
    for(tries in 1:max.tries) {
        T1 <- T2 <- T3 <- T4 <- T5 <- T6 <- T7 <- T8 <- TRUE
        ## Similar to above, but mutator has a single element, not the whole
        ## vector.
        cat("\n u7-ossl: a runif is", runif(1), "\n")
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
                                    onlyCancer = FALSE, detectionProb = NA, seed = NULL)
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
        ## plot(snomSampl("oreoisasabgene", m1.pg1.b)/sum(snomSampl("oreoisasabgene", m1.pg1.b)) ~ 
        ##      pnom("oreoisasabgene", pg1, no, pops)); abline(a = 0, b = 1)
        ##      ## yes, if very large prob for one, it is slightly underestimated
        p.fail <- 1e-3
        T1 <- (chisq.test(snomSampl("oreoisasabgene", m1.pg1.b),
                               p = pnom("oreoisasabgene", pg1, no, pops))$p.value > p.fail)
        if(T1 && T2 && T3 && T4 && T5 && T6 && T7 && T8) break;
    }
    cat(paste("\n done tries", tries, "\n"))
    expect_true(T1 && T2 && T3 && T4 && T5 && T6 && T7 && T8)

})
date()


date()
test_that("Expect freq genotypes, mutator and var mut rates", {
    max.tries <- 4
    for(tries in 1:max.tries) {
        T1 <- T2 <- T3 <- T4 <- T5 <- T6 <- T7 <- T8 <- TRUE
        ## increase mutator, decrease max mu
        ## similar to oss11 in per-gene-mutation-rates but with the mutator
        cat("\n u8-ossl: a runif is", runif(1), "\n")
        pops <- 5000
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
                                    onlyCancer = FALSE, detectionProb = NA, seed = NULL)
        m1.pg1.b$popSummary[, c(1:3, 8:9)]
        ## If numclones is much larger than 2, that signals trouble as you are
        ## smoothing differences between frequencies with oncoSimulSample,
        ## whole pop
        summary(m1.pg1.b$popSummary[, "NumClones"])
        ## Recall that init-mutant tests check always present of initMutant
        ## against a thresholWhole of 1. Here it is slightly different.
        expect_true(smSampl("oreoisasabgene", m1.pg1.b) == pops)
        ## catch a pattern that would make the previous trivially true
        expect_false(sum(m1.pg1.b$popSample) == pops * (lni + 3))
        ## next two, to compare with oss1a
        sort(enom("oreoisasabgene", pg1, no, pops))
        sort(snomSampl("oreoisasabgene", m1.pg1.b))
        ## pnom("oreoisasabgene", pg1, no, pops)
        ## snomSampl("oreoisasabgene", m1.pg1.b)
        ## plot(snomSampl("oreoisasabgene", m1.pg1.b)/sum(snomSampl("oreoisasabgene", m1.pg1.b)) ~ 
        ##      pnom("oreoisasabgene", pg1, no, pops)); abline(a = 0, b = 1)
        ##      ## yes, if very large prob for one, it is slightly underestimated
        p.fail <- 1e-3
        T1 <- (chisq.test(snomSampl("oreoisasabgene", m1.pg1.b),
                               p = pnom("oreoisasabgene", pg1, no, pops))$p.value > p.fail)
        if(T1 && T2 && T3 && T4 && T5 && T6 && T7 && T8) break;
    }
    cat(paste("\n done tries", tries, "\n"))
    expect_true(T1 && T2 && T3 && T4 && T5 && T6 && T7 && T8)

})
date()


date()
test_that("McFL, Expect freq genotypes, mutator and var mut rates", {
    max.tries <- 4
    for(tries in 1:max.tries) {
        T1 <- T2 <- T3 <- T4 <- T5 <- T6 <- T7 <- T8 <- TRUE
        ## We test that mutator does not affect expected frequencies of
        ## mutated genes: they are given by the mutation rate of each gene.
        cat("\n mcfu6-ossl: a runif is", runif(1), "\n")
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
                                    onlyCancer = FALSE, detectionProb = NA, seed = NULL)
        summary(m1.pg1.b$popSummary[, "NumClones"])
        ## m1.pg1.b$popSummary[, c(1:3, 8:9)]
        expect_true(smSampl("oreoisasabgene", m1.pg1.b) == pops)
        enom("oreoisasabgene", pg1, no, pops)
        snomSampl("oreoisasabgene", m1.pg1.b)
        ## plot(snomSampl("oreoisasabgene", m1.pg1.b)/sum(snomSampl("oreoisasabgene", m1.pg1.b)) ~ 
        ##      pnom("oreoisasabgene", pg1, no, pops)); abline(a = 0, b = 1)
        p.fail <- 1e-3
        T1 <- (chisq.test(snomSampl("oreoisasabgene", m1.pg1.b),
                               p = pnom("oreoisasabgene", pg1, no, pops))$p.value > p.fail)
        if(T1 && T2 && T3 && T4 && T5 && T6 && T7 && T8) break;
    }
    cat(paste("\n done tries", tries, "\n"))
    expect_true(T1 && T2 && T3 && T4 && T5 && T6 && T7 && T8)

})
date()


date()
test_that("MCFL Relative ordering of number of clones with mutator effects", {
    max.tries <- 4
    for(tries in 1:max.tries) {
        T1 <- T2 <- T3 <- T4 <- T5 <- T6 <- T7 <- T8 <- TRUE
        ## Stop on size, and do a fast model (stop shortly after)
        ## But note we have only four positions left, so difficult to detect
        cat("\n mcx1-ossl-ossl: a runif is", runif(1), "\n")
        pops <- 160
        mu <- 1e-6
        ft <- 500
        fe <- allFitnessEffects(noIntGenes = c("a" = 0.11,
                                               "b" = 0.11,
                                               "c" = 0.11,
                                               "d" = 0.11))
        fm6 <- allMutatorEffects(noIntGenes = c("a" = 30,
                                                "b" = 30,
                                                "c" = 30,
                                                "d" = 30))
        nc1 <- oncoSimulSample(pops, mu = mu,
                               fe, muEF = fm6, finalTime = ft,
                               mutationPropGrowth = FALSE,
                               sampleEvery = 0.01, thresholdWhole = 1e-20,
                               initSize = 1e5, model = "McFL",
                               detectionSize = 1.11e5,
                               detectionDrivers = 9999, seed = NULL,
                               onlyCancer = FALSE, detectionProb = NA)
        fm8 <- allMutatorEffects(noIntGenes = c("a" = 2,
                                                "b" = 2,
                                                "c" = 2,
                                                "d" = 2))
        nc2 <- oncoSimulSample(pops, mu = mu,
                               fe, muEF = fm8, finalTime = ft,
                               mutationPropGrowth = FALSE,
                               sampleEvery = 0.01, thresholdWhole = 1e-20,
                               initSize = 1e5, model = "McFL",
                               detectionSize = 1.11e5,
                               detectionDrivers = 9999, seed = NULL,
                               onlyCancer = FALSE, detectionProb = NA)
        fm7 <- allMutatorEffects(noIntGenes = c("a" = 1e-3,
                                                "b" = 1e-3,
                                                "c" = 1e-3,
                                                "d" = 1e-3))
        nc3 <- oncoSimulSample(pops, mu = mu,
                               fe, muEF = fm7, finalTime = ft,
                               mutationPropGrowth = FALSE,
                               sampleEvery = 0.01, thresholdWhole = 1e-20,
                               initSize = 1e5, model = "McFL",
                               detectionSize = 1.11e5,
                               detectionDrivers = 9999, seed = NULL,
                               onlyCancer = FALSE, detectionProb = NA)
        expect_true(wilcox.test(nc1$popSummary[, "NumClones"], nc2$popSummary[, "NumClones"],
                                alternative = "greater")$p.value < p.value.threshold)
        expect_true(wilcox.test(nc2$popSummary[, "NumClones"], nc3$popSummary[, "NumClones"],
                                alternative = "greater")$p.value < p.value.threshold)
        nc1$popSummary[, c(1:3, 8:9)]
        nc2$popSummary[, c(1:3, 8:9)]
        nc3$popSummary[, c(1:3, 8:9)]
        T1 <- (t.test(rowSums(nc1$popSample),rowSums(nc2$popSample),
                           alternative = "greater")$p.value < p.value.threshold)
        T2 <- (t.test(rowSums(nc2$popSample),rowSums(nc3$popSample),
                           alternative = "greater")$p.value < p.value.threshold)
        if(T1 && T2 && T3 && T4 && T5 && T6 && T7 && T8) break;
    }
    cat(paste("\n done tries", tries, "\n"))
    expect_true(T1 && T2 && T3 && T4 && T5 && T6 && T7 && T8)

})
date()


date()
test_that("Relative ordering of number of clones with init mutant of mutator effects and s = 0", {
    max.tries <- 4
    for(tries in 1:max.tries) {
        T1 <- T2 <- T3 <- T4 <- T5 <- T6 <- T7 <- T8 <- TRUE
        ## Here stop by time, as s = 0
        cat("\n x2cd-ossl-ossl: a runif is", runif(1), "\n")
        pops <- 80
        ni <- rep(0, 50)
        names(ni) <- c("a", "b", "c", "d", paste0("n", 1:46))
        fe <- allFitnessEffects(noIntGenes = ni)
        fm6 <- allMutatorEffects(noIntGenes = c("a" = .05,
                                                "b" = 1,
                                                "c" = 10,
                                                "d" = 50))
        nca <- oncoSimulSample(pops, fe, muEF = fm6, finalTime =50,
                               mutationPropGrowth = FALSE,
                               initSize = 1e4,
                               initMutant = "a",
                               sampleEvery = 0.01, thresholdWhole = 1e-20,
                               detectionSize = 1e9,
                               detectionDrivers = 9999, seed = NULL,
                               onlyCancer = FALSE, detectionProb = NA)                       
        ncb <- oncoSimulSample(pops, fe, muEF = fm6, finalTime =50,
                               mutationPropGrowth = FALSE,
                               initSize = 1e4,
                               initMutant = "b",
                               sampleEvery = 0.01, thresholdWhole = 1e-20,
                               detectionSize = 1e9,
                               detectionDrivers = 9999, seed = NULL,
                               onlyCancer = FALSE, detectionProb = NA)                       
        ncc <- oncoSimulSample(pops, fe, muEF = fm6, finalTime =50,
                               mutationPropGrowth = FALSE,
                               initSize = 1e4,
                               initMutant = "c",
                               sampleEvery = 0.01, thresholdWhole = 1e-20,
                               detectionSize = 1e9,
                               detectionDrivers = 9999, seed = NULL,
                               onlyCancer = FALSE, detectionProb = NA)                       
        ncd <- oncoSimulSample(pops, fe, muEF = fm6, finalTime =50,
                               mutationPropGrowth = FALSE,
                               initSize = 1e4,
                               initMutant = "d",
                               sampleEvery = 0.01, thresholdWhole = 1e-20,
                               detectionSize = 1e9,
                               detectionDrivers = 9999, seed = NULL,
                               onlyCancer = FALSE, detectionProb = NA)                       
        ## These are the real tests
        T1 <- ( wilcox.test(nca$popSummary[, "NumClones"],
                                 ncb$popSummary[, "NumClones"],
                                 alternative = "less")$p.value < p.value.threshold)
        T2 <- (wilcox.test(ncb$popSummary[, "NumClones"],
                                ncc$popSummary[, "NumClones"],
                                alternative = "less")$p.value < p.value.threshold)
        T3 <- ( wilcox.test(ncc$popSummary[, "NumClones"],
                                 ncd$popSummary[, "NumClones"],
                                 alternative = "less")$p.value < p.value.threshold)
        T4 <- (t.test(rowSums(nca$popSample), rowSums(ncb$popSample),
                           alternative = "less")$p.value < p.value.threshold)
        T5 <- (t.test(rowSums(ncb$popSample), rowSums(ncc$popSample),
                           alternative = "less")$p.value < p.value.threshold)
        T6 <- (t.test(rowSums(ncc$popSample), rowSums(ncd$popSample),
                           alternative = "less")$p.value < p.value.threshold)
        if(T1 && T2 && T3 && T4 && T5 && T6 && T7 && T8) break;
    }
    cat(paste("\n done tries", tries, "\n"))
    expect_true(T1 && T2 && T3 && T4 && T5 && T6 && T7 && T8)

})
date()


date()
test_that("McFL, Expect freq genotypes, mutator and var mut rates, ct mut", {
    max.tries <- 4
    for(tries in 1:max.tries) {
        T1 <- T2 <- T3 <- T4 <- T5 <- T6 <- T7 <- T8 <- TRUE
        ## We test that mutator does not affect expected frequencies of
        ## mutated genes: they are given by the mutation rate of each gene.
        cat("\n ct_mcfu6-ossl: a runif is", runif(1), "\n")
        pops <- 2500
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
        pg1 <- runif(lni + 3, min = 1e-5, max = 1e-5) ## max should not be
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
                                    onlyCancer = FALSE, detectionProb = NA, seed = NULL)
        summary(m1.pg1.b$popSummary[, "NumClones"])
        ## m1.pg1.b$popSummary[, c(1:3, 8:9)]
        expect_true(smSampl("oreoisasabgene", m1.pg1.b) == pops)
        enom("oreoisasabgene", pg1, no, pops)
        snomSampl("oreoisasabgene", m1.pg1.b)
        ## plot(snomSampl("oreoisasabgene", m1.pg1.b)/sum(snomSampl("oreoisasabgene", m1.pg1.b)) ~ 
        ##      pnom("oreoisasabgene", pg1, no, pops)); abline(a = 0, b = 1)
        p.fail <- 1e-3
        T1 <- (chisq.test(snomSampl("oreoisasabgene", m1.pg1.b),
                               p = pnom("oreoisasabgene", pg1, no, pops))$p.value > p.fail)
        if(T1 && T2 && T3 && T4 && T5 && T6 && T7 && T8) break;
    }
    cat(paste("\n done tries", tries, "\n"))
    expect_true(T1 && T2 && T3 && T4 && T5 && T6 && T7 && T8)

})
date()


date()
test_that("McFL, Expect freq genotypes, mutator and var mut rates", {
    max.tries <- 4
    for(tries in 1:max.tries) {
        T1 <- T2 <- T3 <- T4 <- T5 <- T6 <- T7 <- T8 <- TRUE
        ## We test that mutator does not affect expected frequencies of
        ## mutated genes: they are given by the mutation rate of each gene.
        cat("\n mcfu7-ossl: a runif is", runif(1), "\n")
        pops <- 2500
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
                                    onlyCancer = FALSE, detectionProb = NA, seed = NULL)
        summary(m1.pg1.b$popSummary[, "NumClones"])
        ## m1.pg1.b$popSummary[, c(1:3, 8:9)]
        expect_true(smSampl("oreoisasabgene", m1.pg1.b) == pops)
        enom("oreoisasabgene", pg1, no, pops)
        snomSampl("oreoisasabgene", m1.pg1.b)
        ## plot(snomSampl("oreoisasabgene", m1.pg1.b)/sum(snomSampl("oreoisasabgene", m1.pg1.b)) ~ 
        ##      pnom("oreoisasabgene", pg1, no, pops)); abline(a = 0, b = 1)
        p.fail <- 1e-3
        T1 <- (chisq.test(snomSampl("oreoisasabgene", m1.pg1.b),
                               p = pnom("oreoisasabgene", pg1, no, pops))$p.value > p.fail)
        if(T1 && T2 && T3 && T4 && T5 && T6 && T7 && T8) break;
    }
    cat(paste("\n done tries", tries, "\n"))
    expect_true(T1 && T2 && T3 && T4 && T5 && T6 && T7 && T8)

})
date()


date()
test_that("Expect freq genotypes, mutator and var mut rates", {
    max.tries <- 4
    for(tries in 1:max.tries) {
        T1 <- T2 <- T3 <- T4 <- T5 <- T6 <- T7 <- T8 <- TRUE
        ## We test that mutator does not affect expected frequencies of
        ## mutated genes: they are given by the mutation rate of each gene.
        cat("\n sameu6-ossl: a runif is", runif(1), "\n")
        pops <- 1500
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
                                    onlyCancer = FALSE, detectionProb = NA, seed = NULL)
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
        ## plot(snomSampl("oreoisasabgene", m1.pg1.b)/sum(snomSampl("oreoisasabgene", m1.pg1.b)) ~ 
        ##      pnom("oreoisasabgene", pg1, no, pops)); abline(a = 0, b = 1)
        ## yes, if very large prob for one, it is slightly underestimated
        p.fail <- 1e-3
        T1 <- (chisq.test(snomSampl("oreoisasabgene", m1.pg1.b),
                               p = pnom("oreoisasabgene", pg1, no, pops))$p.value > p.fail)
        if(T1 && T2 && T3 && T4 && T5 && T6 && T7 && T8) break;
    }
    cat(paste("\n done tries", tries, "\n"))
    expect_true(T1 && T2 && T3 && T4 && T5 && T6 && T7 && T8)

})
date()


date()
test_that("Expect freq genotypes, mutator and var mut rates", {
    max.tries <- 4
    for(tries in 1:max.tries) {
        T1 <- T2 <- T3 <- T4 <- T5 <- T6 <- T7 <- T8 <- TRUE
        ## We test that mutator does not affect expected frequencies of
        ## mutated genes: they are given by the mutation rate of each gene.
        cat("\n u6-ossl: a runif is", runif(1), "\n")
        pops <- 1500
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
                                    onlyCancer = FALSE, detectionProb = NA, seed = NULL)
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
        ## plot(snomSampl("oreoisasabgene", m1.pg1.b)/sum(snomSampl("oreoisasabgene", m1.pg1.b)) ~ 
        ##      pnom("oreoisasabgene", pg1, no, pops)); abline(a = 0, b = 1)
        ## yes, if very large prob for one, it is slightly underestimated
        p.fail <- 1e-3
        T1 <- (chisq.test(snomSampl("oreoisasabgene", m1.pg1.b),
                               p = pnom("oreoisasabgene", pg1, no, pops))$p.value > p.fail)
        if(T1 && T2 && T3 && T4 && T5 && T6 && T7 && T8) break;
    }
    cat(paste("\n done tries", tries, "\n"))
    expect_true(T1 && T2 && T3 && T4 && T5 && T6 && T7 && T8)
})
date()


date()
test_that("MCFL Relative ordering of number of clones with init mutant of mutator effects and s = 0", {
 max.tries <- 4
    for(tries in 1:max.tries) {
        T1 <- T2 <- T3 <- T4 <- T5 <- T6 <- T7 <- T8 <- TRUE
    ## Stopping on time: all s = 0.
    cat("\n mcx2cd-ossl: a runif is", runif(1), "\n")
    pops <- 40
    ni <- rep(0, 50)
    names(ni) <- c("a", "b", "c", "d", paste0("n", 1:46))
    fe <- allFitnessEffects(noIntGenes = ni)
    fm6 <- allMutatorEffects(noIntGenes = c("a" = .05,
                                            "b" = 1,
                                            "c" = 10,
                                            "d" = 50))
    nca <- oncoSimulSample(pops, fe, muEF = fm6, finalTime =50,
                        mutationPropGrowth = FALSE,
                        initSize = 1e4,
                        initMutant = "a", detectionSize = 1e9,
                        sampleEvery = 0.01, thresholdWhole = 1e-20,
                        onlyCancer = FALSE, detectionProb = NA, model = "McFL")
    ncb <- oncoSimulSample(pops, fe, muEF = fm6, finalTime =50,
                        mutationPropGrowth = FALSE,
                        initSize = 1e4,
                        initMutant = "b", detectionSize = 1e9,
                        sampleEvery = 0.01, thresholdWhole = 1e-20,
                        onlyCancer = FALSE, detectionProb = NA, model = "McFL")
    ncc <- oncoSimulSample(pops, fe, muEF = fm6, finalTime =50,
                        mutationPropGrowth = FALSE,
                        initSize = 1e4,
                        initMutant = "c", detectionSize = 1e9,
                        sampleEvery = 0.01, thresholdWhole = 1e-20,
                        onlyCancer = FALSE, detectionProb = NA, model = "McFL")
    ncd <- oncoSimulSample(pops, fe, muEF = fm6, finalTime =50,
                        mutationPropGrowth = FALSE,
                        initSize = 1e4,
                        initMutant = "d",detectionSize = 1e9,
                        sampleEvery = 0.01, thresholdWhole = 1e-20,
                        onlyCancer = FALSE, detectionProb = NA, model = "McFL")
    ## These are the real tests
    T1 <- ( wilcox.test(nca$popSummary[, "NumClones"],
                             ncb$popSummary[, "NumClones"],
                             alternative = "less")$p.value < p.value.threshold)
    T2 <- (wilcox.test(ncb$popSummary[, "NumClones"],
                            ncc$popSummary[, "NumClones"],
                            alternative = "less")$p.value < p.value.threshold)
    T3 <- ( wilcox.test(ncc$popSummary[, "NumClones"],
                             ncd$popSummary[, "NumClones"],
                             alternative = "less")$p.value < p.value.threshold)
    T4 <- (t.test(rowSums(nca$popSample), rowSums(ncb$popSample),
                       alternative = "less")$p.value < p.value.threshold)
    T5 <- (t.test(rowSums(ncb$popSample), rowSums(ncc$popSample),
                       alternative = "less")$p.value < p.value.threshold)
    T6 <- (t.test(rowSums(ncc$popSample), rowSums(ncd$popSample),
                       alternative = "less")$p.value < p.value.threshold)
         if(T1 && T2 && T3 && T4 && T5 && T6 && T7 && T8) break;
    }
    cat(paste("\n done tries", tries, "\n"))
    expect_true(T1 && T2 && T3 && T4 && T5 && T6 && T7 && T8)

})
date()


date()
test_that("Relative ordering of number of clones with init mutant of mutator effects", {
    max.tries <- 4
    for(tries in 1:max.tries) {
        T1 <- T2 <- T3 <- T4 <- T5 <- T6 <- T7 <- T8 <- TRUE
        ## here we do not fill all positions, except maybe for ncd
        ## We stop on finalTime, not popsize
        cat("\n x2bc-ossl: a runif is", runif(1), "\n")
        pops <- 30
        ni <- rep(0.01, 50)
        names(ni) <- c("a", "b", "c", "d", paste0("n", 1:46))
        fe <- allFitnessEffects(noIntGenes = ni)
        fm6 <- allMutatorEffects(noIntGenes = c("a" = .05,
                                                "b" = 1,
                                                "c" = 10,
                                                "d" = 50))
        nca <- oncoSimulSample(pops, fe, muEF = fm6, finalTime =50,
                               mutationPropGrowth = FALSE,
                               initSize = 1e4,
                               initMutant = "a",
                               sampleEvery = 0.01, thresholdWhole = 1e-20,
                               detectionSize = 1e9,
                               detectionDrivers = 9999, seed = NULL,
                               onlyCancer = FALSE, detectionProb = NA)
        ncb <- oncoSimulSample(pops, fe, muEF = fm6, finalTime =50,
                               mutationPropGrowth = FALSE,
                               initSize = 1e4,
                               initMutant = "b",
                               sampleEvery = 0.01, thresholdWhole = 1e-20,
                               detectionSize = 1e9,
                               detectionDrivers = 9999, seed = NULL,
                               onlyCancer = FALSE, detectionProb = NA)
        ncc <- oncoSimulSample(pops, fe, muEF = fm6, finalTime =50,
                               mutationPropGrowth = FALSE,
                               initSize = 1e4,
                               initMutant = "c",
                               sampleEvery = 0.01, thresholdWhole = 1e-20,
                               detectionSize = 1e9,
                               detectionDrivers = 9999, seed = NULL,
                               onlyCancer = FALSE, detectionProb = NA)
        ncd <- oncoSimulSample(pops, fe, muEF = fm6, finalTime =50,
                               mutationPropGrowth = FALSE,
                               initSize = 1e4,
                               initMutant = "d",
                               sampleEvery = 0.01, thresholdWhole = 1e-20,
                               detectionSize = 1e9,
                               detectionDrivers = 9999, seed = NULL,
                               onlyCancer = FALSE, detectionProb = NA)
        T1 <- ( wilcox.test(nca$popSummary[, "NumClones"],
                                 ncb$popSummary[, "NumClones"],
                                 alternative = "less")$p.value < p.value.threshold)
        T2 <- (wilcox.test(ncb$popSummary[, "NumClones"],
                                ncc$popSummary[, "NumClones"],
                                alternative = "less")$p.value < p.value.threshold)
        T3 <- ( wilcox.test(ncc$popSummary[, "NumClones"],
                                 ncd$popSummary[, "NumClones"],
                                 alternative = "less")$p.value < p.value.threshold)
        T4 <- (t.test(rowSums(nca$popSample), rowSums(ncb$popSample),
                           alternative = "less")$p.value < p.value.threshold)
        T5 <- (t.test(rowSums(ncb$popSample), rowSums(ncc$popSample),
                           alternative = "less")$p.value < p.value.threshold)
        T6 <- (t.test(rowSums(ncc$popSample), rowSums(ncd$popSample),
                           alternative = "less")$p.value < p.value.threshold)
        nca$popSummary[, c(1:3, 8:9)]
        ncb$popSummary[, c(1:3, 8:9)]
        ncc$popSummary[, c(1:3, 8:9)]
        ncd$popSummary[, c(1:3, 8:9)]
        if(T1 && T2 && T3 && T4 && T5 && T6 && T7 && T8) break;
    }
    cat(paste("\n done tries", tries, "\n"))
    expect_true(T1 && T2 && T3 && T4 && T5 && T6 && T7 && T8)

})
date()


date()
test_that("Same mu vector, different mutator; diffs in number muts, larger t", {
    max.tries <- 4
    for(tries in 1:max.tries) {
        T1 <- T2 <- T3 <- T4 <- T5 <- T6 <- T7 <- T8 <- TRUE
        ## reproduction, death, and double and possibly triple mutants. We
        ## decrease init pop size to make this fast.
        cat("\n nm1-ossl: a runif is", runif(1), "\n")
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
                                 seed = NULL, onlyCancer = FALSE, detectionProb = NA)
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
                                  seed = NULL, onlyCancer = FALSE, detectionProb = NA)
        ## number of clones
        T1 <- (wilcox.test(NClonesOSS(pop10), NClonesOSS(pop100),
                                alternative = "less")$p.value < p.value.threshold)
        T2 <- (t.test(rowSums(pop10$popSample), rowSums(pop100$popSample),
                           alternative = "less")$p.value < p.value.threshold)
        if(T1 && T2 && T3 && T4 && T5 && T6 && T7 && T8) break;
    }
    cat(paste("\n done tries", tries, "\n"))
    expect_true(T1 && T2 && T3 && T4 && T5 && T6 && T7 && T8)

})
date()


date()
test_that("McFL: Same mu vector, different mutator; diffs in number muts, larger t", {
    max.tries <- 4
    for(tries in 1:max.tries) {
        T1 <- T2 <- T3 <- T4 <- T5 <- T6 <- T7 <- T8 <- TRUE
        ## reproduction, death, and double and possibly triple mutants. We
        ## decrease init pop size to make this fast.
        cat("\n nm3-ossl: a runif is", runif(1), "\n")
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
                                 seed = NULL, onlyCancer = FALSE, detectionProb = NA)
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
                                  seed = NULL, onlyCancer = FALSE, detectionProb = NA)
        T1 <- (wilcox.test(NClonesOSS(pop10), NClonesOSS(pop100),
                                alternative = "less")$p.value < p.value.threshold)
        T2 <- (t.test(rowSums(pop10$popSample), rowSums(pop100$popSample),
                           alternative = "less")$p.value < p.value.threshold)
        if(T1 && T2 && T3 && T4 && T5 && T6 && T7 && T8) break;
    }
    cat(paste("\n done tries", tries, "\n"))
    expect_true(T1 && T2 && T3 && T4 && T5 && T6 && T7 && T8)

})
date()


date() ## Beware: this uses a lot of RAM without the gc()
test_that("McFL: Mutator modules differences", {
    max.tries <- 4
    for(tries in 1:max.tries) {
        T1 <- T2 <- T3 <- T4 <- T5 <- T6 <- T7 <- T8 <- TRUE
        cat("\n MCFLmmd1-ossl: a runif is", runif(1), "\n")
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
                              onlyCancer = FALSE, detectionProb = NA,
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
                              onlyCancer = FALSE, detectionProb = NA,
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
                              onlyCancer = FALSE, detectionProb = NA,
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
        T1 <- ( t.test( b3$popSummary[, "NumClones"], 
                            b2$popSummary[, "NumClones"], alternative = "greater")$p.value < p.fail)
        T2 <- ( t.test( b2$popSummary[, "NumClones"], 
                            b1$popSummary[, "NumClones"], alternative = "greater")$p.value < p.fail)
        T3 <- ( t.test( rowSums(b3$popSample) ,
                            rowSums(b2$popSample), alternative = "greater")$p.value < p.fail)
        T4 <- ( t.test( rowSums(b2$popSample) ,
                            rowSums(b1$popSample), alternative = "greater")$p.value < p.fail)
        if(T1 && T2 && T3 && T4 && T5 && T6 && T7 && T8) break;
    }
    cat(paste("\n done tries", tries, "\n"))
    expect_true(T1 && T2 && T3 && T4 && T5 && T6 && T7 && T8)

})
date()


###### The following used to be in the non-long file. But that is just an overkill.

date()
test_that("Mutator increases by given factor with per-gene-mut rates: major axis and chi-sq test", {
    ## Two cases: mutator and no mutator, with variable mutation rates.
    ## rates such that rates of no mutator = rates of mutator * mutator.
    ## Why not compare mutlitplication factor keeping mutation rates
    ## constant? Because specially with mutator and large diffs in mut
    ## rates, with oncoSimulSample you undersample variation with
    ## wholePop, etc.
    ## Setings similar to oss11 in per-gene-mutation-rates but with the mutator
    max.tries <- 4
    for(tries in 1:max.tries) {
        T1 <- T2 <- T3 <- T4 <- T5 <- T6 <- T7 <- T8 <- TRUE
        cat("\n AEu8-ossl: a runif is", runif(1), "\n")
        pops <- 200
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
                                    onlyCancer = FALSE, detectionProb = NA, seed = NULL)
        ## m1.pg1.b$popSummary[, c(1:3, 8:9)]
        ## summary(m1.pg1.b$popSummary[, "NumClones"])
        ## Recall that init-mutant tests check always present of initMutant
        ## against a thresholWhole of 1. Here it is slightly different.
        expect_true(smSampl("oreoisasabgene", m1.pg1.b) == pops)
        ## catch a pattern that would make the previous trivially true
        expect_false(sum(m1.pg1.b$popSample) == pops * (lni + 3))
        ## next two, to compare with oss1a
        ## sort(enom("oreoisasabgene", pg1, no, pops))
        ## sort(snomSampl("oreoisasabgene", m1.pg1.b))
        ## Compare with the expected for this scenario
        p.fail <- 1e-3
        T1 <- (chisq.test(snomSampl("oreoisasabgene", m1.pg1.b),
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
                                    onlyCancer = FALSE, detectionProb = NA, seed = NULL)
        ## m1.pg2.b$popSummary[, c(1:3, 8:9)]
        ## summary(m1.pg2.b$popSummary[, "NumClones"])
        ## Recall that init-mutant tests check always present of initMutant
        ## against a thresholWhole of 1. Here it is slightly different.
        expect_true(smSampl("oreoisasabgene", m1.pg2.b) == pops)
        ## catch a pattern that would make the previous trivially true
        expect_false(sum(m1.pg2.b$popSample) == pops * (lni + 3))
        ## next two, to compare with oss1a
        ## sort(enom("oreoisasabgene", pg2, no, pops))
        ## sort(snomSampl("oreoisasabgene", m1.pg2.b))
        p.fail <- 1e-3
        T3 <- (chisq.test(snomSampl("oreoisasabgene", m1.pg2.b),
                          p = pnom("oreoisasabgene", pg2, no, pops))$p.value > p.fail)
        ## Compare the mutator with the no mutator
        T4 <- (chisq.test(snomSampl("oreoisasabgene", m1.pg1.b),
                          snomSampl("oreoisasabgene", m1.pg2.b))$p.value > p.fail)
        y <- sqrt(snomSampl("oreoisasabgene", m1.pg1.b))
        x <- sqrt(snomSampl("oreoisasabgene", m1.pg2.b))
        mma <- smatr::ma(y ~ x, slope.test = 1, elev.test = 0) ## From smatr package, for major axis
        ## intercept not different from 0
        T5 <- (mma$elevtest[[1]]$p > p.fail)
        T6 <- (mma$slopetest[[1]]$p > p.fail)
        if( T1 && T3 && T4 && T5 && T6) break;
    }
    cat(paste("\n done tries", tries, "\n"))
    expect_true(T1 && T2 && T3 && T4 && T5 && T6 && T7 && T8)
})
date()


date()
test_that("Same mu vector, different mutator; diffs in number muts, tiny t", {
    max.tries <- 4
    for(tries in 1:max.tries) {
        T1 <- T2 <- T3 <- T4 <- T5 <- T6 <- T7 <- T8 <- TRUE
        ## Here, there is no reproduction or death. Just mutation. And no double
        ## mutants either.
        ## We test:
        ##  - mutator increases mutation rates as seen in:
        ##        - number of clones created
        ##        - number of total mutation events
        cat("\n nm0-ossl: a runif is", runif(1), "\n")
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
                                 seed = NULL, onlyCancer = FALSE, detectionProb = NA)
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
                                  seed = NULL, onlyCancer = FALSE, detectionProb = NA)
        ## number of total mutations do not make sense with oncoSimulSample,
        ## since we cannot estimate them. we approximate wit sum of
        ## mutations. but that is too thick grain.  number of clones is much cleaner
        T1 <- (wilcox.test(NClonesOSS(pop10), NClonesOSS(pop100),
                                alternative = "less")$p.value < p.value.threshold)
        T2 <- (t.test(rowSums(pop10$popSample), rowSums(pop100$popSample),
                           alternative = "less")$p.value < p.value.threshold)
        if(T1 && T2 && T3 && T4 && T5 && T6 && T7 && T8) break;
    }
    cat(paste("\n done tries", tries, "\n"))
    expect_true(T1 && T2 && T3 && T4 && T5 && T6 && T7 && T8)

})
date()


date()
test_that(" Init with different mutators", {
    max.tries <- 4
    for(tries in 1:max.tries) {
        T1 <- T2 <- T3 <- T4 <- T5 <- T6 <- T7 <- T8 <- TRUE
        cat("\n z2-ossl: a runif is", runif(1), "\n")
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
                                    sampleEvery = 0.01, thresholdWhole = 1e-20,  
                                    seed = NULL, onlyCancer = FALSE, detectionProb = NA)
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
                                    sampleEvery = 0.01, thresholdWhole = 1e-20,
                                    seed = NULL, onlyCancer = FALSE, detectionProb = NA)
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
                                    sampleEvery = 0.01, thresholdWhole = 1e-20,  
                                    seed = NULL, onlyCancer = FALSE, detectionProb = NA)
        T1 <- (wilcox.test(NClonesOSS(m1.pg1.b), NClonesOSS(m1.pg1.a),
                                alternative = "less")$p.value < p.value.threshold)
        T2 <- (wilcox.test(NClonesOSS(m1.pg1.c), NClonesOSS(m1.pg1.b),
                                alternative = "less")$p.value < p.value.threshold)
        T3 <- (t.test(rowSums(m1.pg1.a$popSample),rowSums(m1.pg1.b$popSample),
                           alternative = "greater")$p.value < p.value.threshold)
        T4 <- (t.test(rowSums(m1.pg1.b$popSample),rowSums(m1.pg1.c$popSample),
                           alternative = "greater")$p.value < p.value.threshold)
        if(T1 && T2 && T3 && T4 && T5 && T6 && T7 && T8) break;
    }
    cat(paste("\n done tries", tries, "\n"))
    expect_true(T1 && T2 && T3 && T4 && T5 && T6 && T7 && T8)

})
date()


cat(paste("\n Finished test.mutator-oncoSimulSample-long.R test at", date(), "\n"))

