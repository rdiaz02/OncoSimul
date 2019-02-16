## Repeat tests in test.mutator, using oncoSimulSample.
## This is a concession to extreme paranoia.


cat(paste("\n Starting test.mutator-oncoSimulSample.R test at", date(), "\n"))
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

## These next two tests are probably the strongest (provided we accept
## using the initMutant) as we compare observed with expected and the
## estimated effect of mutator


## Do it with pops small here, for speed, and then with many more in long.
## But we can still fail them just by chance. This is bad. Could use a
## loop and catch it and repeat.




date()
test_that("McFL: Mutator increases by given factor with per-gene-mut rates: major axis and chi-sq test", {
    ## Two cases: mutator and no mutator, with variable mutation rates.
    ## rates such that rates of no mutator = rates of mutator * mutator.
    ## Why not compare mutlitplication factor keeping mutation rates
    ## constant? Because specially with mutator and large diffs in mut
    ## rates, with oncoSimulSample you undersample variation with
    ## wholePop, etc.
    ## Setings similar to oss11 in per-gene-mutation-rates but with the mutator
    max.tries <- 4
    for(tries in 1:max.tries) {
    
    
    cat("\n MCFL: AEu8: a runif is", runif(1), "\n")
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
                           fe,  detectionProb = NA,
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
                           fe,  detectionProb = NA,
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
    ## Compare mutator with no mutator
    T4 <- (chisq.test(snomSampl("oreoisasabgene", m1.pg1.b),
                           snomSampl("oreoisasabgene", m1.pg2.b))$p.value > p.fail)
    y <- sqrt(snomSampl("oreoisasabgene", m1.pg1.b))
    x <- sqrt(snomSampl("oreoisasabgene", m1.pg2.b))
    mma <- smatr::ma(y ~ x, slope.test = 1, elev.test = 0) ## From smatr package, for major axis
    ## intercept not different from 0
    T5 <- (mma$elevtest[[1]]$p > p.fail)
    T6 <- (mma$slopetest[[1]]$p > p.fail)
    ## We could use a lm and do a simultaneous test on both slope and
    ## intercept as. But this is really asking for major axis regression
    ## lm1 <- lm(snomSampl("oreoisasabgene", m1.pg1.b) ~
    ##               snomSampl("oreoisasabgene", m1.pg2.b))
    ## ## test intercept is 0, slope is 1. Not technically fully correct, as
    ## ## X variable has noise. We should do major axis or similar and these
    ## ## are counts.
    ## expect_true(linearHypothesis(lm1, diag(2), c(0, 1))[["Pr(>F)"]][2] >
    ##             p.fail)
    if( T1 && T3 && T4 && T5 && T6) break;
    }
    cat(paste("\n done tries", tries, "\n"))
    expect_true(T1 && T3 && T4 && T5 && T6)
})
date()






date()
test_that("MCFL Relative ordering of number of clones with init mutant of mutator effects", {
    max.tries <- 4
    for(tries in 1:max.tries) {
        T1 <- T2 <- T3 <- T4 <- T5 <- T6 <- T7 <- T8 <- TRUE
        ## Here we stop on  popSize after short model. All have same small s.
        cat("\n mcx2bc: a runif is", runif(1), "\n")
        pops <- 50
        ni <- rep(0.01, 50)
        names(ni) <- c("a", "b", "c", "d", paste0("n", 1:46))
        fe <- allFitnessEffects(noIntGenes = ni)
        fm6 <- allMutatorEffects(noIntGenes = c("a" = .05,
                                                "b" = 1,
                                                "c" = 10,
                                                "d" = 50))
        nca <- oncoSimulSample(pops, fe,  detectionProb = NA, muEF = fm6, finalTime =250,
                               mutationPropGrowth = FALSE,
                               initSize = 1e4,
                               initMutant = "a", detectionSize = 10200,
                               sampleEvery = 0.01, thresholdWhole = 1e-20,
                               onlyCancer = FALSE, model = "McFL")
        ncb <- oncoSimulSample(pops, fe,  detectionProb = NA, muEF = fm6, finalTime =250,
                               mutationPropGrowth = FALSE,
                               initSize = 1e4,
                               initMutant = "b", detectionSize = 10200,
                               sampleEvery = 0.01, thresholdWhole = 1e-20,
                               onlyCancer = FALSE, model = "McFL")
        ncc <- oncoSimulSample(pops, fe,  detectionProb = NA, muEF = fm6, finalTime =250,
                               mutationPropGrowth = FALSE,
                               initSize = 1e4,
                               initMutant = "c",detectionSize = 10200,
                               sampleEvery = 0.01, thresholdWhole = 1e-20,
                               onlyCancer = FALSE, model = "McFL")
        ncd <- oncoSimulSample(pops, fe,  detectionProb = NA, muEF = fm6, finalTime =250,
                               mutationPropGrowth = FALSE,
                               initSize = 1e4,
                               initMutant = "d",detectionSize = 10200,
                               sampleEvery = 0.01, thresholdWhole = 1e-20,
                               onlyCancer = FALSE, model = "McFL")
        T4 <- ( wilcox.test(nca$popSummary[, "NumClones"],
                                 ncb$popSummary[, "NumClones"],
                                 alternative = "less")$p.value < p.value.threshold)
        T5 <- (wilcox.test(ncb$popSummary[, "NumClones"],
                                ncc$popSummary[, "NumClones"],
                                alternative = "less")$p.value < p.value.threshold)
        T6 <- ( wilcox.test(ncc$popSummary[, "NumClones"],
                                 ncd$popSummary[, "NumClones"],
                                 alternative = "less")$p.value < p.value.threshold)
        nca$popSummary[, c(1:3, 8:9)]
        ncb$popSummary[, c(1:3, 8:9)]
        ncc$popSummary[, c(1:3, 8:9)]
        ncd$popSummary[, c(1:3, 8:9)]
        T1 <- (t.test(rowSums(nca$popSample), rowSums(ncb$popSample),
                           alternative = "less")$p.value < p.value.threshold)
        T2 <- (t.test(rowSums(ncb$popSample), rowSums(ncc$popSample),
                           alternative = "less")$p.value < p.value.threshold)
        T3 <- (t.test(rowSums(ncc$popSample), rowSums(ncd$popSample),
                           alternative = "less")$p.value < p.value.threshold)
        if(T1 && T2 && T3 && T4 && T5 && T6 && T7 && T8) break;
    }
    cat(paste("\n done tries", tries, "\n"))
    expect_true(T1 && T2 && T3 && T4 && T5 && T6 && T7 && T8)
})
date()

date()
test_that("Relative ordering of number of clones with mut prop growth and init and scrambled names", {
    max.tries <- 4
    for(tries in 1:max.tries) {
        T1 <- T2 <- T3 <- T4 <- T5 <- T6 <- T7 <- T8 <- TRUE
        ## Stopping on time; s > 0 , but all have same growth rate.
        cat("\n x2ef: a runif is", runif(1), "\n")
        pops <- 10
        ft <- 1
        lni <- 200
        no <- 5e3
        ni <- c(5, 0, rep(0, lni))
        ## scramble around names
        names(ni) <- c("thisistheagene",
                       "thisisthebgene",
                       replicate(lni,
                                 paste(sample(letters, 12), collapse = "")))
        ni <- ni[order(names(ni))]
        fe <- allFitnessEffects(noIntGenes = ni)
        fm1 <- allMutatorEffects(noIntGenes = c("thisistheagene" = 5))
        mpg <- oncoSimulSample(pops, fe,  detectionProb = NA, muEF = fm1,
                               finalTime = ft,
                               mutationPropGrowth = TRUE,
                               initSize = no,
                               sampleEvery = 0.01, thresholdWhole = 1e-20,
                               initMutant = "thisistheagene",
                               detectionSize = 1e9,
                               detectionDrivers = 9999, seed = NULL,
                               onlyCancer = FALSE)
        mnpg <- oncoSimulSample(pops, fe,  detectionProb = NA, muEF = fm1,
                                finalTime = ft,
                                mutationPropGrowth = FALSE,
                                initSize = no,
                                sampleEvery = 0.01, thresholdWhole = 1e-20,
                                initMutant = "thisistheagene",
                                detectionSize = 1e9,
                                detectionDrivers = 9999, seed = NULL,
                                onlyCancer = FALSE)
        pg <- oncoSimulSample(pops, fe,  detectionProb = NA, 
                              finalTime = ft,
                              mutationPropGrowth = TRUE,
                              initSize = no,
                              sampleEvery = 0.01, thresholdWhole = 1e-20,
                              initMutant = "thisistheagene",
                              detectionSize = 1e9,
                              detectionDrivers = 9999, seed = NULL,
                              onlyCancer = FALSE)
        npg <- oncoSimulSample(pops, fe,  detectionProb = NA, 
                               finalTime = ft,
                               mutationPropGrowth = FALSE,
                               initSize = no,
                               sampleEvery = 0.01, thresholdWhole = 1e-20,
                               initMutant = "thisistheagene",
                               detectionSize = 1e9,
                               detectionDrivers = 9999, seed = NULL,
                               onlyCancer = FALSE)
        ## These are the real tests
        T1 <- ( wilcox.test(mpg$popSummary[, "NumClones"], mnpg$popSummary[, "NumClones"],
                                 alternative = "greater")$p.value < p.value.threshold)
        T2 <- (wilcox.test(mpg$popSummary[, "NumClones"], pg$popSummary[, "NumClones"],
                                alternative = "greater")$p.value < p.value.threshold)
        T3 <-  (wilcox.test(mnpg$popSummary[, "NumClones"], npg$popSummary[, "NumClones"],
                                 alternative = "greater")$p.value < p.value.threshold)
        T4 <-  (wilcox.test(pg$popSummary[, "NumClones"], npg$popSummary[, "NumClones"],
                                 alternative = "greater")$p.value < p.value.threshold)
        T5 <- (t.test(rowSums(mpg$popSample),rowSums(mnpg$popSample),
                           alternative = "greater")$p.value < p.value.threshold)
        T6 <- (t.test(rowSums(mpg$popSample),rowSums(pg$popSample),
                           alternative = "greater")$p.value < p.value.threshold)
        T7 <- (t.test(rowSums(mnpg$popSample),rowSums(npg$popSample),
                           alternative = "greater")$p.value < p.value.threshold)
        T8 <- (t.test(rowSums(pg$popSample),rowSums(npg$popSample),
                           alternative = "greater")$p.value < p.value.threshold)
        ## mpg$popSummary[, c(1:3, 8:9)]
        ## mnpg$popSummary[, c(1:3, 8:9)]
        ## pg$popSummary[, c(1:3, 8:9)]
        ## npg$popSummary[, c(1:3, 8:9)]
        if(T1 && T2 && T3 && T4 && T5 && T6 && T7 && T8) break;
    }
    cat(paste("\n done tries", tries, "\n"))
    expect_true(T1 && T2 && T3 && T4 && T5 && T6 && T7 && T8)
})
date()


date()
test_that("McFL: Relative ordering of number of clones with mut prop growth and init and scrambled names", {
    max.tries <- 4
    for(tries in 1:max.tries) {
        T1 <- T2 <- T3 <- T4 <- T5 <- T6 <- T7 <- T8 <- TRUE
        ## Stopping on time; s > 0 but all same growth rate.
        cat("\n x2gh: a runif is", runif(1), "\n")
        pops <- 15
        ft <- 1
        lni <- 200
        no <- 1e3
        ni <- c(5, 0, rep(0, lni))
        ## scramble around names
        names(ni) <- c("thisistheagene",
                       "thisisthebgene",
                       replicate(lni,
                                 paste(sample(letters, 12), collapse = "")))
        ni <- ni[order(names(ni))]
        fe <- allFitnessEffects(noIntGenes = ni)
        fm1 <- allMutatorEffects(noIntGenes = c("thisistheagene" = 5))
        mpg <- oncoSimulSample(pops, fe,  detectionProb = NA, muEF = fm1,
                               finalTime = ft,
                               mutationPropGrowth = TRUE,
                               initSize = no, model = "McFL",
                               initMutant = "thisistheagene",
                               sampleEvery = 0.01, thresholdWhole = 1e-20,
                               detectionSize = 1e9,
                               detectionDrivers = 9999, seed = NULL,
                               onlyCancer = FALSE)
        mnpg <- oncoSimulSample(pops, fe,  detectionProb = NA, muEF = fm1,
                                finalTime = ft,
                                mutationPropGrowth = FALSE,
                                initSize = no, model = "McFL",
                                initMutant = "thisistheagene",
                                sampleEvery = 0.01, thresholdWhole = 1e-20,
                                detectionSize = 1e9,
                                detectionDrivers = 9999, seed = NULL,
                                onlyCancer = FALSE)
        pg <- oncoSimulSample(pops, fe,  detectionProb = NA, 
                              finalTime = ft,
                              mutationPropGrowth = TRUE,
                              initSize = no, model = "McFL",
                              initMutant = "thisistheagene",
                              sampleEvery = 0.01, thresholdWhole = 1e-20,
                              detectionSize = 1e9,
                              detectionDrivers = 9999, seed = NULL,
                              onlyCancer = FALSE)
        npg <- oncoSimulSample(pops, fe,  detectionProb = NA, 
                               finalTime = ft,
                               mutationPropGrowth = FALSE,
                               initSize = no, model = "McFL",
                               initMutant = "thisistheagene",
                               sampleEvery = 0.01, thresholdWhole = 1e-20,
                               detectionSize = 1e9,
                               detectionDrivers = 9999, seed = NULL,
                               onlyCancer = FALSE)
        ## These are the real tests
        T1 <- ( wilcox.test(mpg$popSummary[, "NumClones"], mnpg$popSummary[, "NumClones"],
                                 alternative = "greater")$p.value < p.value.threshold)
        T2 <- (wilcox.test(mpg$popSummary[, "NumClones"], pg$popSummary[, "NumClones"],
                                alternative = "greater")$p.value < p.value.threshold)
        T3 <- ( wilcox.test(mnpg$popSummary[, "NumClones"], npg$popSummary[, "NumClones"],
                                 alternative = "greater")$p.value < p.value.threshold)
        T4 <- ( wilcox.test(pg$popSummary[, "NumClones"], npg$popSummary[, "NumClones"],
                                 alternative = "greater")$p.value < p.value.threshold)
        T5 <- (t.test(rowSums(mpg$popSample),rowSums(mnpg$popSample),
                           alternative = "greater")$p.value < p.value.threshold)
        T6 <- (t.test(rowSums(mpg$popSample),rowSums(pg$popSample),
                           alternative = "greater")$p.value < p.value.threshold)
        T7 <- (t.test(rowSums(mnpg$popSample),rowSums(npg$popSample),
                           alternative = "greater")$p.value < p.value.threshold)
        T8 <- (t.test(rowSums(pg$popSample),rowSums(npg$popSample),
                           alternative = "greater")$p.value < p.value.threshold)
        if(T1 && T2 && T3 && T4 && T5 && T6 && T7 && T8) break;
    }
    cat(paste("\n done tries", tries, "\n"))
    expect_true(T1 && T2 && T3 && T4 && T5 && T6 && T7 && T8)
})
date()



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




date()
test_that("McFL: Same mu vector, different mutator; diffs in number muts, tiny t", {
    max.tries <- 4
    for(tries in 1:max.tries) {
        T1 <- T2 <- T3 <- T4 <- T5 <- T6 <- T7 <- T8 <- TRUE
        ## Here, there is no reproduction or death. Just mutation. And no double
        ## mutants either.
        ## We test:
        ##  - mutator increases mutation rates as seen in:
        ##        - number of clones created
        ##        - number of total mutation events
        cat("\n nm2: a runif is", runif(1), "\n")
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
                                 fe,  detectionProb = NA,
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
                                 seed = NULL, onlyCancer = FALSE)
        pop100 <- oncoSimulSample(pops,
                                  fe,  detectionProb = NA,
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
                                  seed = NULL, onlyCancer = FALSE)
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
test_that(" MCFL Init with different mutators", {
    max.tries <- 4
    for(tries in 1:max.tries) {
        T1 <- T2 <- T3 <- T4 <- T5 <- T6 <- T7 <- T8 <- TRUE
        cat("\n mcz2: a runif is", runif(1), "\n")
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
                                    fe,  detectionProb = NA,
                                    mu = pg1,
                                    muEF = m1,
                                    finalTime = ft,
                                    mutationPropGrowth = FALSE,
                                    initSize = no,
                                    model = "McFL",
                                    initMutant = "hereisoneagene",
                                    detectionSize = 1e9,
                                    detectionDrivers = 9999,
                                    sampleEvery = 0.01, thresholdWhole = 1e-20, 
                                    seed = NULL, onlyCancer = FALSE)
        m1.pg1.b <- oncoSimulSample(pops,
                                    fe,  detectionProb = NA,
                                    mu = pg1,
                                    muEF = m1,
                                    finalTime = ft,
                                    mutationPropGrowth = FALSE,
                                    initSize = no,
                                    model = "McFL",                             
                                    initMutant = "oreoisasabgene",
                                    detectionSize = 1e9,
                                    detectionDrivers = 9999,
                                    sampleEvery = 0.01, thresholdWhole = 1e-20, 
                                    seed = NULL, onlyCancer = FALSE)
        m1.pg1.c <- oncoSimulSample(pops,
                                    fe,  detectionProb = NA,
                                    mu = pg1,
                                    muEF = m1,
                                    finalTime = ft,
                                    mutationPropGrowth = FALSE,
                                    initSize = no,
                                    model = "McFL",                           
                                    initMutant = "nnhsisthecgene",
                                    detectionSize = 1e9,
                                    detectionDrivers = 9999,
                                    sampleEvery = 0.01, thresholdWhole = 1e-20,
                                    seed = NULL, onlyCancer = FALSE)
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

cat(paste("\n Finished test.mutator-oncoSimulSample.R test at", date(), "\n"))



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
