## This test takes about 12 to 13 minutes in my laptop

cat(paste("\n Starting test.mutator-long.R test at", date()))
cat(paste("\n            a runif ", runif(1), "\n"))
## RNGkind("L'Ecuyer-CMRG") ## for the mclapplies
date()


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

## Minifunctions for testing
medianNClones <- function(x) {
    median(summary(x)$NumClones)
}

NClones <- function(x) {
    summary(x)$NumClones
}


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

snom <- function(name, out) {
    ## observed without the init
    cs <- colSums(OncoSimulR:::geneCounts(out))
    ii <- which(names(cs) == name)
    cs <- cs[-ii]
    cs[order(names(cs))]
}

sm <- function(name, out) {
    ## totals for a given gene
    cs <- colSums(OncoSimulR:::geneCounts(out))
    ii <- which(names(cs) == name)
    cs[ii]
}

totalind <- function(out) {
    ## total num indivs
  sum(unlist(lapply(out, function(x) x$TotalPopSize)))  
}

smA <- function(out) {
    ## totals counts for all. So total mutated over all.
    cs <- colSums(OncoSimulR:::geneCounts(out))
    sum(cs)
}

smAPi <- function(out) {
    ## smAnom but divided by number of individuals.
    smA(out)/totalind(out)
}

smAnom <- function(out, name) {
    ## totals counts for all. So total mutated over all.
    ## but remove one, the fixed starting one.
    cs <- colSums(OncoSimulR:::geneCounts(out))
    ii <- which(names(cs) == name)
    cs <- cs[-ii]
    sum(cs)
}

smAnomPi <- function(out, name) {
    ## smAnom but divided by number of individuals.
    smAnom(out, name)/totalind(out)
}

p.value.threshold <- 0.001

cat(paste("\n            a second runif ", runif(1), "\n"))


date()
test_that("MCFL Relative ordering of number of clones with init mutant of mutator effects", {
    max.tries <- 4
    for(tries in 1:max.tries) {
        T1 <- T2 <- T3 <- T4 <- T5 <- T6 <- T7 <- T8 <- TRUE
        ## Here we stop on finalTime, not popSize
        ## Can occasionally blow up with pE.f: pE not finite.
        cat("\n mcx2bc: a runif is", runif(1), "\n")
        pops <- 50
        ni <- rep(0.01, 50)
        names(ni) <- c("a", "b", "c", "d", paste0("n", 1:46))
        fe <- allFitnessEffects(noIntGenes = ni)
        fm6 <- allMutatorEffects(noIntGenes = c("a" = .05,
                                                "b" = 1,
                                                "c" = 10,
                                                "d" = 50))
        nca <- oncoSimulPop(pops, fe, muEF = fm6, finalTime =50,
                            mutationPropGrowth = FALSE,
                            initSize = 1e4,
                            initMutant = "a",
                            sampleEvery = 0.01,
                            keepEvery = 5,
                            onlyCancer = FALSE, detectionProb = NA, seed = NULL, mc.cores = 2, model = "McFL")
        ncb <- oncoSimulPop(pops, fe, muEF = fm6, finalTime =50,
                            mutationPropGrowth = FALSE,
                            initSize = 1e4,
                            initMutant = "b",
                            sampleEvery = 0.01,
                            keepEvery = 5,
                            onlyCancer = FALSE, detectionProb = NA, seed = NULL, mc.cores = 2, model = "McFL")
        ncc <- oncoSimulPop(pops, fe, muEF = fm6, finalTime =50,
                            mutationPropGrowth = FALSE,
                            initSize = 1e4,
                            initMutant = "c",
                            sampleEvery = 0.01,
                            keepEvery = 5,
                            onlyCancer = FALSE, detectionProb = NA, seed = NULL, mc.cores = 2, model = "McFL")
        ncd <- oncoSimulPop(pops, fe, muEF = fm6, finalTime =50,
                            mutationPropGrowth = FALSE,
                            initSize = 1e4,
                            initMutant = "d",
                            sampleEvery = 0.01,
                            keepEvery = 5,
                            onlyCancer = FALSE, detectionProb = NA, seed = NULL, mc.cores = 2, model = "McFL")
        T1 <- ( wilcox.test(summary(nca)$NumClones, 
                                 summary(ncb)$NumClones, alternative = "less")$p.value < p.value.threshold)
        T2 <- (wilcox.test(summary(ncb)$NumClones ,
                                summary(ncc)$NumClones, alternative = "less")$p.value < p.value.threshold)
        T3 <- ( wilcox.test(summary(ncc)$NumClones,
                                 summary(ncd)$NumClones, alternative = "less")$p.value < p.value.threshold)
        summary(nca)[, c(1:3, 8:9)]
        summary(ncb)[, c(1:3, 8:9)]
        summary(ncc)[, c(1:3, 8:9)]
        summary(ncd)[, c(1:3, 8:9)]
        T4 <- (t.test(mutsPerClone(nca), mutsPerClone(ncb), alternative = "less")$p.value < p.value.threshold)
        T5 <- (t.test(mutsPerClone(ncb), mutsPerClone(ncc), alternative = "less")$p.value < p.value.threshold)
        T6 <- (t.test(mutsPerClone(ncc), mutsPerClone(ncd), alternative = "less")$p.value < p.value.threshold)
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
        cat("\n z2: a runif is", runif(1), "\n")
        pops <- 60
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
        m1.pg1.a <- oncoSimulPop(pops,
                                 fe,
                                 mu = pg1,
                                 muEF = m1,
                                 finalTime = ft,
                                 mutationPropGrowth = FALSE,
                                 initSize = no,
                                 initMutant = "hereisoneagene",
                                 sampleEvery = 0.01, keepEvery = 5, seed = NULL,
                                 onlyCancer = FALSE, detectionProb = NA, mc.cores = 2)
        m1.pg1.b <- oncoSimulPop(pops,
                                 fe,
                                 mu = pg1,
                                 muEF = m1,
                                 finalTime = ft,
                                 mutationPropGrowth = FALSE,
                                 initSize = no,
                                 initMutant = "oreoisasabgene",
                                 sampleEvery = 0.01, keepEvery = 5, seed = NULL,
                                 onlyCancer = FALSE, detectionProb = NA, mc.cores = 2)
        m1.pg1.c <- oncoSimulPop(pops,
                                 fe,
                                 mu = pg1,
                                 muEF = m1,
                                 finalTime = ft,
                                 mutationPropGrowth = FALSE,
                                 initSize = no,
                                 initMutant = "nnhsisthecgene",
                                 sampleEvery = 0.01, keepEvery = 5, seed = NULL,
                                 onlyCancer = FALSE, detectionProb = NA, mc.cores = 2)
        T1 <- (wilcox.test(NClones(m1.pg1.a),  NClones(m1.pg1.b), alternative = "greater")$p.value < p.value.threshold)
        T2 <- (wilcox.test(NClones(m1.pg1.b),  NClones(m1.pg1.c), alternative = "greater")$p.value < p.value.threshold)
        T3 <- (t.test(mutsPerClone(m1.pg1.a), mutsPerClone(m1.pg1.b), alternative = "greater")$p.value < p.value.threshold)
        T4 <- (t.test(mutsPerClone(m1.pg1.b), mutsPerClone(m1.pg1.c), alternative = "greater")$p.value < p.value.threshold)
        T5 <- (smAnomPi(m1.pg1.a, "hereisoneagene") >
                    smAnomPi(m1.pg1.b, "oreoisasabgene"))
        T6 <- (smAnomPi(m1.pg1.b, "oreoisasabgene") >
                    smAnomPi(m1.pg1.c, "nnhsisthecgene"))
        summary(m1.pg1.a)[, c(1:3, 8:9)]
        summary(m1.pg1.b)[, c(1:3, 8:9)]
        summary(m1.pg1.c)[, c(1:3, 8:9)]
        if(T1 && T2 && T3 && T4 && T5 && T6 && T7 && T8) break;
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
        cat("\n nm0: a runif is", runif(1), "\n")
        pops <- 40
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
        pop10 <- oncoSimulPop(pops,
                              fe,
                              mu = muvector,
                              muEF = m10,
                              finalTime = ft,
                              mutationPropGrowth = FALSE,
                              initSize = no,
                              initMutant = names(mutator10),
                              onlyCancer = FALSE, detectionProb = NA, seed = NULL, mc.cores = 2)
        pop100 <- oncoSimulPop(pops,
                               fe,
                               mu = muvector,
                               muEF = m100,
                               finalTime = ft,
                               mutationPropGrowth = FALSE,
                               initSize = no,
                               initMutant = names(mutator10),
                               onlyCancer = FALSE, detectionProb = NA, seed = NULL, mc.cores = 2)
        ## number of total mutations
        T1 <- (smAnomPi(pop10, names(mutator10)) < smAnomPi(pop100, names(mutator100)))
        ## number of clones
        T2 <- (wilcox.test(NClones(pop100),  NClones(pop10), alternative = "greater")$p.value < p.value.threshold)
        T3 <- (t.test(mutsPerClone(pop100), mutsPerClone(pop10), alternative = "greater")$p.value < p.value.threshold)
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
        cat("\n nm1: a runif is", runif(1), "\n")
        pops <- 40
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
        pop10 <- oncoSimulPop(pops,
                              fe,
                              mu = muvector,
                              muEF = m10,
                              finalTime = ft,
                              mutationPropGrowth = FALSE,
                              initSize = no,
                              initMutant = names(mutator10),
                              onlyCancer = FALSE, detectionProb = NA, seed = NULL, mc.cores = 2)
        pop100 <- oncoSimulPop(pops,
                               fe,
                               mu = muvector,
                               muEF = m100,
                               finalTime = ft,
                               mutationPropGrowth = FALSE,
                               initSize = no,
                               initMutant = names(mutator10),
                               onlyCancer = FALSE, detectionProb = NA, seed = NULL, mc.cores = 2)
        ## number of total mutations
        T1 <- (smAnomPi(pop10, names(mutator10)) < smAnomPi(pop100, names(mutator100)))
        ## number of clones
        T2 <- (wilcox.test(NClones(pop100),  NClones(pop10), alternative = "greater")$p.value < p.value.threshold)
        T3 <- (t.test(mutsPerClone(pop100), mutsPerClone(pop10), alternative = "greater")$p.value < p.value.threshold)
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
        ## Can occasionally blow up with pE.f: pE not finite.
        cat("\n x2ef: a runif is", runif(1), "\n")
        pops <- 90
        ft <- 1
        lni <- 200
        no <- 5e3
        ni <- c(5, 2, rep(0, lni))
        ## scramble around names
        names(ni) <- c("thisistheagene",
                       "thisisthebgene",
                       replicate(lni,
                                 paste(sample(letters, 12), collapse = "")))
        ni <- ni[order(names(ni))]
        fe <- allFitnessEffects(noIntGenes = ni)
        fm1 <- allMutatorEffects(noIntGenes = c("thisistheagene" = 5))
        mpg <- oncoSimulPop(pops, fe, muEF = fm1,
                            finalTime = ft,
                            mutationPropGrowth = TRUE,
                            initSize = no, sampleEvery = 0.01,
                            initMutant = "thisistheagene",
                            onlyCancer = FALSE, detectionProb = NA, seed = NULL, mc.cores = 2)
        mnpg <- oncoSimulPop(pops, fe, muEF = fm1,
                             finalTime = ft,
                             mutationPropGrowth = FALSE,
                             initSize = no,sampleEvery = 0.01,
                             initMutant = "thisistheagene",
                             onlyCancer = FALSE, detectionProb = NA, seed = NULL, mc.cores = 2)
        pg <- oncoSimulPop(pops, fe, 
                           finalTime = ft,
                           mutationPropGrowth = TRUE,
                           initSize = no, sampleEvery = 0.01,
                           initMutant = "thisistheagene",
                           onlyCancer = FALSE, detectionProb = NA, seed = NULL, mc.cores = 2)
        npg <- oncoSimulPop(pops, fe, 
                            finalTime = ft,
                            mutationPropGrowth = FALSE,
                            initSize = no, sampleEvery = 0.01,
                            initMutant = "thisistheagene",
                            onlyCancer = FALSE, detectionProb = NA, seed = NULL, mc.cores = 2)
        ## These are the real tests
        T1 <- ( wilcox.test(summary(mpg)$NumClones, 
                                 summary(mnpg)$NumClones, alternative = "greater")$p.value < p.value.threshold)
        T2 <- (wilcox.test(summary(mpg)$NumClones, 
                                summary(pg)$NumClones, alternative = "greater")$p.value < p.value.threshold)
        T3 <- ( wilcox.test(summary(mnpg)$NumClones, 
                                 summary(npg)$NumClones, alternative = "greater")$p.value < p.value.threshold)
        T4 <- ( wilcox.test(summary(pg)$NumClones, 
                                 summary(npg)$NumClones, alternative = "greater")$p.value < p.value.threshold)
        T5 <- (t.test(mutsPerClone(mpg), mutsPerClone(mnpg), alternative = "greater")$p.value < p.value.threshold)
        T6 <- (t.test(mutsPerClone(mpg), mutsPerClone(pg), alternative = "greater")$p.value < p.value.threshold)
        T7 <- (t.test(mutsPerClone(mnpg), mutsPerClone(npg), alternative = "greater")$p.value < p.value.threshold)
        T8 <- (t.test(mutsPerClone(pg), mutsPerClone(npg), alternative = "greater")$p.value < p.value.threshold)
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
        cat("\n u8_1: a runif is", runif(1), "\n")
        pops <- 120
        ft <- .0001
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
        mutator1 <- rep(1, lni + 3)
        pg1 <- runif(lni + 3, min = 1e-9, max = 1e-5) ## max should not be
        ## huge here as mutator
        ## is 34. Can get beyond
        ## 1
        names(mutator1) <- sample(names(ni))
        names(pg1) <- sample(names(ni))
        mutator1["oreoisasabgene"] <- 200
        m1 <- allMutatorEffects(noIntGenes = mutator1)
        pg1["hereisoneagene"] <- 1e-3 ## 1e-3
        m1.pg1.b <- oncoSimulPop(pops,
                                 fe,
                                 mu = pg1,
                                 muEF = m1,
                                 finalTime = ft,
                                 mutationPropGrowth = FALSE,
                                 initSize = no,
                                 initMutant ="oreoisasabgene",
                                 onlyCancer = FALSE, detectionProb = NA, seed = NULL, mc.cores = 2)
        expect_true(sm("oreoisasabgene", m1.pg1.b) == totalind(m1.pg1.b))
        ## If you want to see the numbers
        ## enom("oreoisasabgene", pg1)
        ## snom("oreoisasabgene", m1.pg1.b)
        p.fail <- 1e-3
        T1 <- (chisq.test(snom("oreoisasabgene", m1.pg1.b),
                               p = pnom("oreoisasabgene", pg1, no, pops))$p.value > p.fail)
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
        ## Can occasionally blow up with pE.f: pE not finite.
        cat("\n x2cd: a runif is", runif(1), "\n")
        pops <- 50
        ni <- rep(0, 50)
        names(ni) <- c("a", "b", "c", "d", paste0("n", 1:46))
        fe <- allFitnessEffects(noIntGenes = ni)
        fm6 <- allMutatorEffects(noIntGenes = c("a" = .05,
                                                "b" = 1,
                                                "c" = 10,
                                                "d" = 50))
        nca <- oncoSimulPop(pops, fe, muEF = fm6, finalTime =50,
                            mutationPropGrowth = FALSE,
                            initSize = 1e4,
                            initMutant = "a",
                            sampleEvery = 0.01,
                            keepEvery = 5,
                            onlyCancer = FALSE, detectionProb = NA, seed = NULL, mc.cores = 2)
        ncb <- oncoSimulPop(pops, fe, muEF = fm6, finalTime =50,
                            mutationPropGrowth = FALSE,
                            initSize = 1e4,
                            initMutant = "b",
                            sampleEvery = 0.01,
                            keepEvery = 5,
                            onlyCancer = FALSE, detectionProb = NA, seed = NULL, mc.cores = 2)
        ncc <- oncoSimulPop(pops, fe, muEF = fm6, finalTime =50,
                            mutationPropGrowth = FALSE,
                            initSize = 1e4,
                            initMutant = "c",
                            sampleEvery = 0.01,
                            keepEvery = 5,
                            onlyCancer = FALSE, detectionProb = NA, seed = NULL, mc.cores = 2)
        ncd <- oncoSimulPop(pops, fe, muEF = fm6, finalTime =50,
                            mutationPropGrowth = FALSE,
                            initSize = 1e4,
                            initMutant = "d",
                            sampleEvery = 0.01,
                            keepEvery = 5,
                            onlyCancer = FALSE, detectionProb = NA, seed = NULL, mc.cores = 2)
        ## I once saw a weird thing
        expect_true(var(summary(nca)$NumClones) > 1e-4)
        expect_true(var(summary(ncb)$NumClones) > 1e-4)
        expect_true(var(summary(ncc)$NumClones) > 1e-4)
        expect_true(var(summary(ncd)$NumClones) > 1e-4)
        ## These are the real tests
        T1 <- ( wilcox.test(summary(nca)$NumClones,
                                 summary(ncb)$NumClones, alternative = "less")$p.value < p.value.threshold)
        T2 <- (wilcox.test(summary(ncb)$NumClones,
                                summary(ncc)$NumClones, alternative = "less")$p.value < p.value.threshold)
        T3 <- ( wilcox.test(summary(ncc)$NumClones,
                                 summary(ncd)$NumClones, alternative = "less")$p.value < p.value.threshold)
        T4 <- (t.test(mutsPerClone(nca), mutsPerClone(ncb), alternative = "less")$p.value < p.value.threshold)
        T5 <- (t.test(mutsPerClone(ncb), mutsPerClone(ncc), alternative = "less")$p.value < p.value.threshold)
        T6 <- (t.test(mutsPerClone(ncc), mutsPerClone(ncd), alternative = "less")$p.value < p.value.threshold)
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
        ## Here we stop on finalTime, not popSize
        ## Can occasionally blow up with pE.f: pE not finite.
        cat("\n x2bc: a runif is", runif(1), "\n")
        pops <- 60
        ni <- rep(0.01, 50)
        names(ni) <- c("a", "b", "c", "d", paste0("n", 1:46))
        fe <- allFitnessEffects(noIntGenes = ni)
        fm6 <- allMutatorEffects(noIntGenes = c("a" = .05,
                                                "b" = 1,
                                                "c" = 10,
                                                "d" = 50))
        nca <- oncoSimulPop(pops, fe, muEF = fm6, finalTime =50,
                            mutationPropGrowth = FALSE,
                            initSize = 1e4,
                            initMutant = "a",
                            sampleEvery = 0.01,
                            keepEvery = 5,
                            onlyCancer = FALSE, detectionProb = NA, seed = NULL, mc.cores = 2)
        ncb <- oncoSimulPop(pops, fe, muEF = fm6, finalTime =50,
                            mutationPropGrowth = FALSE,
                            initSize = 1e4,
                            initMutant = "b",
                            sampleEvery = 0.01,
                            keepEvery = 5,
                            onlyCancer = FALSE, detectionProb = NA, seed = NULL, mc.cores = 2)
        ncc <- oncoSimulPop(pops, fe, muEF = fm6, finalTime =50,
                            mutationPropGrowth = FALSE,
                            initSize = 1e4,
                            initMutant = "c",
                            sampleEvery = 0.01,
                            keepEvery = 5,
                            onlyCancer = FALSE, detectionProb = NA, seed = NULL, mc.cores = 2)
        ncd <- oncoSimulPop(pops, fe, muEF = fm6, finalTime =50,
                            mutationPropGrowth = FALSE,
                            initSize = 1e4,
                            initMutant = "d",
                            sampleEvery = 0.01,
                            keepEvery = 5,
                            onlyCancer = FALSE, detectionProb = NA, seed = NULL, mc.cores = 2)
        T1 <- ( wilcox.test(summary(nca)$NumClones,
                                 summary(ncb)$NumClones, alternative = "less")$p.value < p.value.threshold)
        T2 <- (wilcox.test(summary(ncb)$NumClones,
                                summary(ncc)$NumClones, alternative = "less")$p.value < p.value.threshold)
        T3 <- ( wilcox.test(summary(ncc)$NumClones,
                                 summary(ncd)$NumClones, alternative = "less")$p.value < p.value.threshold)
        T4 <- (t.test(mutsPerClone(nca), mutsPerClone(ncb), alternative = "less")$p.value < p.value.threshold)
        T5 <- (t.test(mutsPerClone(ncb), mutsPerClone(ncc), alternative = "less")$p.value < p.value.threshold)
        T6 <- (t.test(mutsPerClone(ncc), mutsPerClone(ncd), alternative = "less")$p.value < p.value.threshold)
        summary(nca)[, c(1:3, 8:9)]
        summary(ncb)[, c(1:3, 8:9)]
        summary(ncc)[, c(1:3, 8:9)]
        summary(ncd)[, c(1:3, 8:9)]
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
        ## Can occasionally blow up with pE.f: pE not finite.
        cat("\n mcx1: a runif is", runif(1), "\n")
        pops <- 60
        fe <- allFitnessEffects(noIntGenes = c("a" = 0.12,
                                               "b" = 0.14,
                                               "c" = 0.16,
                                               "d" = 0.11))
        fm6 <- allMutatorEffects(noIntGenes = c("a" = 30,
                                                "b" = 30,
                                                "c" = 30,
                                                "d" = 30))
        nc1 <- oncoSimulPop(pops, fe, muEF = fm6, finalTime =100,
                            mutationPropGrowth = FALSE,
                            sampleEvery = 0.001,
                            keepEvery = 5,
                            initSize = 1e6, mc.cores = 2, model = "McFL",
                            onlyCancer = FALSE, detectionProb = NA, seed = NULL)
        fm8 <- allMutatorEffects(noIntGenes = c("a" = 1,
                                                "b" = 1,
                                                "c" = 1,
                                                "d" = 1))
        nc2 <- oncoSimulPop(pops, fe, muEF = fm8, finalTime =100,
                            mutationPropGrowth = FALSE,
                            sampleEvery = 0.001,
                            keepEvery = 5,
                            initSize = 1e6, mc.cores = 2, model = "McFL",
                            onlyCancer = FALSE, detectionProb = NA, seed = NULL)
        fm7 <- allMutatorEffects(noIntGenes = c("a" = 1e-3,
                                                "b" = 1e-3,
                                                "c" = 1e-3,
                                                "d" = 1e-3))
        nc3 <- oncoSimulPop(pops, fe, muEF = fm7, finalTime =100,
                            mutationPropGrowth = FALSE,
                            sampleEvery = 0.001,
                            keepEvery = 5,
                            initSize = 1e6, mc.cores = 2, model = "McFL",
                            onlyCancer = FALSE, detectionProb = NA, seed = NULL)
        T1 <- (wilcox.test(summary(nc1)$NumClones, summary(nc2)$NumClones, alternative = "greater")$p.value < p.value.threshold)
        T2 <- (wilcox.test(summary(nc2)$NumClones, summary(nc3)$NumClones, alternative = "greater")$p.value < p.value.threshold)
        summary(nc1)[, c(1:3, 8:9)]
        summary(nc2)[, c(1:3, 8:9)]
        summary(nc3)[, c(1:3, 8:9)]
        T3 <- (t.test(mutsPerClone(nc1), mutsPerClone(nc2), alternative = "greater")$p.value < p.value.threshold)
        T4 <- (t.test(mutsPerClone(nc2), mutsPerClone(nc3), alternative = "greater")$p.value < p.value.threshold)
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
        ## Can occasionally blow up with pE.f: pE not finite.
        cat("\n mcx2cd: a runif is", runif(1), "\n")
        pops <- 60
        ni <- rep(0, 50)
        names(ni) <- c("a", "b", "c", "d", paste0("n", 1:46))
        fe <- allFitnessEffects(noIntGenes = ni)
        fm6 <- allMutatorEffects(noIntGenes = c("a" = .05,
                                                "b" = 1,
                                                "c" = 10,
                                                "d" = 50))
        nca <- oncoSimulPop(pops, fe, muEF = fm6, finalTime =50,
                            mutationPropGrowth = FALSE,
                            initSize = 1e4,
                            initMutant = "a",
                            sampleEvery = 0.01,
                            keepEvery = 5,
                            onlyCancer = FALSE, detectionProb = NA, seed = NULL, mc.cores = 2, model = "McFL")
        ncb <- oncoSimulPop(pops, fe, muEF = fm6, finalTime =50,
                            mutationPropGrowth = FALSE,
                            initSize = 1e4,
                            initMutant = "b",
                            sampleEvery = 0.01,
                            keepEvery = 5,
                            onlyCancer = FALSE, detectionProb = NA, seed = NULL, mc.cores = 2, model = "McFL")
        ncc <- oncoSimulPop(pops, fe, muEF = fm6, finalTime =50,
                            mutationPropGrowth = FALSE,
                            initSize = 1e4,
                            initMutant = "c",
                            sampleEvery = 0.01,
                            keepEvery = 5,
                            onlyCancer = FALSE, detectionProb = NA, seed = NULL, mc.cores = 2, model = "McFL")
        ncd <- oncoSimulPop(pops, fe, muEF = fm6, finalTime =50,
                            mutationPropGrowth = FALSE,
                            initSize = 1e4,
                            initMutant = "d",
                            sampleEvery = 0.01,
                            keepEvery = 5,
                            onlyCancer = FALSE, detectionProb = NA, seed = NULL, mc.cores = 2, model = "McFL")
        ## I once saw a weird thing
        expect_true(var(summary(nca)$NumClones) > 1e-4)
        expect_true(var(summary(ncb)$NumClones) > 1e-4)
        expect_true(var(summary(ncc)$NumClones) > 1e-4)
        expect_true(var(summary(ncd)$NumClones) > 1e-4)
        ## These are the real tests
        T1 <- ( wilcox.test(summary(nca)$NumClones,
                                 summary(ncb)$NumClones, alternative = "less")$p.value < p.value.threshold)
        T2 <- (wilcox.test(summary(ncb)$NumClones,
                                summary(ncc)$NumClones, alternative = "less")$p.value < p.value.threshold)
        T3 <- ( wilcox.test(summary(ncc)$NumClones,
                                 summary(ncd)$NumClones, alternative = "less")$p.value < p.value.threshold)
        T4 <- (t.test(mutsPerClone(nca), mutsPerClone(ncb), alternative = "less")$p.value < p.value.threshold)
        T5 <- (t.test(mutsPerClone(ncb), mutsPerClone(ncc), alternative = "less")$p.value < p.value.threshold)
        T6 <- (t.test(mutsPerClone(ncc), mutsPerClone(ncd), alternative = "less")$p.value < p.value.threshold)
        if(T1 && T2 && T3 && T4 && T5 && T6 && T7 && T8) break;
    }
    cat(paste("\n done tries", tries, "\n"))
    expect_true(T1 && T2 && T3 && T4 && T5 && T6 && T7 && T8)

})
date()


## Slow (~ 5 seconds) but tests modules of mutator nicely.

date() 
test_that("Mutator modules differences", {
    max.tries <- 4
    for(tries in 1:max.tries) {
        T1 <- T2 <- T3 <- T4 <- T5 <- T6 <- T7 <- T8 <- TRUE
        cat("\n mmd1: a runif is", runif(1), "\n")
        reps <- 60
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
        b1 <- oncoSimulPop(reps,
                           f1,
                           mu = mu,
                           muEF = mut1,
                           onlyCancer = FALSE, detectionProb = NA,
                           initSize = no,
                           finalTime = ft,
                           sampleEvery = 0.01,
                           keepEvery = 5,
                           seed = NULL, mc.cores = 2
                           )
        gc()
        b2 <- oncoSimulPop(reps,
                           f1,
                           mu = mu,
                           muEF = mut2,
                           onlyCancer = FALSE, detectionProb = NA,
                           initSize = no,
                           finalTime = ft,
                           sampleEvery = 0.01,
                           keepEvery = 5,
                           seed = NULL, mc.cores = 2
                           )
        gc()
        b3 <- oncoSimulPop(reps,
                           f1,
                           mu = mu,
                           muEF = mut3,
                           onlyCancer = FALSE, detectionProb = NA,
                           initSize = no,
                           finalTime = ft,
                           sampleEvery = 0.01,
                           keepEvery = 5,
                           seed = NULL, mc.cores = 2
                           )
        gc()
        summary(b3)[, c(1:3, 8:9)]
        summary(b2)[, c(1:3, 8:9)]
        summary(b1)[, c(1:3, 8:9)]
        ## mean(mutsPerClone(b3))
        ## mean(mutsPerClone(b2))
        ## mean(mutsPerClone(b1))
        ## This is, of course, affected by sampling only at end: we do not see
        ## the many intermediate events.
        T1 <- ( wilcox.test(summary(b3)$NumClones, 
                                 summary(b2)$NumClones, alternative = "greater")$p.value < p.value.threshold)
        T2 <- ( wilcox.test(summary(b2)$NumClones, 
                                 summary(b1)$NumClones, alternative = "greater")$p.value < p.value.threshold)
        T3 <- (t.test(mutsPerClone(b3), mutsPerClone(b2), alternative = "greater")$p.value < p.value.threshold)
        T4 <- (t.test(mutsPerClone(b2), mutsPerClone(b1), alternative = "greater")$p.value < p.value.threshold)
        if(T1 && T2 && T3 && T4 && T5 && T6 && T7 && T8) break;
    }
    cat(paste("\n done tries", tries, "\n"))
    expect_true(T1 && T2 && T3 && T4 && T5 && T6 && T7 && T8)

})
date()


date() 
test_that("McFL: Mutator modules differences", {
    max.tries <- 4
    for(tries in 1:max.tries) {
        T1 <- T2 <- T3 <- T4 <- T5 <- T6 <- T7 <- T8 <- TRUE
        cat("\n MCFLmmd1: a runif is", runif(1), "\n")
        reps <- 60
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
        b1 <- oncoSimulPop(reps,
                           f1,
                           mu = mu,
                           muEF = mut1,
                           onlyCancer = FALSE, detectionProb = NA,
                           initSize = no,
                           finalTime = ft,
                           sampleEvery = 0.01,
                           keepEvery = 5, model = "McFL",
                           seed = NULL, mc.cores = 2
                           )
        gc()
        b2 <- oncoSimulPop(reps,
                           f1,
                           mu = mu,
                           muEF = mut2,
                           onlyCancer = FALSE, detectionProb = NA,
                           initSize = no,
                           finalTime = ft,
                           sampleEvery = 0.01,
                           keepEvery = 5, model = "McFL",
                           seed = NULL, mc.cores = 2
                           )
        gc()
        b3 <- oncoSimulPop(reps,
                           f1,
                           mu = mu,
                           muEF = mut3,
                           onlyCancer = FALSE, detectionProb = NA,
                           initSize = no,
                           finalTime = ft,
                           sampleEvery = 0.01,
                           keepEvery = 5, model = "McFL",
                           seed = NULL, mc.cores = 2
                           )
        gc()
        summary(b3)[, c(1:3, 8:9)]
        summary(b2)[, c(1:3, 8:9)]
        summary(b1)[, c(1:3, 8:9)]
        ## mean(mutsPerClone(b3))
        ## mean(mutsPerClone(b2))
        ## mean(mutsPerClone(b1))
        ## This is, of course, affected by sampling only at end: we do not see
        ## the many intermediate events.
        T1 <- ( wilcox.test(summary(b3)$NumClones, 
                                 summary(b2)$NumClones, alternative = "greater")$p.value < p.value.threshold)
        T2 <- ( wilcox.test(summary(b2)$NumClones, 
                                 summary(b1)$NumClones, alternative = "greater")$p.value < p.value.threshold)
        T3 <- (t.test(mutsPerClone(b3), mutsPerClone(b2), alternative = "greater")$p.value < p.value.threshold)
        T4 <- (t.test(mutsPerClone(b2), mutsPerClone(b1), alternative = "greater")$p.value < p.value.threshold)
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
        ## Can occasionally blow up with pE.f: pE not finite.
        cat("\n x1-910: a runif is", runif(1), "\n")
        pops <- 80
        fe <- allFitnessEffects(noIntGenes = c("a" = 0.12,
                                               "b" = 0.14,
                                               "c" = 0.16,
                                               "d" = 0.11))
        fm6 <- allMutatorEffects(noIntGenes = c("a" = 5,
                                                "b" = 10,
                                                "c" = 12,
                                                "d" = 14))
        nc1 <- oncoSimulPop(pops, fe, muEF = fm6, finalTime =250,
                            mutationPropGrowth = FALSE,
                            sampleEvery = 0.001,
                            keepEvery = 5,
                            initSize = 1e6, mc.cores = 2,
                            onlyCancer = FALSE, detectionProb = NA, seed = NULL)
        fm8 <- allMutatorEffects(noIntGenes = c("a" = 1,
                                                "b" = 1,
                                                "c" = 1,
                                                "d" = 1))
        nc2 <- oncoSimulPop(pops, fe, muEF = fm8, finalTime =250,
                            mutationPropGrowth = FALSE,
                            sampleEvery = 0.001,
                            keepEvery = 5,
                            initSize = 1e6, mc.cores = 2,
                            onlyCancer = FALSE, detectionProb = NA, seed = NULL)
        fm7 <- allMutatorEffects(noIntGenes = c("a" = 1e-3,
                                                "b" = 1e-3,
                                                "c" = 1e-3,
                                                "d" = 1e-3))
        nc3 <- oncoSimulPop(pops, fe, muEF = fm7, finalTime =250,
                            mutationPropGrowth = FALSE,
                            sampleEvery = 0.001,
                            keepEvery = 5,
                            initSize = 1e6, mc.cores = 2,
                            onlyCancer = FALSE, detectionProb = NA, seed = NULL)
        T1 <- (wilcox.test(summary(nc1)$NumClones, summary(nc2)$NumClones, alternative = "greater")$p.value < p.value.threshold)
        T2 <- (wilcox.test(summary(nc2)$NumClones, summary(nc3)$NumClones, alternative = "greater")$p.value < p.value.threshold)
        T3 <- (t.test(mutsPerClone(nc1), mutsPerClone(nc2), alternative = "greater")$p.value < p.value.threshold)
        T4 <- (t.test(mutsPerClone(nc2), mutsPerClone(nc3), alternative = "greater")$p.value < p.value.threshold)
        summary(nc1)[, c(1:3, 8:9)]
        summary(nc2)[, c(1:3, 8:9)]
        summary(nc3)[, c(1:3, 8:9)]
        if(T1 && T2 && T3 && T4 && T5 && T6 && T7 && T8) break;
    }
    cat(paste("\n done tries", tries, "\n"))
    expect_true(T1 && T2 && T3 && T4 && T5 && T6 && T7 && T8)

})
date()


## very slow, because huge number of clones. But tests several phenomena comprehensively.
## same with McFL below

date()
test_that("per-gene-mut rates and mutator", {
    max.tries <- 4
    for(tries in 1:max.tries) {
        T1 <- T2 <- T3 <- T4 <- T5 <- T6 <- T7 <- T8 <- TRUE
        cat("\n oss11-x970: a runif is", runif(1), "\n")
        ng <- 10
        ni <- rep(0, ng)
        m1 <- runif(ng, min = 1e-7, max = 5e-6)
        m2 <- runif(ng, min = 1e-5, max = 1e-4)
        names(ni) <- names(m1) <- names(m2) <- c(replicate(ng,
                                                           paste(sample(letters, 12), collapse = "")))
        fe1 <- allFitnessEffects(noIntGenes = ni)
        ft <- 50
        no <- 5e5 
        reps <- 100
        gn <- paste(names(ni), collapse = ", ")
        ## MUs used to be 25 and 100. Way too slow.
        mutator1 <- allMutatorEffects(epistasis = c("MU" = 15),
                                      geneToModule = c("MU" = gn))
        mutator2 <- allMutatorEffects(epistasis = c("MU" = 30),
                                      geneToModule = c("MU" = gn))
        m1.mutator0 <- oncoSimulPop(reps,
                                    fe1,
                                    mu = m1,
                                    onlyCancer = FALSE, detectionProb = NA,
                                    initSize = no,
                                    finalTime = ft,
                                    sampleEvery = 0.01,
                                    keepEvery = 5,
                                    seed = NULL, mc.cores = 2
                                    )
        m1.mutator1 <- oncoSimulPop(reps,
                                    fe1,
                                    mu = m1,
                                    muEF = mutator1,
                                    onlyCancer = FALSE, detectionProb = NA,
                                    initSize = no,
                                    finalTime = ft,
                                    sampleEvery = 0.01,
                                    keepEvery = 5,
                                    seed = NULL, mc.cores = 2
                                    )
        runif(1)
        m1.mutator2 <- oncoSimulPop(reps,
                                    fe1,
                                    mu = m1,
                                    muEF = mutator2,
                                    onlyCancer = FALSE, detectionProb = NA,
                                    initSize = no,
                                    finalTime = ft,
                                    sampleEvery = 0.01,
                                    keepEvery = 5,
                                    seed = NULL, mc.cores = 2
                                    )
        m2.mutator0 <- oncoSimulPop(reps,
                                    fe1,
                                    mu = m2,
                                    onlyCancer = FALSE, detectionProb = NA,
                                    initSize = no,
                                    finalTime = ft,
                                    sampleEvery = 0.01,
                                    keepEvery = 5,
                                    seed = NULL, mc.cores = 2
                                    )
        m2.mutator1 <- oncoSimulPop(reps,
                                    fe1,
                                    mu = m2,
                                    muEF = mutator1,
                                    onlyCancer = FALSE, detectionProb = NA,
                                    initSize = no,
                                    finalTime = ft,
                                    sampleEvery = 0.01,
                                    keepEvery = 5,
                                    seed = NULL, mc.cores = 2
                                    )
        m2.mutator2 <- oncoSimulPop(reps,
                                    fe1,
                                    mu = m2,
                                    muEF = mutator2,
                                    onlyCancer = FALSE, detectionProb = NA,
                                    initSize = no,
                                    finalTime = ft,
                                    sampleEvery = 0.01,
                                    keepEvery = 5,
                                    seed = NULL, mc.cores = 2
                                    )
        summary(m1.mutator0)[, c(1:3, 8:9)]
        summary(m1.mutator1)[, c(1:3, 8:9)]
        summary(m1.mutator2)[, c(1:3, 8:9)]
        summary(m2.mutator0)[, c(1:3, 8:9)]
        summary(m2.mutator1)[, c(1:3, 8:9)]
        summary(m2.mutator2)[, c(1:3, 8:9)]
        ## Mutator increases if larger mutator and compared to no mutator
        ## within levels of per-gene mutation rates
        T1 <- ( wilcox.test(summary(m1.mutator2)$NumClones, 
                                 summary(m1.mutator1)$NumClones, alternative = "greater")$p.value < p.value.threshold)
        T2 <- ( wilcox.test(summary(m1.mutator1)$NumClones, 
                                 summary(m1.mutator0)$NumClones, alternative = "greater")$p.value < p.value.threshold)
        T3 <- ( wilcox.test(summary(m2.mutator2)$NumClones, 
                                 summary(m2.mutator1)$NumClones, alternative = "greater")$p.value < p.value.threshold)
        T4 <- ( wilcox.test(summary(m2.mutator1)$NumClones, 
                                 summary(m2.mutator0)$NumClones, alternative = "greater")$p.value < p.value.threshold)
        T5 <- (t.test(mutsPerClone(m1.mutator2), mutsPerClone(m1.mutator1), alternative = "greater")$p.value < p.value.threshold)
        T6 <- (t.test(mutsPerClone(m1.mutator1), mutsPerClone(m1.mutator0), alternative = "greater")$p.value < p.value.threshold)
        T7 <- (t.test(mutsPerClone(m2.mutator2), mutsPerClone(m2.mutator1), alternative = "greater")$p.value < p.value.threshold)
        T8 <- (t.test(mutsPerClone(m2.mutator1), mutsPerClone(m2.mutator0), alternative = "greater")$p.value < p.value.threshold)
        ## Increases in mutation rates increase clones, etc, within levels of
        ## mutator.
        Ta <- ( wilcox.test(summary(m1.mutator0)$NumClones,
                                 summary(m2.mutator0)$NumClones, alternative = "less")$p.value < p.value.threshold)
        Tb <- ( wilcox.test(summary(m1.mutator1)$NumClones,
                                 summary(m2.mutator1)$NumClones, alternative = "less")$p.value < p.value.threshold)
        Tc <- ( wilcox.test(summary(m1.mutator2)$NumClones,
                                 summary(m2.mutator2)$NumClones, alternative = "less")$p.value < p.value.threshold)
        Td <- (t.test(mutsPerClone(m1.mutator0), mutsPerClone(m2.mutator0), alternative = "less")$p.value < p.value.threshold)
        Te <- (t.test(mutsPerClone(m1.mutator1), mutsPerClone(m2.mutator1), alternative = "less")$p.value < p.value.threshold)
        Tf <- (t.test(mutsPerClone(m1.mutator2), mutsPerClone(m2.mutator2), alternative = "less")$p.value < p.value.threshold)
        if(T1 && T2 && T3 && T4 && T5 && T6 && T7 && T8 &&
           Ta && Tb && Tc && Td && Te && Tf) break;
    }
    cat(paste("\n done tries", tries, "\n"))
    expect_true(T1 && T2 && T3 && T4 && T5 && T6 && T7 && T8 &&
                Ta && Tb && Tc && Td && Te && Tf)

})
date()


date()
test_that("McFL: per-gene-mut rates and mutator", {
    max.tries <- 4
    for(tries in 1:max.tries) {
        T1 <- T2 <- T3 <- T4 <- T5 <- T6 <- T7 <- T8 <- TRUE
        cat("\n mcfloss11: a runif is", runif(1), "\n")
        ng <- 10
        ni <- rep(0, ng)
        m1 <- runif(ng, min = 1e-7, max = 5e-6)
        m2 <- runif(ng, min = 1e-5, max = 1e-4)
        names(ni) <- names(m1) <- names(m2) <- c(replicate(ng,
                                                           paste(sample(letters, 12), collapse = "")))
        fe1 <- allFitnessEffects(noIntGenes = ni)
        ft <- 50
        no <- 5e5 
        reps <- 150
        gn <- paste(names(ni), collapse = ", ")
        mutator1 <- allMutatorEffects(epistasis = c("MU" = 15),
                                      geneToModule = c("MU" = gn))
        mutator2 <- allMutatorEffects(epistasis = c("MU" = 30),
                                      geneToModule = c("MU" = gn))
        m1.mutator0 <- oncoSimulPop(reps,
                                    fe1,
                                    mu = m1,
                                    onlyCancer = FALSE, detectionProb = NA,
                                    initSize = no,
                                    finalTime = ft,
                                    sampleEvery = 0.01,
                                    keepEvery = 5,
                                    seed = NULL, mc.cores = 2, model = "McFL"
                                    )
        m1.mutator1 <- oncoSimulPop(reps,
                                    fe1,
                                    mu = m1,
                                    muEF = mutator1,
                                    onlyCancer = FALSE, detectionProb = NA,
                                    initSize = no,
                                    finalTime = ft,
                                    sampleEvery = 0.01,
                                    keepEvery = 5,
                                    seed = NULL, mc.cores = 2, model = "McFL"
                                    )
        m1.mutator2 <- oncoSimulPop(reps,
                                    fe1,
                                    mu = m1,
                                    muEF = mutator2,
                                    onlyCancer = FALSE, detectionProb = NA,
                                    initSize = no,
                                    finalTime = ft,
                                    sampleEvery = 0.01,
                                    keepEvery = 5,
                                    seed = NULL, mc.cores = 2, model = "McFL"
                                    )
        cat("\n starting m2\n")
        m2.mutator0 <- oncoSimulPop(reps,
                                    fe1,
                                    mu = m2,
                                    onlyCancer = FALSE, detectionProb = NA,
                                    initSize = no,
                                    finalTime = ft,
                                    sampleEvery = 0.01,
                                    keepEvery = 5,
                                    seed = NULL, mc.cores = 2, model = "McFL"
                                    )
        m2.mutator1 <- oncoSimulPop(reps,
                                    fe1,
                                    mu = m2,
                                    muEF = mutator1,
                                    onlyCancer = FALSE, detectionProb = NA,
                                    initSize = no,
                                    finalTime = ft,
                                    sampleEvery = 0.01,
                                    keepEvery = 5,
                                    seed = NULL, mc.cores = 2, model = "McFL"
                                    )
        m2.mutator2 <- oncoSimulPop(reps,
                                    fe1,
                                    mu = m2,
                                    muEF = mutator2,
                                    onlyCancer = FALSE, detectionProb = NA,
                                    initSize = no,
                                    finalTime = ft,
                                    sampleEvery = 0.01,
                                    keepEvery = 5,
                                    seed = NULL, mc.cores = 2, model = "McFL"
                                    )
        summary(m1.mutator0)[, c(1:3, 8:9)]
        summary(m1.mutator1)[, c(1:3, 8:9)]
        summary(m1.mutator2)[, c(1:3, 8:9)]
        summary(m2.mutator0)[, c(1:3, 8:9)]
        summary(m2.mutator1)[, c(1:3, 8:9)]
        summary(m2.mutator2)[, c(1:3, 8:9)]
        ## Mutator increases if larger mutator and compared to no mutator
        ## within levels of per-gene mutation rates
        T1 <- ( wilcox.test(summary(m1.mutator2)$NumClones, 
                                 summary(m1.mutator1)$NumClones, alternative = "greater")$p.value < p.value.threshold)
        T2 <- ( wilcox.test(summary(m1.mutator1)$NumClones, 
                                 summary(m1.mutator0)$NumClones, alternative = "greater")$p.value < p.value.threshold)
        T3 <- ( wilcox.test(summary(m2.mutator2)$NumClones, 
                                 summary(m2.mutator1)$NumClones, alternative = "greater")$p.value < p.value.threshold)
        T4 <- ( wilcox.test(summary(m2.mutator1)$NumClones, 
                                 summary(m2.mutator0)$NumClones, alternative = "greater")$p.value < p.value.threshold)
        T5 <- (t.test(mutsPerClone(m1.mutator2), mutsPerClone(m1.mutator1), alternative = "greater")$p.value < p.value.threshold)
        T6 <- (t.test(mutsPerClone(m1.mutator1), mutsPerClone(m1.mutator0), alternative = "greater")$p.value < p.value.threshold)
        T7 <- (t.test(mutsPerClone(m2.mutator2), mutsPerClone(m2.mutator1), alternative = "greater")$p.value < p.value.threshold)
        T8 <- (t.test(mutsPerClone(m2.mutator1), mutsPerClone(m2.mutator0), alternative = "greater")$p.value < p.value.threshold)
        ## Increases in mutation rates increase clones, etc, within levels of
        ## mutator.
        Ta <- ( wilcox.test(summary(m1.mutator0)$NumClones,
                                 summary(m2.mutator0)$NumClones, alternative = "less")$p.value < p.value.threshold)
        Tb <- ( wilcox.test(summary(m1.mutator1)$NumClones,
                                 summary(m2.mutator1)$NumClones, alternative = "less")$p.value < p.value.threshold)
        Tc <- ( wilcox.test(summary(m1.mutator2)$NumClones,
                                 summary(m2.mutator2)$NumClones, alternative = "less")$p.value < p.value.threshold)
        Td <- (t.test(mutsPerClone(m1.mutator0), mutsPerClone(m2.mutator0), alternative = "less")$p.value < p.value.threshold)
        Te <- (t.test(mutsPerClone(m1.mutator1), mutsPerClone(m2.mutator1), alternative = "less")$p.value < p.value.threshold)
        Tf <- (t.test(mutsPerClone(m1.mutator2), mutsPerClone(m2.mutator2), alternative = "less")$p.value < p.value.threshold)
        if(T1 && T2 && T3 && T4 && T5 && T6 && T7 && T8 &&
           Ta && Tb && Tc && Td && Te && Tf) break;
    }
    cat(paste("\n done tries", tries, "\n"))
    expect_true(T1 && T2 && T3 && T4 && T5 && T6 && T7 && T8 &&
           Ta && Tb && Tc && Td && Te && Tf)

})
date()


date()
test_that("Mutator, several modules differences, McFL", {
    max.tries <- 4
    for(tries in 1:max.tries) {
        T1 <- T2 <- T3 <- T4 <- T5 <- T6 <- T7 <- T8 <- TRUE
        cat("\n mmdSMMC1: a runif is", runif(1), "\n")
        reps <- 50
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
        b1 <- oncoSimulPop(reps,
                           f1,
                           mu = mu,
                           muEF = mut1,
                           onlyCancer = FALSE, detectionProb = NA,
                           initSize = no,
                           finalTime = ft,
                           sampleEvery = 0.01, 
                           detectionSize = 1e9,
                           detectionDrivers = 9999,
                           seed = NULL, mc.cores = 2, model = "McFL"
                           )
        gc()
        b2 <- oncoSimulPop(reps,
                           f1,
                           mu = mu,
                           muEF = mut2,
                           onlyCancer = FALSE, detectionProb = NA,
                           initSize = no,
                           finalTime = ft,
                           sampleEvery = 0.01, 
                           detectionSize = 1e9,
                           detectionDrivers = 9999,
                           seed = NULL, mc.cores = 2, model = "McFL"
                           )
        gc()
        summary(b2)[, c(1:3, 8:9)]
        summary(b1)[, c(1:3, 8:9)]
        mean(mutsPerClone(b2))
        mean(mutsPerClone(b1))
        ## This is, of course, affected by sampling only at end: we do not see
        ## the many intermediate events.
        T1 <- ( wilcox.test(summary(b2)$NumClones, 
                                 summary(b1)$NumClones, alternative = "greater")$p.value < p.value.threshold)
        T2 <- (t.test(mutsPerClone(b2), mutsPerClone(b1), alternative = "greater")$p.value < p.value.threshold)
        if(T1 && T2 && T3 && T4 && T5 && T6 && T7 && T8) break;
    }
    cat(paste("\n done tries", tries, "\n"))
    expect_true(T1 && T2 && T3 && T4 && T5 && T6 && T7 && T8)

})
date()


date()
test_that("McFL: Expect freq genotypes, mutator and var mut rates", {
    max.tries <- 4
    for(tries in 1:max.tries) {
        T1 <- T2 <- T3 <- T4 <- T5 <- T6 <- T7 <- T8 <- TRUE
        ## increase mutator
        cat("\n u8: a runif is", runif(1), "\n")
        pops <- 300
        ft <- .0001
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
        mutator1 <- c("oreoisasabgene" = 50)
        m1 <- allMutatorEffects(noIntGenes = mutator1)
        pg1["hereisoneagene"] <- 1e-3 ## to compare with a laarge one
        m1.pg1.b <- oncoSimulPop(pops,
                                 fe,
                                 mu = pg1,
                                 muEF = m1,
                                 model = "McFL",
                                 finalTime = ft,
                                 mutationPropGrowth = FALSE,
                                 initSize = no,
                                 initMutant ="oreoisasabgene",
                                 onlyCancer = FALSE, detectionProb = NA, seed = NULL, mc.cores = 2)
        expect_true(sm("oreoisasabgene", m1.pg1.b) == totalind(m1.pg1.b))
        ## enom("oreoisasabgene", pg1)
        ## snom("oreoisasabgene", m1.pg1.b)
        p.fail <- 1e-3
        T1 <- (chisq.test(snom("oreoisasabgene", m1.pg1.b),
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
        ## Similar to above, but mutator has a single element, not the whole
        ## vector.
        cat("\n u7_1: a runif is", runif(1), "\n")
        pops <- 300
        ft <- .0001
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
        mutator1 <- c("oreoisasabgene" = 50)
        m1 <- allMutatorEffects(noIntGenes = mutator1)
        pg1["hereisoneagene"] <- 1e-3 ## to compare with a laarge one
        m1.pg1.b <- oncoSimulPop(pops,
                                 fe,
                                 mu = pg1,
                                 muEF = m1,
                                 finalTime = ft,
                                 mutationPropGrowth = FALSE,
                                 initSize = no,
                                 initMutant ="oreoisasabgene",
                                 onlyCancer = FALSE, detectionProb = NA, seed = NULL, mc.cores = 2)
        expect_true(sm("oreoisasabgene", m1.pg1.b) == totalind(m1.pg1.b))
        enom("oreoisasabgene", pg1, no, pops)
        snom("oreoisasabgene", m1.pg1.b)
        p.fail <- 1e-3
        T1 <- (chisq.test(snom("oreoisasabgene", m1.pg1.b),
                               p = pnom("oreoisasabgene", pg1, no, pops))$p.value > p.fail)
        if(T1 && T2 && T3 && T4 && T5 && T6 && T7 && T8) break;
    }
    cat(paste("\n done tries", tries, "\n"))
    expect_true(T1 && T2 && T3 && T4 && T5 && T6 && T7 && T8)

})
date()


date()
test_that("Mutator and mutPropGrowth", {
    max.tries <- 4
    for(tries in 1:max.tries) {
        T1 <- T2 <- T3 <- T4 <- T5 <- T6 <- T7 <- T8 <- TRUE
        ## we stop on size
        cat("\n mmpg: a runif is", runif(1), "\n")
        reps <- 200
        no <- 5e3
        ds <- 7e4
        ft <- 1000 
        mu <- 1e-5
        ln <- 60
        m1 <- 1 ## if this is too large, easy to get it to blow.
        m50 <- 50
        gna <- paste0("a", 1:ln)
        gna <- paste(gna, collapse = ", ")
        f1 <- allFitnessEffects(epistasis = c("A" = 0.4),
                                geneToModule = c("A" = gna))
        mut1 <- allMutatorEffects(epistasis = c("M" = m1),
                                  geneToModule = c("M" = gna))
        mut50 <- allMutatorEffects(epistasis = c("M" = m50),
                                   geneToModule = c("M" = gna))
        m1.pg <- oncoSimulPop(reps,
                              f1,
                              mu = mu,
                              muEF = mut1,
                              mutationPropGrowth = TRUE,
                              onlyCancer = FALSE, detectionProb = NA,
                              initSize = no,
                              finalTime = ft,
                              detectionSize = ds,
                              sampleEvery = 0.01,
                              keepEvery = 5,
                              seed = NULL, mc.cores = 2
                              )
        ## gc()
        m1.npg <- oncoSimulPop(reps,
                               f1,
                               mu = mu,
                               muEF = mut1,
                               mutationPropGrowth = FALSE,
                               onlyCancer = FALSE, detectionProb = NA,
                               initSize = no,
                               finalTime = ft,
                               detectionSize = ds,
                               sampleEvery = 0.01,
                               keepEvery = 5,
                               seed = NULL, mc.cores = 2
                               )
        ##gc()
        m50.pg <- oncoSimulPop(reps,
                               f1,
                               mu = mu,
                               muEF = mut50,
                               mutationPropGrowth = TRUE,
                               onlyCancer = FALSE, detectionProb = NA,
                               initSize = no,
                               finalTime = ft,
                               detectionSize = ds,
                               sampleEvery = 0.01,
                               keepEvery = 5,
                               seed = NULL, mc.cores = 2
                               )
        ## gc()
        m50.npg <- oncoSimulPop(reps,
                                f1,
                                mu = mu,
                                muEF = mut50,
                                mutationPropGrowth = FALSE,
                                onlyCancer = FALSE, detectionProb = NA,
                                initSize = no,
                                finalTime = ft,
                                detectionSize = ds,
                                sampleEvery = 0.01,
                                keepEvery = 5,
                                seed = NULL, mc.cores = 2
                                )
        ##gc()
        summary(m1.pg)[, c(1:3, 8:9)]
        summary(m50.pg)[, c(1:3, 8:9)]
        summary(m1.npg)[, c(1:3, 8:9)]
        summary(m50.npg)[, c(1:3, 8:9)]
        ## Over mutator, as we have mutPropGrowth, clones, etc, increase
        T1 <- ( wilcox.test(summary(m1.pg)$NumClones,
                                 summary(m1.npg)$NumClones,
                                 alternative = "greater")$p.value < p.value.threshold)
        T2 <- (t.test(mutsPerClone(m1.pg),
                           mutsPerClone(m1.npg),
                           alternative = "greater")$p.value < p.value.threshold)
        T3 <- ( wilcox.test(summary(m50.pg)$NumClones,
                                 summary(m50.npg)$NumClones,
                                 alternative = "greater")$p.value < p.value.threshold)
        T4 <- (t.test(mutsPerClone(m50.pg),
                           mutsPerClone(m50.npg),
                           alternative = "greater")$p.value < p.value.threshold)
        ## Over mutPropGrowth, as we increase mutator, clones, etc, increase
        T5 <- ( wilcox.test(summary(m50.pg)$NumClones,
                                 summary(m1.pg)$NumClones,
                                 alternative = "greater")$p.value < p.value.threshold)
        T6 <- (t.test(mutsPerClone(m50.pg),
                           mutsPerClone(m1.pg),
                           alternative = "greater")$p.value < p.value.threshold)
        T7 <- ( wilcox.test(summary(m50.npg)$NumClones,
                                 summary(m1.npg)$NumClones,
                                 alternative = "greater")$p.value < p.value.threshold)
        T8 <- (t.test(mutsPerClone(m50.npg),
                           mutsPerClone(m1.npg),
                           alternative = "greater")$p.value < p.value.threshold)
        if(T1 && T2 && T3 && T4 && T5 && T6 && T7 && T8) break;
    }
    cat(paste("\n done tries", tries, "\n"))
    expect_true(T1 && T2 && T3 && T4 && T5 && T6 && T7 && T8)

})
date()



cat(paste("\n Finished test.mutator-long.R test at", date(), "\n"))
































## ### These are too convoluted and trying to pick too tiny diffs. Substitute
## ###  by diffs initMutant and by per-gene-mut rates and mutator

## date()
## test_that("Num clones: Mutator and var mut rates and init and really scrambled names", {
##     
##     
##     cat("\n z2: a runif is", runif(1), "\n")
##     pops <- 40
##     ft <- .005
##     lni <- 50
##     no <- 1e7
##     ni <- c(0, 0, 0, rep(0, lni))
##     ## scramble around names
##     names(ni) <- c("hereisoneagene",
##                    "oreoisasabgene",
##                    "nnhsisthecgene",
##                    replicate(lni,
##                              paste(sample(letters, 12), collapse = "")))
##     ni <- ni[order(names(ni))]
##     fe <- allFitnessEffects(noIntGenes = ni)
##     mutator1 <- mutator2 <- rep(1, lni + 3)
##     pg1 <- pg2 <- rep(5e-6, lni + 3)
##     ## scramble names of mutator and per-gene too
##     names(mutator1) <- sample(names(ni))
##     names(mutator2) <- sample(names(ni))
##     names(pg1) <- sample(names(ni))
##     names(pg2) <- sample(names(ni))
##     mutator1["hereisoneagene"] <- 100
##     mutator2["hereisoneagene"] <- 1
##     mutator1["oreoisasabgene"] <- 1
##     mutator2["oreoisasabgene"] <- 0.01
##     mutator1["nnhsisthecgene"] <- 0.01
##     mutator2["nnhsisthecgene"] <- 100
##     m1 <- allMutatorEffects(noIntGenes = mutator1)
##     m2 <- allMutatorEffects(noIntGenes = mutator2)
##     pg1["hereisoneagene"] <- 1e-3
##     pg2["hereisoneagene"] <- 1e-14
##     pg1["oreoisasabgene"] <- 1e-7
##     pg2["oreoisasabgene"] <- 1e-7
##     pg1["nnhsisthecgene"] <- 1e-14
##     pg2["nnhsisthecgene"] <- 1e-3
##     m1.pg1.a <- oncoSimulPop(pops,
##                            fe,
##                            mu = pg1,
##                            muEF = m1,
##                            finalTime = ft,
##                            mutationPropGrowth = FALSE,
##                            initSize = no,
##                            initMutant = "hereisoneagene",
##                            onlyCancer = FALSE, detectionProb = NA, mc.cores = 2)
##     m1.pg2.a <- oncoSimulPop(pops,
##                            fe,
##                            mu = pg2,
##                            muEF = m1,
##                            finalTime = ft,
##                            mutationPropGrowth = FALSE,
##                            initSize = no,
##                            initMutant = "hereisoneagene",
##                            onlyCancer = FALSE, detectionProb = NA, mc.cores = 2)
##     m1.pg1.c <- oncoSimulPop(pops,
##                            fe,
##                            mu = pg1,
##                            muEF = m1,
##                            finalTime = ft,
##                            mutationPropGrowth = FALSE,
##                            initSize = no,
##                            initMutant = "nnhsisthecgene",
##                            onlyCancer = FALSE, detectionProb = NA, mc.cores = 2)
##     m1.pg2.c <- oncoSimulPop(pops,
##                            fe,
##                            mu = pg2,
##                            muEF = m1,
##                            finalTime = ft,
##                            mutationPropGrowth = FALSE,
##                            initSize = no,
##                            initMutant = "nnhsisthecgene",
##                            onlyCancer = FALSE, detectionProb = NA, mc.cores = 2)
##     m1.pg1.b <- oncoSimulPop(pops,
##                            fe,
##                            mu = pg1,
##                            muEF = m1,
##                            finalTime = ft,
##                            mutationPropGrowth = FALSE,
##                            initSize = no,
##                            initMutant = "oreoisasabgene",
##                            onlyCancer = FALSE, detectionProb = NA, mc.cores = 2)
##     expect_true(medianNClones(m1.pg1.a) < medianNClones(m1.pg2.a))
##     ## Too tiny diffs
##     ## expect_true(medianNClones(m1.pg1.c) > medianNClones(m1.pg2.c))
##     expect_true(medianNClones(m1.pg1.a) > medianNClones(m1.pg1.c))
##     expect_true(medianNClones(m1.pg2.a) > medianNClones(m1.pg2.c))
##     expect_true(medianNClones(m1.pg1.a) > medianNClones(m1.pg1.b))
##     m2.pg1.a <- oncoSimulPop(pops,
##                            fe,
##                            mu = pg1,
##                            muEF = m2,
##                            finalTime = ft,
##                            mutationPropGrowth = FALSE,
##                            initSize = no,
##                            initMutant = "hereisoneagene",
##                            onlyCancer = FALSE, detectionProb = NA, mc.cores = 2)
##     m2.pg2.a <- oncoSimulPop(pops,
##                            fe,
##                            mu = pg2,
##                            muEF = m2,
##                            finalTime = ft,
##                            mutationPropGrowth = FALSE,
##                            initSize = no,
##                            initMutant = "hereisoneagene",
##                            onlyCancer = FALSE, detectionProb = NA, mc.cores = 2)
##     m2.pg1.c <- oncoSimulPop(pops,
##                            fe,
##                            mu = pg1,
##                            muEF = m2,
##                            finalTime = ft,
##                            mutationPropGrowth = FALSE,
##                            initSize = no,
##                            initMutant = "nnhsisthecgene",
##                            onlyCancer = FALSE, detectionProb = NA, mc.cores = 2)
##     m2.pg2.c <- oncoSimulPop(pops,
##                            fe,
##                            mu = pg2,
##                            muEF = m2,
##                            finalTime = ft,
##                            mutationPropGrowth = FALSE,
##                            initSize = no,
##                            initMutant = "nnhsisthecgene",
##                            onlyCancer = FALSE, detectionProb = NA, mc.cores = 2)
##     ## too tiny
##     ## expect_true(medianNClones(m2.pg1.a) < medianNClones(m2.pg2.a))
##     expect_true(medianNClones(m2.pg1.c) > medianNClones(m2.pg2.c))
##     expect_true(medianNClones(m2.pg1.a) < medianNClones(m2.pg1.c))
##     expect_true(medianNClones(m2.pg2.a) < medianNClones(m2.pg2.c))
##     expect_true(medianNClones(m1.pg1.c) < medianNClones(m2.pg1.c))
##     expect_true(medianNClones(m1.pg1.c) < medianNClones(m2.pg2.c))
##     expect_true(medianNClones(m2.pg1.a) < medianNClones(m1.pg1.a))
##     expect_true(medianNClones(m2.pg1.a) < medianNClones(m1.pg2.a))
##     expect_true(medianNClones(m2.pg2.a) < medianNClones(m2.pg1.c))
##     expect_true(medianNClones(m1.pg2.a) > medianNClones(m1.pg1.c))
##     expect_true(smAnomPi(m2.pg1.a, "hereisoneagene") <
##                 smAnomPi(m2.pg1.c, "nnhsisthecgene"))
##     expect_true(smAnomPi(m2.pg1.a, "hereisoneagene") <
##                 smAnomPi(m1.pg1.a, "hereisoneagene"))
## })
## date()

## ## This is slow and very fragile, because unless we want to increase pops
## ## a lots, we can often miss differences in
## ## medianNClones(m1.pg1.a) < medianNClones(m1.pg2.a)
## ## and
## ## medianNClones(m1.pg1.c) > medianNClones(m1.pg2.c)
## ## and
## ## medianNClones(m2.pg1.c) > medianNClones(m2.pg2.c)
## date()
## test_that("McFL: Num clones: Mutator and var mut rates and init and really scrambled names", {
##     pseed <- sample(300:400, 1)
##     
##     cat("\n x3: a runif is", runif(1), "\n")
##     pops <- 80
##     ft <- .002
##     lni <- 20
##     no <- 3e5
##     ni <- c(0, 0, 0, rep(0, lni))
##     ## scramble around names
##     names(ni) <- c("hereisoneagene",
##                    "oreoisasabgene",
##                    "nnhsisthecgene",
##                    replicate(lni,
##                              paste(sample(letters, 12), collapse = "")))
##     ni <- ni[order(names(ni))]
##     fe <- allFitnessEffects(noIntGenes = ni)
##     mutator1 <- mutator2 <- rep(1, lni + 3)
##     pg1 <- pg2 <-  rep(1e-4, lni + 3)
##     ## this: runif(lni + 3, min = 1e-6, max = 1e-5) 
##     ## adds too much variation in mutation rates between runs
##     ## scramble names of mutator and per-gene too
##     names(mutator1) <- sample(names(ni))
##     names(mutator2) <- sample(names(ni))
##     names(pg1) <- sample(names(ni))
##     names(pg2) <- sample(names(ni))
##     mutator1["hereisoneagene"] <- 50
##     mutator2["hereisoneagene"] <- 5 # 1
##     mutator1["oreoisasabgene"] <- 1
##     mutator2["oreoisasabgene"] <- 0.01
##     mutator1["nnhsisthecgene"] <- 1 # 0.5
##     mutator2["nnhsisthecgene"] <- 50
##     m1 <- allMutatorEffects(noIntGenes = mutator1)
##     m2 <- allMutatorEffects(noIntGenes = mutator2)
##     pg1["hereisoneagene"] <- 5e-3
##     pg2["hereisoneagene"] <- 5e-15
##     pg1["oreoisasabgene"] <- 1e-7
##     pg2["oreoisasabgene"] <- 1e-7
##     pg1["nnhsisthecgene"] <- 5e-15
##     pg2["nnhsisthecgene"] <- 5e-3
##     m1.pg1.a <- oncoSimulPop(pops,
##                            fe,
##                            mu = pg1,
##                            muEF = m1,
##                            finalTime = ft, model = "McFL",
##                            mutationPropGrowth = FALSE,
##                            initSize = no,
##                            initMutant = "hereisoneagene",
##                            onlyCancer = FALSE, detectionProb = NA, mc.cores = 2)
##     m1.pg2.a <- oncoSimulPop(pops,
##                            fe,
##                            mu = pg2,
##                            muEF = m1,
##                            finalTime = ft, model = "McFL",
##                            mutationPropGrowth = FALSE,
##                            initSize = no,
##                            initMutant = "hereisoneagene",
##                            onlyCancer = FALSE, detectionProb = NA, mc.cores = 2)
##     m1.pg1.c <- oncoSimulPop(pops,
##                            fe,
##                            mu = pg1,
##                            muEF = m1,
##                            finalTime = ft, model = "McFL",
##                            mutationPropGrowth = FALSE,
##                            initSize = no,
##                            initMutant = "nnhsisthecgene",
##                            onlyCancer = FALSE, detectionProb = NA, mc.cores = 2)
##     m1.pg2.c <- oncoSimulPop(pops,
##                            fe,
##                            mu = pg2,
##                            muEF = m1,
##                            finalTime = ft, model = "McFL",
##                            mutationPropGrowth = FALSE,
##                            initSize = no,
##                            initMutant = "nnhsisthecgene",
##                            onlyCancer = FALSE, detectionProb = NA, mc.cores = 2)
##     m1.pg1.b <- oncoSimulPop(pops,
##                            fe,
##                            mu = pg1,
##                            muEF = m1,
##                            finalTime = ft, model = "McFL",
##                            mutationPropGrowth = FALSE,
##                            initSize = no,
##                            initMutant = "oreoisasabgene",
##                            onlyCancer = FALSE, detectionProb = NA, mc.cores = 2)
##     expect_true(medianNClones(m1.pg1.a) < medianNClones(m1.pg2.a))
##     ## Too tiny diffs
##     ## expect_true(medianNClones(m1.pg1.c) > medianNClones(m1.pg2.c))
##     expect_true(medianNClones(m1.pg1.a) > medianNClones(m1.pg1.c))
##     expect_true(medianNClones(m1.pg2.a) > medianNClones(m1.pg2.c))
##     expect_true(medianNClones(m1.pg1.a) > medianNClones(m1.pg1.b))
##     m2.pg1.a <- oncoSimulPop(pops,
##                            fe,
##                            mu = pg1,
##                            muEF = m2,
##                            finalTime = ft, model = "McFL",
##                            mutationPropGrowth = FALSE,
##                            initSize = no,
##                            initMutant = "hereisoneagene",
##                            onlyCancer = FALSE, detectionProb = NA, mc.cores = 2)
##     m2.pg2.a <- oncoSimulPop(pops,
##                            fe,
##                            mu = pg2,
##                            muEF = m2,
##                            finalTime = ft, model = "McFL",
##                            mutationPropGrowth = FALSE,
##                            initSize = no,
##                            initMutant = "hereisoneagene",
##                            onlyCancer = FALSE, detectionProb = NA, mc.cores = 2)
##     m2.pg1.c <- oncoSimulPop(pops,
##                            fe,
##                            mu = pg1,
##                            muEF = m2,
##                            finalTime = ft, model = "McFL",
##                            mutationPropGrowth = FALSE,
##                            initSize = no,
##                            initMutant = "nnhsisthecgene",
##                            onlyCancer = FALSE, detectionProb = NA, mc.cores = 2)
##     m2.pg2.c <- oncoSimulPop(pops,
##                            fe,
##                            mu = pg2,
##                            muEF = m2,
##                            finalTime = ft, model = "McFL",
##                            mutationPropGrowth = FALSE,
##                            initSize = no,
##                            initMutant = "nnhsisthecgene",
##                            onlyCancer = FALSE, detectionProb = NA, mc.cores = 2)
##     ## too small.
##     ## expect_true(medianNClones(m2.pg1.a) < medianNClones(m2.pg2.a))
##     expect_true(medianNClones(m2.pg1.c) > medianNClones(m2.pg2.c))
##     expect_true(medianNClones(m2.pg1.a) < medianNClones(m2.pg1.c))
##     expect_true(medianNClones(m2.pg2.a) < medianNClones(m2.pg2.c))
##     expect_true(medianNClones(m1.pg1.c) < medianNClones(m2.pg1.c))
##     expect_true(medianNClones(m1.pg1.c) < medianNClones(m2.pg2.c))
##     expect_true(medianNClones(m2.pg1.a) < medianNClones(m1.pg1.a))
##     expect_true(medianNClones(m2.pg1.a) < medianNClones(m1.pg2.a))
##     expect_true(medianNClones(m2.pg2.a) < medianNClones(m2.pg1.c))
##     expect_true(medianNClones(m1.pg2.a) > medianNClones(m1.pg1.c))
##     expect_true(smAnomPi(m2.pg1.a, "hereisoneagene") <
##                 smAnomPi(m2.pg1.c, "nnhsisthecgene"))
##     expect_true(smAnomPi(m2.pg1.a, "hereisoneagene") <
##                 smAnomPi(m1.pg1.a, "hereisoneagene"))
## })
## date()






## ## Tiny diffs often missed. Again, culprits are often
## ## * Not expected: medianNClones(m1.pg1.a) < medianNClones(m1.pg2.a) isn't true.
## ## * Not expected: medianNClones(m1.pg1.c) > medianNClones(m1.pg2.c) isn't true.
## ## * Not expected: medianNClones(m2.pg1.c) > medianNClones(m2.pg2.c) isn't true.

## date()
## test_that("Num clones: Mutator and var mut rates and init and really scrambled names, mutPropgrowth", {
##     pseed <- sample(1:100, 1)
##     
##     cat("\n x4: a runif is", runif(1), "\n")
##     pops <- 100 ## the more, the better, but takes long.
##     ft <- .001
##     lni <- 30  ## 200
##     no <- 1e7
##     ni <- c(0, 0, 0, rep(0, lni))
##     ## scramble around names
##     names(ni) <- c("hereisoneagene",
##                    "oreoisasabgene",
##                    "nnhsisthecgene",
##                    replicate(lni,
##                              paste(sample(letters, 12), collapse = "")))
##     ni <- ni[order(names(ni))]
##     fe <- allFitnessEffects(noIntGenes = ni)
##     mutator1 <- mutator2 <- rep(1, lni + 3)
##     pg1 <- pg2 <- rep(1e-5, lni + 3)
##     ## scramble names of mutator and per-gene too
##     names(mutator1) <- sample(names(ni))
##     names(mutator2) <- sample(names(ni))
##     names(pg1) <- sample(names(ni))
##     names(pg2) <- sample(names(ni))
##     mutator1["hereisoneagene"] <- 100
##     mutator2["hereisoneagene"] <- 1
##     mutator1["oreoisasabgene"] <- 1
##     mutator2["oreoisasabgene"] <- 0.01
##     mutator1["nnhsisthecgene"] <- 0.5
##     mutator2["nnhsisthecgene"] <- 100
##     m1 <- allMutatorEffects(noIntGenes = mutator1)
##     m2 <- allMutatorEffects(noIntGenes = mutator2)
##     pg1["hereisoneagene"] <- 1e-3
##     pg2["hereisoneagene"] <- 1e-14
##     pg1["oreoisasabgene"] <- 1e-7
##     pg2["oreoisasabgene"] <- 1e-7
##     pg1["nnhsisthecgene"] <- 1e-14
##     pg2["nnhsisthecgene"] <- 1e-3
##     m1.pg1.a <- oncoSimulPop(pops,
##                            fe,
##                            mu = pg1,
##                            muEF = m1,
##                            finalTime = ft,
##                            mutationPropGrowth = TRUE,
##                            initSize = no,
##                            initMutant = "hereisoneagene",
##                            onlyCancer = FALSE, detectionProb = NA, seed = NULL, mc.cores = 2)
##     m1.pg2.a <- oncoSimulPop(pops,
##                            fe,
##                            mu = pg2,
##                            muEF = m1,
##                            finalTime = ft,
##                            mutationPropGrowth = TRUE,
##                            initSize = no,
##                            initMutant = "hereisoneagene",
##                            onlyCancer = FALSE, detectionProb = NA, seed = NULL, mc.cores = 2)
##     m1.pg1.c <- oncoSimulPop(pops,
##                            fe,
##                            mu = pg1,
##                            muEF = m1,
##                            finalTime = ft,
##                            mutationPropGrowth = TRUE,
##                            initSize = no,
##                            initMutant = "nnhsisthecgene",
##                            onlyCancer = FALSE, detectionProb = NA, seed = NULL, mc.cores = 2)
##     m1.pg2.c <- oncoSimulPop(pops,
##                            fe,
##                            mu = pg2,
##                            muEF = m1,
##                            finalTime = ft,
##                            mutationPropGrowth = TRUE,
##                            initSize = no,
##                            initMutant = "nnhsisthecgene",
##                            onlyCancer = FALSE, detectionProb = NA, seed = NULL, mc.cores = 2)
##     m1.pg1.b <- oncoSimulPop(pops,
##                            fe,
##                            mu = pg1,
##                            muEF = m1,
##                            finalTime = ft,
##                            mutationPropGrowth = TRUE,
##                            initSize = no,
##                            initMutant = "oreoisasabgene",
##                            onlyCancer = FALSE, detectionProb = NA, seed = NULL, mc.cores = 2)
##     expect_true(medianNClones(m1.pg1.a) < medianNClones(m1.pg2.a))
##     ## next can fail sometimes. Similar to why
##     ## expect_true(medianNClones(m2.pg1.a) < medianNClones(m2.pg2.a)) can fail
##     ## but here, differences might even be smaller. So commented out.
##     ## expect_true(medianNClones(m1.pg1.c) > medianNClones(m1.pg2.c))
##     expect_true(medianNClones(m1.pg1.a) > medianNClones(m1.pg1.c))
##     expect_true(medianNClones(m1.pg2.a) > medianNClones(m1.pg2.c))
##     expect_true(medianNClones(m1.pg1.a) > medianNClones(m1.pg1.b))
##     m2.pg1.a <- oncoSimulPop(pops,
##                            fe,
##                            mu = pg1,
##                            muEF = m2,
##                            finalTime = ft,
##                            mutationPropGrowth = TRUE,
##                            initSize = no,
##                            initMutant = "hereisoneagene",
##                            onlyCancer = FALSE, detectionProb = NA, seed = NULL, mc.cores = 2)
##     m2.pg2.a <- oncoSimulPop(pops,
##                            fe,
##                            mu = pg2,
##                            muEF = m2,
##                            finalTime = ft,
##                            mutationPropGrowth = TRUE,
##                            initSize = no,
##                            initMutant = "hereisoneagene",
##                            onlyCancer = FALSE, detectionProb = NA, seed = NULL, mc.cores = 2)
##     m2.pg1.c <- oncoSimulPop(pops,
##                            fe,
##                            mu = pg1,
##                            muEF = m2,
##                            finalTime = ft,
##                            mutationPropGrowth = TRUE,
##                            initSize = no,
##                            initMutant = "nnhsisthecgene",
##                            onlyCancer = FALSE, detectionProb = NA, seed = NULL, mc.cores = 2)
##     m2.pg2.c <- oncoSimulPop(pops,
##                            fe,
##                            mu = pg2,
##                            muEF = m2,
##                            finalTime = ft,
##                            mutationPropGrowth = TRUE,
##                            initSize = no,
##                            initMutant = "nnhsisthecgene",
##                            onlyCancer = FALSE, detectionProb = NA, seed = NULL, mc.cores = 2)
##     ## next test sometimes fails; it can, as we use a small difference in total mut. rate.
##     expect_true(medianNClones(m2.pg1.a) < medianNClones(m2.pg2.a))
##     expect_true(medianNClones(m2.pg1.c) > medianNClones(m2.pg2.c))
##     expect_true(medianNClones(m2.pg1.a) < medianNClones(m2.pg1.c))
##     expect_true(medianNClones(m2.pg2.a) < medianNClones(m2.pg2.c))
##     expect_true(medianNClones(m1.pg1.c) < medianNClones(m2.pg1.c))
##     expect_true(medianNClones(m1.pg1.c) < medianNClones(m2.pg2.c))
##     expect_true(medianNClones(m2.pg1.a) < medianNClones(m1.pg1.a))
##     expect_true(medianNClones(m2.pg1.a) < medianNClones(m1.pg2.a))
##     expect_true(medianNClones(m2.pg2.a) < medianNClones(m2.pg1.c))
##     expect_true(medianNClones(m1.pg2.a) > medianNClones(m1.pg1.c))
##     expect_true(smAnomPi(m2.pg1.a, "hereisoneagene") <
##                 smAnomPi(m2.pg1.c, "nnhsisthecgene"))
##     expect_true(smAnomPi(m2.pg1.a, "hereisoneagene") <
##                 smAnomPi(m1.pg1.a, "hereisoneagene"))
## })
## date()

## date()
## test_that("McFL: Num clones: Mutator and var mut rates and init and really scrambled names, mutPropGrowt", {
##     pseed <- sample(1:100, 1)
##     
##     cat("\n x5: a runif is", runif(1), "\n")
##     pops <- 200
##     ft <- .051
##     lni <- 200
##     no <- 1e5
##     ni <- c(0, 0, 0, rep(0, lni))
##     ## scramble around names
##     names(ni) <- c("hereisoneagene",
##                    "oreoisasabgene",
##                    "nnhsisthecgene",
##                    replicate(lni,
##                              paste(sample(letters, 12), collapse = "")))
##     ni <- ni[order(names(ni))]
##     fe <- allFitnessEffects(noIntGenes = ni)
##     mutator1 <- mutator2 <- rep(1, lni + 3)
##     pg1 <- pg2 <- rep(1e-6, lni + 3)
##     ## scramble names of mutator and per-gene too
##     names(mutator1) <- sample(names(ni))
##     names(mutator2) <- sample(names(ni))
##     names(pg1) <- sample(names(ni))
##     names(pg2) <- sample(names(ni))
##     mutator1["hereisoneagene"] <- 100
##     mutator2["hereisoneagene"] <- 1
##     mutator1["oreoisasabgene"] <- 1
##     mutator2["oreoisasabgene"] <- 0.01
##     mutator1["nnhsisthecgene"] <- 0.01
##     mutator2["nnhsisthecgene"] <- 100
##     m1 <- allMutatorEffects(noIntGenes = mutator1)
##     m2 <- allMutatorEffects(noIntGenes = mutator2)
##     pg1["hereisoneagene"] <- 1e-3
##     pg2["hereisoneagene"] <- 1e-14
##     pg1["oreoisasabgene"] <- 1e-7
##     pg2["oreoisasabgene"] <- 1e-7
##     pg1["nnhsisthecgene"] <- 1e-14
##     pg2["nnhsisthecgene"] <- 1e-3
##     m1.pg1.a <- oncoSimulPop(pops,
##                            fe,
##                            mu = pg1,
##                            muEF = m1,
##                            finalTime = ft, model = "McFL",
##                            mutationPropGrowth = TRUE,
##                            initSize = no,
##                            initMutant = "hereisoneagene",
##                            onlyCancer = FALSE, detectionProb = NA, seed = NULL, mc.cores = 2)
##     m1.pg2.a <- oncoSimulPop(pops,
##                            fe,
##                            mu = pg2,
##                            muEF = m1,
##                            finalTime = ft, model = "McFL",
##                            mutationPropGrowth = TRUE,
##                            initSize = no,
##                            initMutant = "hereisoneagene",
##                            onlyCancer = FALSE, detectionProb = NA, seed = NULL, mc.cores = 2)
##     m1.pg1.c <- oncoSimulPop(pops,
##                            fe,
##                            mu = pg1,
##                            muEF = m1,
##                            finalTime = ft, model = "McFL",
##                            mutationPropGrowth = TRUE,
##                            initSize = no,
##                            initMutant = "nnhsisthecgene",
##                            onlyCancer = FALSE, detectionProb = NA, seed = NULL, mc.cores = 2)
##     m1.pg2.c <- oncoSimulPop(pops,
##                            fe,
##                            mu = pg2,
##                            muEF = m1,
##                            finalTime = ft, model = "McFL",
##                            mutationPropGrowth = TRUE,
##                            initSize = no,
##                            initMutant = "nnhsisthecgene",
##                            onlyCancer = FALSE, detectionProb = NA, seed = NULL, mc.cores = 2)
##     m1.pg1.b <- oncoSimulPop(pops,
##                            fe,
##                            mu = pg1,
##                            muEF = m1,
##                            finalTime = ft, model = "McFL",
##                            mutationPropGrowth = TRUE,
##                            initSize = no,
##                            initMutant = "oreoisasabgene",
##                            onlyCancer = FALSE, detectionProb = NA, seed = NULL, mc.cores = 2)
##     expect_true(medianNClones(m1.pg1.a) < medianNClones(m1.pg2.a))
##     ## Again, too tiny diffs to be detectable with reasonable time
##     ## expect_true(medianNClones(m1.pg1.c) > medianNClones(m1.pg2.c))
##     expect_true(medianNClones(m1.pg1.a) > medianNClones(m1.pg1.c))
##     expect_true(medianNClones(m1.pg2.a) > medianNClones(m1.pg2.c))
##     expect_true(medianNClones(m1.pg1.a) > medianNClones(m1.pg1.b))
##     m2.pg1.a <- oncoSimulPop(pops,
##                            fe,
##                            mu = pg1,
##                            muEF = m2,
##                            finalTime = ft, model = "McFL",
##                            mutationPropGrowth = TRUE,
##                            initSize = no,
##                            initMutant = "hereisoneagene",
##                            onlyCancer = FALSE, detectionProb = NA, seed = NULL, mc.cores = 2)
##     m2.pg2.a <- oncoSimulPop(pops,
##                            fe,
##                            mu = pg2,
##                            muEF = m2,
##                            finalTime = ft, model = "McFL",
##                            mutationPropGrowth = TRUE,
##                            initSize = no,
##                            initMutant = "hereisoneagene",
##                            onlyCancer = FALSE, detectionProb = NA, seed = NULL, mc.cores = 2)
##     m2.pg1.c <- oncoSimulPop(pops,
##                            fe,
##                            mu = pg1,
##                            muEF = m2,
##                            finalTime = ft, model = "McFL",
##                            mutationPropGrowth = TRUE,
##                            initSize = no,
##                            initMutant = "nnhsisthecgene",
##                            onlyCancer = FALSE, detectionProb = NA, seed = NULL, mc.cores = 2)
##     m2.pg2.c <- oncoSimulPop(pops,
##                            fe,
##                            mu = pg2,
##                            muEF = m2,
##                            finalTime = ft, model = "McFL",
##                            mutationPropGrowth = TRUE,
##                            initSize = no,
##                            initMutant = "nnhsisthecgene",
##                            onlyCancer = FALSE, detectionProb = NA, seed = NULL, mc.cores = 2)
##     expect_true(medianNClones(m2.pg1.a) < medianNClones(m2.pg2.a))
##     expect_true(medianNClones(m2.pg1.c) > medianNClones(m2.pg2.c))
##     expect_true(medianNClones(m2.pg1.a) < medianNClones(m2.pg1.c))
##     expect_true(medianNClones(m2.pg2.a) < medianNClones(m2.pg2.c))
##     expect_true(medianNClones(m1.pg1.c) < medianNClones(m2.pg1.c))
##     expect_true(medianNClones(m1.pg1.c) < medianNClones(m2.pg2.c))
##     expect_true(medianNClones(m2.pg1.a) < medianNClones(m1.pg1.a))
##     expect_true(medianNClones(m2.pg1.a) < medianNClones(m1.pg2.a))
##     expect_true(medianNClones(m2.pg2.a) < medianNClones(m2.pg1.c))
##     expect_true(medianNClones(m1.pg2.a) > medianNClones(m1.pg1.c))
##     ## some of these differences are actually tiny, as the largest clone
##     ## is the initial one, with a single mutation, which is the same for
##     ## both in many cases
##     ## expect_true(smA(m2.pg1.a) < smA(m2.pg1.c))
##     ## expect_true(smA(m2.pg1.a) < smA(m1.pg1.a))
##     expect_true(smAnomPi(m2.pg1.a, "hereisoneagene") <
##                 smAnomPi(m2.pg1.c, "nnhsisthecgene"))
##     expect_true(smAnomPi(m2.pg1.a, "hereisoneagene") <
##                 smAnomPi(m1.pg1.a, "hereisoneagene"))
## })
## date()


## ## This is probably getting silly. Too many of the same thing.
## ## And done more cleanly in "Init with different mutators"
## test_that("Yet another test of ordering, add scrambling", {
##     pseed <- sample(6789:26789, 1)
##     
##     cat("\n u5: a runif is", runif(1), "\n")
##     pops <- 10
##     ft <- 1
##     lni <- 100 
##     no <- 1e5
##     ni <- c(0.1, 0.01, 0.1, rep(0, lni))
##     ## scramble around names
##     names(ni) <- c("hereisoneagene",
##                    "oreoisasabgene",
##                    "nnhsisthecgene",
##                    replicate(lni,
##                              paste(sample(letters, 12), collapse = "")))
##     ni <- ni[order(names(ni))]
##     fe <- allFitnessEffects(noIntGenes = ni)
##     pg1 <- runif(lni + 3, min = 1e-8, max = 1e-6) 
##     names(pg1) <- sample(names(ni))
##     m1 <- allMutatorEffects(noIntGenes = c("oreoisasabgene" = 5))
##     m1.pg1.b <- oncoSimulPop(pops,
##                            fe,
##                            mu = pg1,
##                            muEF = m1,
##                            finalTime = ft,
##                            mutationPropGrowth = FALSE,
##                            initSize = no,
##                            initMutant ="oreoisasabgene",
##                            onlyCancer = FALSE, detectionProb = NA, seed = NULL, mc.cores = 2)
##     expect_true(sm("oreoisasabgene", m1.pg1.b) == totalind(m1.pg1.b))
##     m1.pg1.b
##     m2 <- allMutatorEffects(noIntGenes = c("oreoisasabgene" = 20))
##     m2.pg1.b <- oncoSimulPop(pops,
##                            fe,
##                            mu = pg1,
##                            muEF = m2,
##                            finalTime = ft,
##                            mutationPropGrowth = FALSE,
##                            initSize = no,
##                            initMutant ="oreoisasabgene",
##                            onlyCancer = FALSE, detectionProb = NA, seed = NULL, mc.cores = 2)
##     expect_true(sm("oreoisasabgene", m2.pg1.b) == totalind(m2.pg1.b))
##     m2.pg1.b
##     expect_true(medianNClones(m2.pg1.b) > medianNClones(m1.pg1.b))
##     ## the next is swamped again by the initMutant gene
##     ## expect_true(smA(m2.pg1.b) > smA(m1.pg1.b))
##     expect_true(smAnomPi(m2.pg1.b, "oreoisasabgene") >
##                 smAnomPi(m1.pg1.b, "oreoisasabgene"))
## })
## date()



## ## Too convoluted, and sampling at end covered by oncoSimulSample
## date() 
## test_that("Mutator, several modules differences, sample at end only, and smaller effect", {
##     ## we also change the effect (smaller) and finalTime (longer)
##     ## you can increase ln, but then occasionally will take very long.
##     ## or increase my, or increase m1, etc.
##     pseed <- sample(99999999, 1)
##     
##     cat("\n l-mmd2: a runif is", runif(1), "\n")
##     reps <- 10
##     no <- 5e3
##     ft <- 1000 ## you need it large enough to get enough hits
##     mu <- 1e-5
##     ln <- 20 
##     m1 <- 2 ## if this is too large, easy to get it to blow.
##     ni <- rep(0, 3 * ln)
##     gna <- paste0("a", 1:ln)
##     gnb <- paste0("b", 1:ln)
##     gnc <- paste0("c", 1:ln)
##     names(ni) <- c(gna, gnb, gnc)
##     gn1 <- paste(c(gna, gnb, gnc), collapse = ", ")
##     gna <- paste(gna, collapse = ", ")
##     gnb <- paste(gnb, collapse = ", ")
##     gnc <- paste(gnc, collapse = ", ")
##     mut1 <- allMutatorEffects(epistasis = c("A" = m1),
##                               geneToModule = c("A" = gn1))
##     mut2 <- allMutatorEffects(epistasis = c("A" = m1,
##                                             "B" = m1,
##                                             "C" = m1),
##                               geneToModule = c("A" = gna,
##                                                "B" = gnb,
##                                                "C" = gnc))
##     f1 <- allFitnessEffects(noIntGenes = ni)
##     b1 <- oncoSimulPop(reps,
##                        f1,
##                        mu = mu,
##                        muEF = mut1,
##                        onlyCancer = FALSE, detectionProb = NA,
##                        initSize = no,
##                        finalTime = ft,
##                        seed = NULL, keepEvery = ft
##                        )
##     gc()
##     b2 <- oncoSimulPop(reps,
##                        f1,
##                        mu = mu,
##                        muEF = mut2,
##                        onlyCancer = FALSE, detectionProb = NA,
##                        initSize = no,
##                        finalTime = ft,
##                        seed = NULL, keepEvery = ft
##                        )
##     gc()
##     summary(b2)[, c(1:3, 8:9)]
##     summary(b1)[, c(1:3, 8:9)]
##     mean(mutsPerClone(b2))
##     mean(mutsPerClone(b1))
##     ## This is, of course, affected by sampling only at end: we do not see
##     ## the many intermediate events.
##     expect_true( median(summary(b2)$NumClones) >
##                  median(summary(b1)$NumClones))
##     expect_true( mean(mutsPerClone(b2)) >
##                  mean(mutsPerClone(b1)))
## })
## date()

