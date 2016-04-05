cat(paste("\n Starting test.mutator-long.R test at", date()))
RNGkind("L'Ecuyer-CMRG") ## for the mclapplies


date()

RNGkind("L'Ecuyer-CMRG") ## for the mclapplies

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

date()
test_that("MCFL Relative ordering of number of clones with init mutant of mutator effects", {
    ## Here we stop on finalTime, not popSize
    ## Can occasionally blow up with pE.f: pE not finite.
    pseed <-sample(9999999, 1)
    set.seed(pseed)
    cat("\n mcx2bc: the seed is", pseed, "\n")
    pops <- 10
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
                        onlyCancer = FALSE, mc.cores = 2, model = "McFL")
    ncb <- oncoSimulPop(pops, fe, muEF = fm6, finalTime =50,
                        mutationPropGrowth = FALSE,
                        initSize = 1e4,
                        initMutant = "b",
                        sampleEvery = 0.01,
                        keepEvery = 5,
                        onlyCancer = FALSE, mc.cores = 2, model = "McFL")
    ncc <- oncoSimulPop(pops, fe, muEF = fm6, finalTime =50,
                        mutationPropGrowth = FALSE,
                        initSize = 1e4,
                        initMutant = "c",
                        sampleEvery = 0.01,
                        keepEvery = 5,
                        onlyCancer = FALSE, mc.cores = 2, model = "McFL")
    ncd <- oncoSimulPop(pops, fe, muEF = fm6, finalTime =50,
                        mutationPropGrowth = FALSE,
                        initSize = 1e4,
                        initMutant = "d",
                        sampleEvery = 0.01,
                        keepEvery = 5,
                        onlyCancer = FALSE, mc.cores = 2, model = "McFL")
    expect_true( wilcox.test(summary(nca)$NumClones, 
                 summary(ncb)$NumClones, alternative = "less")$p.value < p.fail)
    expect_true(wilcox.test(summary(ncb)$NumClones ,
                summary(ncc)$NumClones, alternative = "less")$p.value < p.fail)
    expect_true( wilcox.test(summary(ncc)$NumClones,
                 summary(ncd)$NumClones, alternative = "less")$p.value < p.fail)
    summary(nca)[, c(1:3, 8:9)]
    summary(ncb)[, c(1:3, 8:9)]
    summary(ncc)[, c(1:3, 8:9)]
    summary(ncd)[, c(1:3, 8:9)]
    expect_true(t.test(mutsPerClone(nca), mutsPerClone(ncb), alternative = "less")$p.value < p.fail)
    expect_true(t.test(mutsPerClone(ncb), mutsPerClone(ncc), alternative = "less")$p.value < p.fail)
    expect_true(t.test(mutsPerClone(ncc), mutsPerClone(ncd), alternative = "less")$p.value < p.fail)
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
    m1.pg1.a <- oncoSimulPop(pops,
                           fe,
                           mu = pg1,
                           muEF = m1,
                           finalTime = ft,
                           mutationPropGrowth = FALSE,
                           initSize = no,
                           initMutant = "hereisoneagene",
                           sampleEvery = 0.01, keepEvery = 5, seed = NULL,
                           onlyCancer = FALSE, mc.cores = 2)
    m1.pg1.b <- oncoSimulPop(pops,
                             fe,
                             mu = pg1,
                             muEF = m1,
                             finalTime = ft,
                             mutationPropGrowth = FALSE,
                             initSize = no,
                             initMutant = "oreoisasabgene",
                             sampleEvery = 0.01, keepEvery = 5, seed = NULL,
                             onlyCancer = FALSE, mc.cores = 2)
    m1.pg1.c <- oncoSimulPop(pops,
                           fe,
                           mu = pg1,
                           muEF = m1,
                           finalTime = ft,
                           mutationPropGrowth = FALSE,
                           initSize = no,
                           initMutant = "nnhsisthecgene",
                           sampleEvery = 0.01, keepEvery = 5, seed = NULL,
                           onlyCancer = FALSE, mc.cores = 2)
    expect_true(wilcox.test(NClones(m1.pg1.a),  NClones(m1.pg1.b), alternative = "greater")$p.value < p.fail)
    expect_true(wilcox.test(NClones(m1.pg1.b),  NClones(m1.pg1.c), alternative = "greater")$p.value < p.fail)
    expect_true(mean(mutsPerClone(m1.pg1.a)) > mean(mutsPerClone(m1.pg1.b)))
    expect_true(mean(mutsPerClone(m1.pg1.b)) > mean(mutsPerClone(m1.pg1.c)))
    expect_true(smAnomPi(m1.pg1.a, "hereisoneagene") >
                smAnomPi(m1.pg1.b, "oreoisasabgene"))
    expect_true(smAnomPi(m1.pg1.b, "oreoisasabgene") >
                smAnomPi(m1.pg1.c, "nnhsisthecgene"))
    summary(m1.pg1.a)[, c(1:3, 8:9)]
    summary(m1.pg1.b)[, c(1:3, 8:9)]
    summary(m1.pg1.c)[, c(1:3, 8:9)]
})
date()


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
    pops <- 6
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
                        onlyCancer = FALSE, mc.cores = 2)
    pop100 <- oncoSimulPop(pops,
                        fe,
                        mu = muvector,
                        muEF = m100,
                        finalTime = ft,
                        mutationPropGrowth = FALSE,
                        initSize = no,
                        initMutant = names(mutator10),
                        onlyCancer = FALSE, mc.cores = 2)
    ## number of total mutations
    expect_true(smAnomPi(pop10, names(mutator10)) < smAnomPi(pop100, names(mutator100)))
    ## number of clones
    expect_true(wilcox.test(NClones(pop100),  NClones(pop10), alternative = "greater")$p.value < p.fail)
    expect_true(mean(mutsPerClone(pop10)) < mean(mutsPerClone(pop100)))
    
})
date()

date()
test_that("Same mu vector, different mutator; diffs in number muts, larger t", {
    
    ## reproduction, death, and double and possibly triple mutants. We
    ## decrease init pop size to make this fast.
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n nm1: the seed is", pseed, "\n")
    pops <- 6
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
                        onlyCancer = FALSE, mc.cores = 2)
    pop100 <- oncoSimulPop(pops,
                        fe,
                        mu = muvector,
                        muEF = m100,
                        finalTime = ft,
                        mutationPropGrowth = FALSE,
                        initSize = no,
                        initMutant = names(mutator10),
                        onlyCancer = FALSE, mc.cores = 2)
    ## number of total mutations
    expect_true(smAnomPi(pop10, names(mutator10)) < smAnomPi(pop100, names(mutator100)))
    ## number of clones
    expect_true(wilcox.test(NClones(pop100),  NClones(pop10), alternative = "greater")$p.value < p.fail)
    expect_true(mean(mutsPerClone(pop10)) < mean(mutsPerClone(pop100)))

})


date()
test_that("Relative ordering of number of clones with mut prop growth and init and scrambled names", {
    ## Can occasionally blow up with pE.f: pE not finite.
    pseed <- sample(99999999, 1)
    set.seed(pseed)
    cat("\n x2ef: the seed is", pseed, "\n")
    pops <- 30
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
                        initSize = no,
                        initMutant = "thisistheagene",
                        onlyCancer = FALSE, mc.cores = 2)
    mnpg <- oncoSimulPop(pops, fe, muEF = fm1,
                         finalTime = ft,
                         mutationPropGrowth = FALSE,
                         initSize = no,
                         initMutant = "thisistheagene",
                         onlyCancer = FALSE, mc.cores = 2)
    pg <- oncoSimulPop(pops, fe, 
                       finalTime = ft,
                       mutationPropGrowth = TRUE,
                       initSize = no,
                       initMutant = "thisistheagene",
                       onlyCancer = FALSE, mc.cores = 2)
    npg <- oncoSimulPop(pops, fe, 
                        finalTime = ft,
                        mutationPropGrowth = FALSE,
                        initSize = no,
                        initMutant = "thisistheagene",
                        onlyCancer = FALSE, mc.cores = 2)
      ## I once saw a weird thing
    expect_true(var(summary(mpg)$NumClones) > 1e-4)
    expect_true(var(summary(mnpg)$NumClones) > 1e-4)
    expect_true(var(summary(pg)$NumClones) > 1e-4)
    expect_true(var(summary(npg)$NumClones) > 1e-4)
    ## These are the real tests
    expect_true( wilcox.test(summary(mpg)$NumClones, 
                 summary(mnpg)$NumClones, alternative = "greater")$p.value < p.fail)
    expect_true(wilcox.test(summary(mpg)$NumClones, 
                summary(pg)$NumClones), alternative = "greater")$p.value < p.fail)
    expect_true( wilcox.test(summary(mnpg)$NumClones, 
                 summary(npg)$NumClones), alternative = "greater")$p.value < p.fail)
    expect_true( wilcox.test(summary(pg)$NumClones, 
                 summary(npg)$NumClones), alternative = "greater")$p.value < p.fail)
    expect_true(mean(mutsPerClone(mpg)) > mean(mutsPerClone(mnpg)))
    expect_true(mean(mutsPerClone(mpg)) > mean(mutsPerClone(pg)))
    expect_true(mean(mutsPerClone(mnpg)) > mean(mutsPerClone(npg)))
    expect_true(mean(mutsPerClone(pg)) > mean(mutsPerClone(npg)))
})
date()


date()
test_that("Expect freq genotypes, mutator and var mut rates", {
    ## increase mutator, decrease max mu
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n u8_1: the seed is", pseed, "\n")
    pops <- 100
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
                           onlyCancer = FALSE, seed = NULL, mc.cores = 2)
    expect_true(sm("oreoisasabgene", m1.pg1.b) == totalind(m1.pg1.b))
    ## If you want to see the numbers
    ## enom("oreoisasabgene", pg1)
    ## snom("oreoisasabgene", m1.pg1.b)
    p.fail <- 1e-3
    expect_true(chisq.test(snom("oreoisasabgene", m1.pg1.b),
                           p = pnom("oreoisasabgene", pg1, no, pops))$p.value > p.fail)
})
date()


date()
test_that("Relative ordering of number of clones with init mutant of mutator effects and s = 0", {
    ## Can occasionally blow up with pE.f: pE not finite.
    pseed <-sample(9999999, 1)
    set.seed(pseed)
    cat("\n x2cd: the seed is", pseed, "\n")
    pops <- 30
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
                        onlyCancer = FALSE, mc.cores = 2)
    ncb <- oncoSimulPop(pops, fe, muEF = fm6, finalTime =50,
                        mutationPropGrowth = FALSE,
                        initSize = 1e4,
                        initMutant = "b",
                        sampleEvery = 0.01,
                        keepEvery = 5,
                        onlyCancer = FALSE, mc.cores = 2)
    ncc <- oncoSimulPop(pops, fe, muEF = fm6, finalTime =50,
                        mutationPropGrowth = FALSE,
                        initSize = 1e4,
                        initMutant = "c",
                        sampleEvery = 0.01,
                        keepEvery = 5,
                        onlyCancer = FALSE, mc.cores = 2)
    ncd <- oncoSimulPop(pops, fe, muEF = fm6, finalTime =50,
                        mutationPropGrowth = FALSE,
                        initSize = 1e4,
                        initMutant = "d",
                        sampleEvery = 0.01,
                        keepEvery = 5,
                        onlyCancer = FALSE, mc.cores = 2)
    ## I once saw a weird thing
    expect_true(var(summary(nca)$NumClones) > 1e-4)
    expect_true(var(summary(ncb)$NumClones) > 1e-4)
    expect_true(var(summary(ncc)$NumClones) > 1e-4)
    expect_true(var(summary(ncd)$NumClones) > 1e-4)
    ## These are the real tests
    expect_true( wilcox.test(summary(nca)$NumClones,
                 summary(ncb)$NumClones, alternative = "less")$p.value < p.fail)
    expect_true(wilcox.test(summary(ncb)$NumClones,
                summary(ncc)$NumClones, alternative = "less")$p.value < p.fail)
    expect_true( wilcox.test(summary(ncc)$NumClones,
                 summary(ncd)$NumClones, alternative = "less")$p.value < p.fail)
    expect_true(mean(mutsPerClone(nca)) < mean(mutsPerClone(ncb)))
    expect_true(mean(mutsPerClone(ncb)) < mean(mutsPerClone(ncc)))
    expect_true(mean(mutsPerClone(ncc)) < mean(mutsPerClone(ncd)))
})
date()



date()
test_that("Relative ordering of number of clones with init mutant of mutator effects", {
    ## Here we stop on finalTime, not popSize
    ## Can occasionally blow up with pE.f: pE not finite.
    pseed <-sample(9999999, 1)
    set.seed(pseed)
    cat("\n x2bc: the seed is", pseed, "\n")
    pops <- 40
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
                        onlyCancer = FALSE, mc.cores = 2)
    ncb <- oncoSimulPop(pops, fe, muEF = fm6, finalTime =50,
                        mutationPropGrowth = FALSE,
                        initSize = 1e4,
                        initMutant = "b",
                        sampleEvery = 0.01,
                        keepEvery = 5,
                        onlyCancer = FALSE, mc.cores = 2)
    ncc <- oncoSimulPop(pops, fe, muEF = fm6, finalTime =50,
                        mutationPropGrowth = FALSE,
                        initSize = 1e4,
                        initMutant = "c",
                        sampleEvery = 0.01,
                        keepEvery = 5,
                        onlyCancer = FALSE, mc.cores = 2)
    ncd <- oncoSimulPop(pops, fe, muEF = fm6, finalTime =50,
                        mutationPropGrowth = FALSE,
                        initSize = 1e4,
                        initMutant = "d",
                        sampleEvery = 0.01,
                        keepEvery = 5,
                        onlyCancer = FALSE, mc.cores = 2)
    expect_true( wilcox.test(summary(nca)$NumClones,
                 summary(ncb)$NumClones, alternative = "less")$p.value < p.fail)
    expect_true(wilcox.test(summary(ncb)$NumClones,
                summary(ncc)$NumClones, alternative = "less")$p.value < p.fail)
    expect_true( wilcox.test(summary(ncc)$NumClones,
                 summary(ncd)$NumClones, alternative = "less")$p.value < p.fail)
    expect_true(mean(mutsPerClone(nca)) < mean(mutsPerClone(ncb)))
    expect_true(mean(mutsPerClone(ncb)) < mean(mutsPerClone(ncc)))
    expect_true(mean(mutsPerClone(ncc)) < mean(mutsPerClone(ncd)))
    summary(nca)[, c(1:3, 8:9)]
    summary(ncb)[, c(1:3, 8:9)]
    summary(ncc)[, c(1:3, 8:9)]
    summary(ncd)[, c(1:3, 8:9)]
})
date()

date()
test_that("MCFL Relative ordering of number of clones with mutator effects", {
    ## Can occasionally blow up with pE.f: pE not finite.
    pseed <-sample(9999999, 1)
    set.seed(pseed)
    cat("\n mcx1: the seed is", pseed, "\n")
    pops <- 40
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
                        sampleEvery = 0.01,
                        keepEvery = 5,
                        initSize = 1e6, mc.cores = 2, model = "McFL",
                        onlyCancer = FALSE)
    fm8 <- allMutatorEffects(noIntGenes = c("a" = 1,
                                            "b" = 1,
                                            "c" = 1,
                                            "d" = 1))
    nc2 <- oncoSimulPop(pops, fe, muEF = fm8, finalTime =100,
                        mutationPropGrowth = FALSE,
                        sampleEvery = 0.01,
                        keepEvery = 5,
                        initSize = 1e6, mc.cores = 2, model = "McFL",
                        onlyCancer = FALSE)
    fm7 <- allMutatorEffects(noIntGenes = c("a" = 1e-6,
                                            "b" = 1e-6,
                                            "c" = 1e-6,
                                            "d" = 1e-6))
    nc3 <- oncoSimulPop(pops, fe, muEF = fm7, finalTime =100,
                        mutationPropGrowth = FALSE,
                        sampleEvery = 0.01,
                        keepEvery = 5,
                        initSize = 1e6, mc.cores = 2, model = "McFL",
                        onlyCancer = FALSE)
    expect_true(wilcox.test(summary(nc1)$NumClones) > median(summary(nc2)$NumClones))
    expect_true(wilcox.test(summary(nc2)$NumClones) > median(summary(nc3)$NumClones))
    summary(nc1)[, c(1:3, 8:9)]
    summary(nc2)[, c(1:3, 8:9)]
    summary(nc3)[, c(1:3, 8:9)]
    expect_true(mean(mutsPerClone(nc1)) > mean(mutsPerClone(nc2)))
    expect_true(mean(mutsPerClone(nc2)) > mean(mutsPerClone(nc3)))
})
date()

date()
test_that("MCFL Relative ordering of number of clones with init mutant of mutator effects and s = 0", {
    ## Can occasionally blow up with pE.f: pE not finite.
    pseed <-sample(9999999, 1)
    set.seed(pseed)
    cat("\n mcx2cd: the seed is", pseed, "\n")
    pops <- 40
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
                        onlyCancer = FALSE, mc.cores = 2, model = "McFL")
    ncb <- oncoSimulPop(pops, fe, muEF = fm6, finalTime =50,
                        mutationPropGrowth = FALSE,
                        initSize = 1e4,
                        initMutant = "b",
                        sampleEvery = 0.01,
                        keepEvery = 5,
                        onlyCancer = FALSE, mc.cores = 2, model = "McFL")
    ncc <- oncoSimulPop(pops, fe, muEF = fm6, finalTime =50,
                        mutationPropGrowth = FALSE,
                        initSize = 1e4,
                        initMutant = "c",
                        sampleEvery = 0.01,
                        keepEvery = 5,
                        onlyCancer = FALSE, mc.cores = 2, model = "McFL")
    ncd <- oncoSimulPop(pops, fe, muEF = fm6, finalTime =50,
                        mutationPropGrowth = FALSE,
                        initSize = 1e4,
                        initMutant = "d",
                        sampleEvery = 0.01,
                        keepEvery = 5,
                        onlyCancer = FALSE, mc.cores = 2, model = "McFL")
    ## I once saw a weird thing
    expect_true(var(summary(nca)$NumClones) > 1e-4)
    expect_true(var(summary(ncb)$NumClones) > 1e-4)
    expect_true(var(summary(ncc)$NumClones) > 1e-4)
    expect_true(var(summary(ncd)$NumClones) > 1e-4)
    ## These are the real tests
    expect_true( wilcox.test(summary(nca)$NumClones,
                 summary(ncb)$NumClones, alternative = "less")$p.value < p.fail)
    expect_true(wilcox.test(summary(ncb)$NumClones,
                summary(ncc)$NumClones, alternative = "less")$p.value < p.fail)
    expect_true( wilcox.test(summary(ncc)$NumClones,
                 summary(ncd)$NumClones, alternative = "less")$p.value < p.fail)
    expect_true(mean(mutsPerClone(nca)) < mean(mutsPerClone(ncb)))
    expect_true(mean(mutsPerClone(ncb)) < mean(mutsPerClone(ncc)))
    expect_true(mean(mutsPerClone(ncc)) < mean(mutsPerClone(ncd)))
})
date()


## FIXME move later to long
## Slow (~ 5 seconds) but tests modules of mutator nicely.
date() ## Beware: this uses a lot of RAM without the gc()
test_that("Mutator modules differences", {
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n mmd1: the seed is", pseed, "\n")
    reps <- 20
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
                       onlyCancer = FALSE,
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
                       onlyCancer = FALSE,
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
                       onlyCancer = FALSE,
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
    expect_true( wilcox.test(summary(b3)$NumClones, 
                 summary(b2)$NumClones, alternative = "greater")$p.value < p.fail)
    expect_true( wilcox.test(summary(b2)$NumClones, 
                 summary(b1)$NumClones, alternative = "greater")$p.value < p.fail)
    expect_true( mean(mutsPerClone(b3)) >
                 mean(mutsPerClone(b2)))
    expect_true( mean(mutsPerClone(b2)) >
                 mean(mutsPerClone(b1)))
})
date()




date() ## Beware: this uses a lot of RAM without the gc()
test_that("McFL: Mutator modules differences", {
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n MCFLmmd1: the seed is", pseed, "\n")
    reps <- 20
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
                       onlyCancer = FALSE,
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
                       onlyCancer = FALSE,
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
                       onlyCancer = FALSE,
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
    expect_true( wilcox.test(summary(b3)$NumClones, 
                 summary(b2)$NumClones, alternative = "greater")$p.value < p.fail)
    expect_true( wilcox.test(summary(b2)$NumClones, 
                 summary(b1)$NumClones, alternative = "greater")$p.value < p.fail)
    expect_true( mean(mutsPerClone(b3)) >
                 mean(mutsPerClone(b2)))
    expect_true( mean(mutsPerClone(b2)) >
                 mean(mutsPerClone(b1)))
})
date()


date()
test_that("Relative ordering of number of clones with mutator effects", {
    ## Can occasionally blow up with pE.f: pE not finite.
    pseed <-sample(9999999, 1)
    set.seed(pseed)
    cat("\n x1: the seed is", pseed, "\n")
    pops <- 60
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
                        sampleEvery = 0.01,
                        keepEvery = 5,
                        initSize = 1e6, mc.cores = 2,
                        onlyCancer = FALSE)
    fm8 <- allMutatorEffects(noIntGenes = c("a" = 1,
                                            "b" = 1,
                                            "c" = 1,
                                            "d" = 1))
    nc2 <- oncoSimulPop(pops, fe, muEF = fm8, finalTime =250,
                        mutationPropGrowth = FALSE,
                        sampleEvery = 0.01,
                        keepEvery = 5,
                        initSize = 1e6, mc.cores = 2,
                        onlyCancer = FALSE)
    fm7 <- allMutatorEffects(noIntGenes = c("a" = 1e-6,
                                            "b" = 1e-6,
                                            "c" = 1e-6,
                                            "d" = 1e-6))
    nc3 <- oncoSimulPop(pops, fe, muEF = fm7, finalTime =250,
                        mutationPropGrowth = FALSE,
                        sampleEvery = 0.01,
                        keepEvery = 5,
                        initSize = 1e6, mc.cores = 2,
                        onlyCancer = FALSE)
    expect_true(wilcox.test(summary(nc1)$NumClones) > median(summary(nc2)$NumClones))
    expect_true(wilcox.test(summary(nc2)$NumClones) > median(summary(nc3)$NumClones))
    expect_true(mean(mutsPerClone(nc1)) > mean(mutsPerClone(nc2)))
    expect_true(mean(mutsPerClone(nc2)) > mean(mutsPerClone(nc3)))
    summary(nc1)[, c(1:3, 8:9)]
    summary(nc2)[, c(1:3, 8:9)]
    summary(nc3)[, c(1:3, 8:9)]
})
date()


## FIXME: move to long, later, and increase reps.
## very slow, because huge number of clones. But tests several phenomena comprehensively.
## same with McFL below
date()
test_that("per-gene-mut rates and mutator", {
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n oss11: the seed is", pseed, "\n")
    ng <- 10
    ni <- rep(0, ng)
    m1 <- runif(ng, min = 1e-7, max = 5e-6)
    m2 <- runif(ng, min = 1e-5, max = 1e-4)
    names(ni) <- names(m1) <- names(m2) <- c(replicate(ng,
                                 paste(sample(letters, 12), collapse = "")))
    fe1 <- allFitnessEffects(noIntGenes = ni)
    ft <- 50
    no <- 5e5 
    reps <- 60
    gn <- paste(names(ni), collapse = ", ")
    ## MUs used to be 25 and 100. Way too slow.
    mutator1 <- allMutatorEffects(epistasis = c("MU" = 15),
                                  geneToModule = c("MU" = gn))
    mutator2 <- allMutatorEffects(epistasis = c("MU" = 30),
                                  geneToModule = c("MU" = gn))
    m1.mutator0 <- oncoSimulPop(reps,
                           fe1,
                           mu = m1,
                           onlyCancer = FALSE,
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
                           onlyCancer = FALSE,
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
                           onlyCancer = FALSE,
                           initSize = no,
                           finalTime = ft,
                           sampleEvery = 0.01,
                           keepEvery = 5,
                           seed = NULL, mc.cores = 2
                           )
    m2.mutator0 <- oncoSimulPop(reps,
                           fe1,
                           mu = m2,
                           onlyCancer = FALSE,
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
                           onlyCancer = FALSE,
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
                           onlyCancer = FALSE,
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
    expect_true( wilcox.test(summary(m1.mutator2)$NumClones, 
                 summary(m1.mutator1)$NumClones, alternative = "greater")$p.value < p.fail)
    expect_true( wilcox.test(summary(m1.mutator1)$NumClones, 
                 summary(m1.mutator0)$NumClones, alternative = "greater")$p.value < p.fail)
    expect_true( wilcox.test(summary(m2.mutator2)$NumClones, 
                 summary(m2.mutator1)$NumClones, alternative = "greater")$p.value < p.fail)
    expect_true( wilcox.test(summary(m2.mutator1)$NumClones, 
                 summary(m2.mutator0)$NumClones, alternative = "greater")$p.value < p.fail)
    expect_true( mean(mutsPerClone(m1.mutator2)) >
                 mean(mutsPerClone(m1.mutator1)))
    expect_true( mean(mutsPerClone(m1.mutator1)) >
                 mean(mutsPerClone(m1.mutator0)))
    expect_true( mean(mutsPerClone(m2.mutator2)) >
                 mean(mutsPerClone(m2.mutator1)))
    expect_true( mean(mutsPerClone(m2.mutator1)) >
                 mean(mutsPerClone(m2.mutator0)))
    ## Increases in mutation rates increase clones, etc, within levels of
    ## mutator.
    expect_true( wilcox.test(summary(m1.mutator0)$NumClones,
                 summary(m2.mutator0, alternative = "less")$p.value < p.fail)NumClones))
    expect_true( wilcox.test(summary(m1.mutator1)$NumClones,
                 summary(m2.mutator1, alternative = "less")$p.value < p.fail)NumClones))
    expect_true( wilcox.test(summary(m1.mutator2)$NumClones,
                 summary(m2.mutator2, alternative = "less")$p.value < p.fail)NumClones))
    expect_true( mean(mutsPerClone(m1.mutator0)) <
                 mean(mutsPerClone(m2.mutator0)))
    expect_true( mean(mutsPerClone(m1.mutator1)) <
                 mean(mutsPerClone(m2.mutator1)))
    expect_true( mean(mutsPerClone(m1.mutator2)) <
                 mean(mutsPerClone(m2.mutator2)))
})
date()


date()
test_that("McFL: per-gene-mut rates and mutator", {
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n mcfloss11: the seed is", pseed, "\n")
    ng <- 10
    ni <- rep(0, ng)
    m1 <- runif(ng, min = 1e-7, max = 5e-6)
    m2 <- runif(ng, min = 1e-5, max = 1e-4)
    names(ni) <- names(m1) <- names(m2) <- c(replicate(ng,
                                 paste(sample(letters, 12), collapse = "")))
    fe1 <- allFitnessEffects(noIntGenes = ni)
    ft <- 50
    no <- 5e5 
    reps <- 60
    gn <- paste(names(ni), collapse = ", ")
    mutator1 <- allMutatorEffects(epistasis = c("MU" = 15),
                                  geneToModule = c("MU" = gn))
    mutator2 <- allMutatorEffects(epistasis = c("MU" = 30),
                                  geneToModule = c("MU" = gn))
    m1.mutator0 <- oncoSimulPop(reps,
                           fe1,
                           mu = m1,
                           onlyCancer = FALSE,
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
                           onlyCancer = FALSE,
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
                           onlyCancer = FALSE,
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
                           onlyCancer = FALSE,
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
                           onlyCancer = FALSE,
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
                           onlyCancer = FALSE,
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
    expect_true( wilcox.test(summary(m1.mutator2)$NumClones, 
                 summary(m1.mutator1)$NumClones, alternative = "greater")$p.value < p.fail)
    expect_true( wilcox.test(summary(m1.mutator1)$NumClones, 
                 summary(m1.mutator0)$NumClones, alternative = "greater")$p.value < p.fail)
    expect_true( wilcox.test(summary(m2.mutator2)$NumClones, 
                 summary(m2.mutator1)$NumClones, alternative = "greater")$p.value < p.fail)
    expect_true( wilcox.test(summary(m2.mutator1)$NumClones, 
                 summary(m2.mutator0)$NumClones, alternative = "greater")$p.value < p.fail)
    expect_true( mean(mutsPerClone(m1.mutator2)) >
                 mean(mutsPerClone(m1.mutator1)))
    expect_true( mean(mutsPerClone(m1.mutator1)) >
                 mean(mutsPerClone(m1.mutator0)))
    expect_true( mean(mutsPerClone(m2.mutator2)) >
                 mean(mutsPerClone(m2.mutator1)))
    expect_true( mean(mutsPerClone(m2.mutator1)) >
                 mean(mutsPerClone(m2.mutator0)))
    ## Increases in mutation rates increase clones, etc, within levels of
    ## mutator.
    expect_true( wilcox.test(summary(m1.mutator0)$NumClones,
                 summary(m2.mutator0, alternative = "less")$p.value < p.fail)NumClones))
    expect_true( wilcox.test(summary(m1.mutator1)$NumClones,
                 summary(m2.mutator1, alternative = "less")$p.value < p.fail)NumClones))
    expect_true( wilcox.test(summary(m1.mutator2)$NumClones,
                 summary(m2.mutator2, alternative = "less")$p.value < p.fail)NumClones))
    expect_true( mean(mutsPerClone(m1.mutator0)) <
                 mean(mutsPerClone(m2.mutator0)))
    expect_true( mean(mutsPerClone(m1.mutator1)) <
                 mean(mutsPerClone(m2.mutator1)))
    expect_true( mean(mutsPerClone(m1.mutator2)) <
                 mean(mutsPerClone(m2.mutator2)))
})
date()


date() 
test_that("Mutator, several modules differences, McFL", {
    pseed <- sample(99999999, 1)
    pseed <- 91339980
    set.seed(pseed)
    cat("\n mmdSMMC1: the seed is", pseed, "\n")
    reps <- 30
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
                       onlyCancer = FALSE,
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
                       onlyCancer = FALSE,
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
    expect_true( wilcox.test(summary(b2)$NumClones, 
                 summary(b1)$NumClones, alternative = "greater")$p.value < p.fail)
    expect_true( mean(mutsPerClone(b2)) >
                 mean(mutsPerClone(b1)))
})
date()




date()
test_that("McFL: Expect freq genotypes, mutator and var mut rates", {
    ## increase mutator
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n u8: the seed is", pseed, "\n")
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
                           onlyCancer = FALSE, seed = NULL, mc.cores = 2)
    expect_true(sm("oreoisasabgene", m1.pg1.b) == totalind(m1.pg1.b))
    ## enom("oreoisasabgene", pg1)
    ## snom("oreoisasabgene", m1.pg1.b)
    p.fail <- 1e-3
    expect_true(chisq.test(snom("oreoisasabgene", m1.pg1.b),
                           p = pnom("oreoisasabgene", pg1, no, pops))$p.value > p.fail)
})
date()




date()
test_that("Expect freq genotypes, mutator and var mut rates", {
    ## Similar to above, but mutator has a single element, not the whole
    ## vector.
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n u7_1: the seed is", pseed, "\n")
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
                           onlyCancer = FALSE, seed = NULL, mc.cores = 2)
    expect_true(sm("oreoisasabgene", m1.pg1.b) == totalind(m1.pg1.b))
    enom("oreoisasabgene", pg1, no, pops)
    snom("oreoisasabgene", m1.pg1.b)
    p.fail <- 1e-3
    expect_true(chisq.test(snom("oreoisasabgene", m1.pg1.b),
                           p = pnom("oreoisasabgene", pg1, no, pops))$p.value > p.fail)
})
date()


cat(paste("\n Finished test.mutator-long.R test at", date(), "\n"))







