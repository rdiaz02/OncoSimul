### Since some tests are slow and some tests are very fragile, for now I
### leave date() and seed()

cat(paste("\n Starting test.mutator.R test at", date()))


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




test_that("This should not crash", {
    ## This used to crash because of not nulling the empty mutator effects
    fe <- allFitnessEffects(epistasis = c("a : b" = 0.3,
                                          "b : c" = 0.5),
                            noIntGenes = c("e" = 0.1),
                            drvNames = c("a", "b", "c"))
    moo <- rep(1e-5, 4)
    names(moo) <- c("a", "b", "c", "e")
    expect_output(print(oncoSimulIndiv(fe,
                                       mu = moo,
                                       finalTime = 1,
                                       mutationPropGrowth = FALSE,
                                       initSize = 1000,
                                       onlyCancer = FALSE)),
                  "Individual OncoSimul",
                  fixed = TRUE)
    muvar2 <- c("U" = 1e-5, "z" = 1e-5, "e" = 1e-5, "m" = 1e-5, "D" = 1e-5)
    ni1 <- rep(0.02, 5)
    names(ni1) <- names(muvar2)
    fe1 <- allFitnessEffects(noIntGenes = ni1)
    no <- 1e5
    reps <- 10
    bb <- oncoSimulIndiv(fe1, mu = muvar2, onlyCancer = FALSE,
                         initSize = no,
                         finalTime = 1,
                         seed = NULL
                         )
    expect_output(print(bb),
                  "Individual OncoSimul",
                  fixed = TRUE)
    expect_output(print(oncoSimulPop(4,
                                     fe1,
                                     mu = muvar2,
                                     onlyCancer = FALSE,
                                     initSize = 1000,
                                     finalTime = 1,
                                     seed = NULL, mc.cores = 2)),
                  "Population of OncoSimul",
                  fixed = TRUE)
    expect_output(print(oncoSimulPop(4,
                                     fe,
                                     mu = moo,
                                     onlyCancer = FALSE,
                                     initSize = 1000,
                                     finalTime = 1,
                                     seed = NULL, mc.cores = 2)),
                  "Population of OncoSimul",
                  fixed = TRUE)
})


test_that("eval fitness and mut OK", {
    fe <- allFitnessEffects(epistasis = c("a : b" = 0.3,
                                          "b : c" = 0.5),
                            noIntGenes = c("e" = 0.1))
    fm <- allMutatorEffects(noIntGenes = c("a" = 10,
                                           "c" = 5))
    expect_output(ou <- evalGenotypeFitAndMut("a", fe, fm, verbose = TRUE),
                  "10", fixed = TRUE)
    expect_identical(ou, c(1, 10))
    expect_identical(evalGenotypeFitAndMut("b", fe, fm),
                     c(1, 1))
    expect_identical(evalGenotypeFitAndMut("e", fe, fm),
                     c(1.1, 1))
    expect_identical(evalGenotypeFitAndMut("b, e", fe, fm),
                     c(1.1, 1))
    expect_identical(evalGenotypeFitAndMut("a, b, e", fe, fm),
                     c(1.3 * 1.1, 10))
    expect_identical(evalGenotypeFitAndMut("a, b, c, e", fe, fm),
                     c(1.3 * 1.5 * 1.1, 10 * 5))
})

test_that("expect output oncoSimulIndiv", {
    fe <- allFitnessEffects(noIntGenes = c("a" = 0.2,
                                           "c" = 0.4,
                                           "d" = 0.6,
                                           "e" = 0.1))
    fm <- allMutatorEffects(noIntGenes = c("a" = 10,
                                           "c" = 5))
    expect_output(oncoSimulIndiv(fe, muEF = fm, sampleEvery = 0.01,
                                 keepEvery = 5),
                  "Individual OncoSimul trajectory",
                  fixed = TRUE)
    expect_output(oncoSimulIndiv(fe, muEF = fm, sampleEvery = 0.01,
                                 keepEvery = 5),
                  "Individual OncoSimul trajectory",
                  fixed = TRUE)
})

test_that("eval mut genotypes", {
    fe <- allFitnessEffects(epistasis = c("a : b" = 0.3,
                                          "b : c" = 0.5),
                            noIntGenes = c("e" = 0.1),
                            drvNames = c(letters[1:3]))
    fm <- allMutatorEffects(noIntGenes = c("a" = 10,
                                           "c" = 5))
    expect_identical(evalAllGenotypesMut(fm)[, 2],
                     c(10, 5, 50))
    expect_identical(evalGenotypeMut("a", fm),
                     10)
    expect_identical(evalGenotypeMut("c", fm),
                     5)
    expect_identical(evalGenotypeMut("c, a", fm),
                     50)
    expect_identical(evalGenotypeMut("a, c", fm),
                     50)
    expect_identical(evalGenotypeMut("a > c", fm),
                     50)
    expect_identical(evalGenotypeMut("c > a", fm),
                     50)
    expect_error(evalGenotypeMut("b", fm),
                 "genotype contains NAs or genes not in fitnessEffects/mutatorEffects",
                 fixed = TRUE)
})

test_that("we evaluate the WT", {
    ## Is fitness of wildtype always 0? Really? Evaluate it.
    ## It is: see evalGenotypeFitness
    fe <- allFitnessEffects(epistasis = c("a : b" = 0.3,
                                          "b : c" = 0.5),
                            noIntGenes = c("e" = 0.1))
    expect_warning(ou <- OncoSimulR:::evalRGenotype(vector(mode = "integer",
                                                           length = 0),
                                                    fe, TRUE, FALSE,
                                                    "evalGenotype"),
                   "WARNING: you have evaluated fitness/mutator status of a genotype of length zero",
                   fixed = TRUE)
    expect_identical(ou, 1)
})


test_that("we evaluate the WT, 2", {
    fe <- allFitnessEffects(epistasis = c("a : b" = 0.3,
                                          "b : c" = 0.5),
                            noIntGenes = c("e" = 0.1))
    fm <- OncoSimulR:::allMutatorEffects(noIntGenes = c("a" = 10,
                                                        "c" = 5))
    expect_warning(ou2 <- OncoSimulR:::evalRGenotypeAndMut(
                       vector(mode = "integer", length = 0),
                       fe,
                       fm,
                       OncoSimulR:::matchGeneIDs(fm, fe)$Reduced,
                       TRUE, FALSE),
                   "WARNING: you have evaluated fitness of a genotype of length zero.",
                   fixed = TRUE)
    expect_identical(ou2, c(1, 1))
})
    


test_that("evaluating genotype and mutator", {
    fe <- allFitnessEffects(epistasis = c("a : b" = 0.3,
                                          "b : c" = 0.5),
                            noIntGenes = c("e" = 0.1),
                            drvNames = letters[1:3])
    fm <- allMutatorEffects(noIntGenes = c("a" = 10,
                                           "c" = 5))
    ou <- evalAllGenotypesFitAndMut(fe, fm, order = FALSE)
    expect_true(all(dplyr::filter(ou, Genotype == "a")[, c(2, 3)] ==
                    c(1, 10)))
    expect_true(all(dplyr::filter(ou, Genotype == "b")[, c(2, 3)] ==
                    c(1, 1)))
    expect_true(all(dplyr::filter(ou, Genotype == "c")[, c(2, 3)] ==
                    c(1, 5)))
    expect_true(all(dplyr::filter(ou, Genotype == "e")[, c(2, 3)] ==
                    c(1.1, 1)))
    expect_true(all(dplyr::filter(ou, Genotype == "a, b, c")[, c(2, 3)] ==
                    c(1.3 * 1.5, 10 * 5)))
    expect_true(all(dplyr::filter(ou, Genotype == "b, c, e")[, c(2, 3)] ==
                    c(1.5 * 1.1, 5)))
    expect_true(all(dplyr::filter(ou, Genotype == "a, b, c, e")[, c(2, 3)] ==
                    c(1.3 * 1.5 * 1.1, 10 * 5)))
    oo <- evalAllGenotypesFitAndMut(fe, fm)
    expect_true(all(dplyr::filter(oo, Genotype == "a")[, c(2, 3)] ==
                    c(1, 10)))
    expect_true(all(dplyr::filter(oo, Genotype == "b")[, c(2, 3)] ==
                    c(1, 1)))
    expect_true(all(dplyr::filter(oo, Genotype == "c")[, c(2, 3)] ==
                    c(1, 5)))
    expect_true(all(dplyr::filter(oo, Genotype == "e")[, c(2, 3)] ==
                    c(1.1, 1)))
    expect_true(all(dplyr::filter(oo, Genotype == "a > b > c")[, c(2, 3)] ==
                    c(1.3 * 1.5, 10 * 5)))
    expect_true(all(dplyr::filter(oo, Genotype == "b > c > e")[, c(2, 3)] ==
                    c(1.5 * 1.1, 5)))
    expect_true(all(dplyr::filter(oo, Genotype == "e > b > c")[, c(2, 3)] ==
                    c(1.5 * 1.1, 5)))
    expect_true(all(dplyr::filter(oo, Genotype == "a > b > c > e")[, c(2, 3)] ==
                    c(1.3 * 1.5 * 1.1, 10 * 5)))
})

test_that("fails if genes in mutator not in fitness", {
    fe <- allFitnessEffects(noIntGenes = c("a" = 0.12,
                                           "c" = 0.14,
                                           "d" = 0.16,
                                           "e" = 0.11))
    fm4 <- allMutatorEffects(noIntGenes = c("a" = .010,
                                            "b" = .03,
                                            "d" = .08,
                                            "c" = .05))
    expect_error(oncoSimulIndiv(fe, muEF = fm4),
                 "Genes in mutatorEffects not present in fitnessEffects",
                 fixed = TRUE)
})

test_that("We cannot pass mutator/fitness objects to the wrong functions", {
    fe2 <- allFitnessEffects(noIntGenes =
                                 c(a1 = 0.1, a2 = 0.2,
                                   b1 = 0.01, b2 = 0.3, b3 = 0.2,
                                   c1 = 0.3, c2 = -0.2))
    fm2 <- allMutatorEffects(epistasis = c("A" = 5,
                                           "B" = 10,
                                           "C" = 3),
                             geneToModule = c("A" = "a1, a2",
                                              "B" = "b1, b2, b3",
                                              "C" = "c1, c2"))
    expect_error(evalAllGenotypesMut(fe2),
                 "You are trying to get the mutator effects of a fitness specification.",
                 fixed = TRUE)
    expect_error(evalAllGenotypes(fm2),
                 "You are trying to get the fitness of a mutator specification.",
                 fixed = TRUE)
    expect_error(evalGenotypeMut("a1, b2", fe2),
                 "You are trying to get the mutator effects of a fitness specification.",
                 fixed = TRUE)
    expect_error(evalGenotype("a2, c2", fm2),
                 "You are trying to get the fitness of a mutator specification.",
                 fixed = TRUE)
})

## FIXME: about 20 seconds. Move to long tests later
date()
test_that("Relative ordering of number of clones with mutator effects", {
    ## Can occasionally blow up with pE.f: pE not finite.
    pseed <-sample(9999999, 1)
    set.seed(pseed)
    cat("\n x1: the seed is", pseed, "\n")
    pops <- 20
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
    expect_true(median(summary(nc1)$NumClones) > median(summary(nc2)$NumClones))
    expect_true(median(summary(nc2)$NumClones) > median(summary(nc3)$NumClones))
    expect_true(mean(mutsPerClone(nc1)) > mean(mutsPerClone(nc2)))
    expect_true(mean(mutsPerClone(nc2)) > mean(mutsPerClone(nc3)))
    summary(nc1)[, c(1:3, 8:9)]
    summary(nc2)[, c(1:3, 8:9)]
    summary(nc3)[, c(1:3, 8:9)]
})
date()



date()
test_that("Relative ordering of number of clones with init mutant of mutator effects", {
    ## Here we stop on finalTime, not popSize
    ## Can occasionally blow up with pE.f: pE not finite.
    pseed <-sample(9999999, 1)
    set.seed(pseed)
    cat("\n x2bc: the seed is", pseed, "\n")
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
    expect_true( median(summary(nca)$NumClones) <
                 median(summary(ncb)$NumClones))
    expect_true(median(summary(ncb)$NumClones) <
                median(summary(ncc)$NumClones) )
    expect_true( median(summary(ncc)$NumClones) <
                 median(summary(ncd)$NumClones) )
    expect_true(mean(mutsPerClone(nca)) < mean(mutsPerClone(ncb)))
    expect_true(mean(mutsPerClone(ncb)) < mean(mutsPerClone(ncc)))
    expect_true(mean(mutsPerClone(ncc)) < mean(mutsPerClone(ncd)))
    summary(nca)[, c(1:3, 8:9)]
    summary(ncb)[, c(1:3, 8:9)]
    summary(ncc)[, c(1:3, 8:9)]
    summary(ncd)[, c(1:3, 8:9)]
})
date()


test_that("Relative ordering of number of clones with init mutant of mutator effects and s = 0", {
    ## Can occasionally blow up with pE.f: pE not finite.
    pseed <-sample(9999999, 1)
    set.seed(pseed)
    cat("\n x2cd: the seed is", pseed, "\n")
    pops <- 10
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
    expect_true( median(summary(nca)$NumClones) <
                 median(summary(ncb)$NumClones))
    expect_true(median(summary(ncb)$NumClones) <
                median(summary(ncc)$NumClones) )
    expect_true( median(summary(ncc)$NumClones) <
                 median(summary(ncd)$NumClones) )
    expect_true(mean(mutsPerClone(nca)) < mean(mutsPerClone(ncb)))
    expect_true(mean(mutsPerClone(ncb)) < mean(mutsPerClone(ncc)))
    expect_true(mean(mutsPerClone(ncc)) < mean(mutsPerClone(ncd)))
})

test_that("MCFL Relative ordering of number of clones with mutator effects", {
    ## Can occasionally blow up with pE.f: pE not finite.
    pseed <-sample(9999999, 1)
    set.seed(pseed)
    cat("\n mcx1: the seed is", pseed, "\n")
    pops <- 10
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
    expect_true(median(summary(nc1)$NumClones) > median(summary(nc2)$NumClones))
    expect_true(median(summary(nc2)$NumClones) > median(summary(nc3)$NumClones))
    summary(nc1)[, c(1:3, 8:9)]
    summary(nc2)[, c(1:3, 8:9)]
    summary(nc3)[, c(1:3, 8:9)]
    expect_true(mean(mutsPerClone(nc1)) > mean(mutsPerClone(nc2)))
    expect_true(mean(mutsPerClone(nc2)) > mean(mutsPerClone(nc3)))
})
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
    expect_true( median(summary(nca)$NumClones) <
                 median(summary(ncb)$NumClones))
    expect_true(median(summary(ncb)$NumClones) <
                median(summary(ncc)$NumClones) )
    expect_true( median(summary(ncc)$NumClones) <
                 median(summary(ncd)$NumClones) )
    summary(nca)[, c(1:3, 8:9)]
    summary(ncb)[, c(1:3, 8:9)]
    summary(ncc)[, c(1:3, 8:9)]
    summary(ncd)[, c(1:3, 8:9)]
    expect_true(mean(mutsPerClone(nca)) < mean(mutsPerClone(ncb)))
    expect_true(mean(mutsPerClone(ncb)) < mean(mutsPerClone(ncc)))
    expect_true(mean(mutsPerClone(ncc)) < mean(mutsPerClone(ncd)))
})



test_that("MCFL Relative ordering of number of clones with init mutant of mutator effects and s = 0", {
    ## Can occasionally blow up with pE.f: pE not finite.
    pseed <-sample(9999999, 1)
    set.seed(pseed)
    cat("\n mcx2cd: the seed is", pseed, "\n")
    pops <- 20
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
    expect_true( median(summary(nca)$NumClones) <
                 median(summary(ncb)$NumClones))
    expect_true(median(summary(ncb)$NumClones) <
                median(summary(ncc)$NumClones) )
    expect_true( median(summary(ncc)$NumClones) <
                 median(summary(ncd)$NumClones) )
    expect_true(mean(mutsPerClone(nca)) < mean(mutsPerClone(ncb)))
    expect_true(mean(mutsPerClone(ncb)) < mean(mutsPerClone(ncc)))
    expect_true(mean(mutsPerClone(ncc)) < mean(mutsPerClone(ncd)))
})



test_that("Relative ordering of number of clones with mut prop growth and init and scrambled names", {
    ## Can occasionally blow up with pE.f: pE not finite.
    pseed <- sample(99999999, 1)
    set.seed(pseed)
    cat("\n x2ef: the seed is", pseed, "\n")
    pops <- 10
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
    expect_true( median(summary(mpg)$NumClones) >
                 median(summary(mnpg)$NumClones))
    expect_true(median(summary(mpg)$NumClones) >
                median(summary(pg)$NumClones) )
    expect_true( median(summary(mnpg)$NumClones) >
                 median(summary(npg)$NumClones) )
    expect_true( median(summary(pg)$NumClones) >
                 median(summary(npg)$NumClones) )
    expect_true(mean(mutsPerClone(mpg)) > mean(mutsPerClone(mnpg)))
    expect_true(mean(mutsPerClone(mpg)) > mean(mutsPerClone(pg)))
    expect_true(mean(mutsPerClone(mnpg)) > mean(mutsPerClone(npg)))
    expect_true(mean(mutsPerClone(pg)) > mean(mutsPerClone(npg)))
})


test_that("McFL: Relative ordering of number of clones with mut prop growth and init and scrambled names", {
    ## Can occasionally blow up with pE.f: pE not finite.
    pseed <-sample(9999999, 1)
    set.seed(pseed)
    cat("\n x2gh: the seed is", pseed, "\n")
    pops <- 10
    ft <- 1
    lni <- 200
    no <- 1e3
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
                        initSize = no, model = "McFL",
                        initMutant = "thisistheagene",
                        onlyCancer = FALSE, mc.cores = 2)
    mnpg <- oncoSimulPop(pops, fe, muEF = fm1,
                         finalTime = ft,
                         mutationPropGrowth = FALSE,
                         initSize = no, model = "McFL",
                         initMutant = "thisistheagene",
                         onlyCancer = FALSE, mc.cores = 2)
    pg <- oncoSimulPop(pops, fe, 
                       finalTime = ft,
                       mutationPropGrowth = TRUE,
                       initSize = no, model = "McFL",
                       initMutant = "thisistheagene",
                       onlyCancer = FALSE, mc.cores = 2)
    npg <- oncoSimulPop(pops, fe, 
                        finalTime = ft,
                        mutationPropGrowth = FALSE,
                        initSize = no, model = "McFL",
                        initMutant = "thisistheagene",
                        onlyCancer = FALSE, mc.cores = 2)
      ## I once saw a weird thing
    expect_true(var(summary(mpg)$NumClones) > 1e-4)
    expect_true(var(summary(mnpg)$NumClones) > 1e-4)
    expect_true(var(summary(pg)$NumClones) > 1e-4)
    expect_true(var(summary(npg)$NumClones) > 1e-4)
    ## These are the real tests
    expect_true( median(summary(mpg)$NumClones) >
                 median(summary(mnpg)$NumClones))
    expect_true(median(summary(mpg)$NumClones) >
                median(summary(pg)$NumClones) )
    expect_true( median(summary(mnpg)$NumClones) >
                 median(summary(npg)$NumClones) )
    expect_true( median(summary(pg)$NumClones) >
                 median(summary(npg)$NumClones) )
    expect_true(mean(mutsPerClone(mpg)) > mean(mutsPerClone(mnpg)))
    expect_true(mean(mutsPerClone(mpg)) > mean(mutsPerClone(pg)))
    expect_true(mean(mutsPerClone(mnpg)) > mean(mutsPerClone(npg)))
    expect_true(mean(mutsPerClone(pg)) > mean(mutsPerClone(npg)))
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
test_that("Expect freq genotypes, mutator and var mut rates", {
    ## We test that mutator does not affect expected frequencies of
    ## mutated genes: they are given by the mutation rate of each gene.
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n u6: the seed is", pseed, "\n")
    pops <- 40
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
    ## of course, passing a mutator of 1 makes everything slow.
    mutator1 <- rep(1, lni + 3)
    ## pg1 <- rep(1e-5, lni + 3)
    pg1 <- runif(lni + 3, min = 1e-7, max = 1e-4) ## max should not be
                                                  ## huge here as mutator
                                                  ## is 34. Can get beyond
                                                  ## 1
    names(mutator1) <- sample(names(ni))
    names(pg1) <- sample(names(ni))
    mutator1["oreoisasabgene"] <- 34    ## 53
    m1 <- allMutatorEffects(noIntGenes = mutator1)
    ## have something with much larger mutation rate
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
    enom("oreoisasabgene", pg1, no, pops)
    snom("oreoisasabgene", m1.pg1.b)
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
    cat("\n u7: the seed is", pseed, "\n")
    pops <- 200
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

date()
test_that("Expect freq genotypes, mutator and var mut rates", {
    ## increase mutator, decrease max mu
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n u8: the seed is", pseed, "\n")
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
test_that("McFL, Expect freq genotypes, mutator and var mut rates", {
    
    ## We test that mutator does not affect expected frequencies of
    ## mutated genes: they are given by the mutation rate of each gene.
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n mcfu6: the seed is", pseed, "\n")
    pops <- 40
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
    m1.pg1.b <- oncoSimulPop(pops,
                           fe,
                           mu = pg1,
                           muEF = m1,
                           finalTime = ft,
                           mutationPropGrowth = FALSE,
                           initSize = no,
                           model = "McFL",
                           initMutant ="oreoisasabgene",
                           onlyCancer = FALSE, seed = NULL, mc.cores = 2)
    expect_true(sm("oreoisasabgene", m1.pg1.b) == totalind(m1.pg1.b))
    enom("oreoisasabgene", pg1, no, pops)
    snom("oreoisasabgene", m1.pg1.b)
    p.fail <- 1e-3
    expect_true(chisq.test(snom("oreoisasabgene", m1.pg1.b),
                           p = pnom("oreoisasabgene", pg1, no, pops))$p.value > p.fail)
    summary(m1.pg1.b)[, c(1:3, 8:9)]
    
})
date()



date()
test_that("McFL: Expect freq genotypes, mutator and var mut rates", {
    ## increase mutator
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n u8: the seed is", pseed, "\n")
    pops <- 200
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
    expect_true(medianNClones(pop10) < medianNClones(pop100))
    expect_true(mean(mutsPerClone(pop10)) < mean(mutsPerClone(pop100)))
    
})


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
    expect_true(medianNClones(pop10) < medianNClones(pop100))
    expect_true(mean(mutsPerClone(pop10)) < mean(mutsPerClone(pop100)))

})





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
                        model = "McFL",
                        finalTime = ft,
                        mutationPropGrowth = FALSE,
                        initSize = no,
                        initMutant = names(mutator10),
                        onlyCancer = FALSE, mc.cores = 2)
    pop100 <- oncoSimulPop(pops,
                        fe,
                        mu = muvector,
                        muEF = m100,
                        model = "McFL",                        
                        finalTime = ft,
                        mutationPropGrowth = FALSE,
                        initSize = no,
                        initMutant = names(mutator10),
                        onlyCancer = FALSE, mc.cores = 2)
    ## number of total mutations
    expect_true(smAnomPi(pop10, names(mutator10)) < smAnomPi(pop100, names(mutator100)))
    ## number of clones
    expect_true(medianNClones(pop10) < medianNClones(pop100))
    expect_true(mean(mutsPerClone(pop10)) < mean(mutsPerClone(pop100)))

})
date()


date()
test_that("McFL: Same mu vector, different mutator; diffs in number muts, larger t", {
    ## reproduction, death, and double and possibly triple mutants. We
    ## decrease init pop size to make this fast.
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n nm3: the seed is", pseed, "\n")
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
                        model = "McFL",                        
                        finalTime = ft,
                        mutationPropGrowth = FALSE,
                        initSize = no,
                        initMutant = names(mutator10),
                        onlyCancer = FALSE, mc.cores = 2)
    pop100 <- oncoSimulPop(pops,
                        fe,
                        mu = muvector,
                        muEF = m100,
                        model = "McFL",                        
                        finalTime = ft,
                        mutationPropGrowth = FALSE,
                        initSize = no,
                        initMutant = names(mutator10),
                        onlyCancer = FALSE, mc.cores = 2)
    ## number of total mutations
    expect_true(smAnomPi(pop10, names(mutator10)) < smAnomPi(pop100, names(mutator100)))
    ## number of clones
    expect_true(medianNClones(pop10) < medianNClones(pop100))
        expect_true(mean(mutsPerClone(pop10)) < mean(mutsPerClone(pop100)))
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
    expect_true(medianNClones(m1.pg1.a) > medianNClones(m1.pg1.b))
    expect_true(medianNClones(m1.pg1.b) > medianNClones(m1.pg1.c))
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
test_that(" MCFL Init with different mutators", {
    
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
                           model = "McFL",
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
                             model = "McFL",                             
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
                           model = "McFL",                           
                           initMutant = "nnhsisthecgene",
                           sampleEvery = 0.01, keepEvery = 5, seed = NULL,
                           onlyCancer = FALSE, mc.cores = 2)
    expect_true(medianNClones(m1.pg1.a) > medianNClones(m1.pg1.b))
    expect_true(medianNClones(m1.pg1.b) > medianNClones(m1.pg1.c))
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
    reps <- 40
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
    expect_true( median(summary(m1.mutator2)$NumClones) >
                 median(summary(m1.mutator1)$NumClones))
    expect_true( median(summary(m1.mutator1)$NumClones) >
                 median(summary(m1.mutator0)$NumClones))
    expect_true( median(summary(m2.mutator2)$NumClones) >
                 median(summary(m2.mutator1)$NumClones))
    expect_true( median(summary(m2.mutator1)$NumClones) >
                 median(summary(m2.mutator0)$NumClones))
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
    expect_true( median(summary(m1.mutator0)$NumClones) <
                 median(summary(m2.mutator0)$NumClones))
    expect_true( median(summary(m1.mutator1)$NumClones) <
                 median(summary(m2.mutator1)$NumClones))
    expect_true( median(summary(m1.mutator2)$NumClones) <
                 median(summary(m2.mutator2)$NumClones))
    expect_true( mean(mutsPerClone(m1.mutator0)) <
                 mean(mutsPerClone(m2.mutator0)))
    expect_true( mean(mutsPerClone(m1.mutator1)) <
                 mean(mutsPerClone(m2.mutator1)))
    expect_true( mean(mutsPerClone(m1.mutator2)) <
                 mean(mutsPerClone(m2.mutator2)))
})
date()



## FIXME: very slow, move to long later
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
    reps <- 40
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
    expect_true( median(summary(m1.mutator2)$NumClones) >
                 median(summary(m1.mutator1)$NumClones))
    expect_true( median(summary(m1.mutator1)$NumClones) >
                 median(summary(m1.mutator0)$NumClones))
    expect_true( median(summary(m2.mutator2)$NumClones) >
                 median(summary(m2.mutator1)$NumClones))
    expect_true( median(summary(m2.mutator1)$NumClones) >
                 median(summary(m2.mutator0)$NumClones))
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
    expect_true( median(summary(m1.mutator0)$NumClones) <
                 median(summary(m2.mutator0)$NumClones))
    expect_true( median(summary(m1.mutator1)$NumClones) <
                 median(summary(m2.mutator1)$NumClones))
    expect_true( median(summary(m1.mutator2)$NumClones) <
                 median(summary(m2.mutator2)$NumClones))
    expect_true( mean(mutsPerClone(m1.mutator0)) <
                 mean(mutsPerClone(m2.mutator0)))
    expect_true( mean(mutsPerClone(m1.mutator1)) <
                 mean(mutsPerClone(m2.mutator1)))
    expect_true( mean(mutsPerClone(m1.mutator2)) <
                 mean(mutsPerClone(m2.mutator2)))
})
date()



## FIXME move later to long
## Slow (~ 5 seconds) but tests modules of mutator nicely.
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
    expect_true( median(summary(b3)$NumClones) >
                 median(summary(b2)$NumClones))
    expect_true( median(summary(b2)$NumClones) >
                 median(summary(b1)$NumClones))
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
    expect_true( median(summary(b3)$NumClones) >
                 median(summary(b2)$NumClones))
    expect_true( median(summary(b2)$NumClones) >
                 median(summary(b1)$NumClones))
    expect_true( mean(mutsPerClone(b3)) >
                 mean(mutsPerClone(b2)))
    expect_true( mean(mutsPerClone(b2)) >
                 mean(mutsPerClone(b1)))
})
date()

    

date() 
test_that("Mutator, several modules differences", {
    pseed <- sample(99999999, 1)
    set.seed(pseed)
    cat("\n mmd1: the seed is", pseed, "\n")
    reps <- 10
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
                       seed = NULL, mc.cores = 2
                       )
    gc()
    ## summary(b2)[, c(1:3, 8:9)]
    ## summary(b1)[, c(1:3, 8:9)]
    ## mean(mutsPerClone(b2))
    ## mean(mutsPerClone(b1))
    ## This is, of course, affected by sampling only at end: we do not see
    ## the many intermediate events.
    expect_true( median(summary(b2)$NumClones) >
                 median(summary(b1)$NumClones))
    expect_true( mean(mutsPerClone(b2)) >
                 mean(mutsPerClone(b1)))
})
date()



date() 
test_that("Mutator, several modules differences, McFL", {
    pseed <- sample(99999999, 1)
    pseed <- 91339980
    set.seed(pseed)
    cat("\n mmd1: the seed is", pseed, "\n")
    reps <- 10
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
    expect_true( median(summary(b2)$NumClones) >
                 median(summary(b1)$NumClones))
    expect_true( mean(mutsPerClone(b2)) >
                 mean(mutsPerClone(b1)))
})
date()


date() 
test_that("Mutator, several modules differences, fitness eval", {
    ## the basis of what we do below, but fewer genes here
    ln <- 2 
    m1 <- 5
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
    evalAllGenotypesFitAndMut(f1, mut1, order = FALSE)
    evalAllGenotypesFitAndMut(f1, mut2, order = FALSE)
    ## FIXME: complete these evals
})

    



## 1.
### Use rT and order and epist and modules in fitness
##  and different epist and modules in mutator

## 2.  Like above, but add no interaction. Only some of fit. in mutator.


## 3. 4., Like above but have all genes be in both.

## 5. Like 2, but same modules in fitness and mutator.
## 6. Like 4, but same modules in fitness and mutator.

## 7. Check fail in 1., 2., 3., 4., where some genes in mutator not in
## fitness.

## 8. check fail if mutator and fitness not both given in the FitAndMut
## functions.


## docs:
##    - help
##  -fignete
##  - finish docs



######################################################################
######################################################################

## ## Some checks of C++ code with mutator and per-gene mutation rates

######################################################################
######################################################################



## set.seed(2) 
## ft <- 600
## lni <- 10
## no <- 5e3
## ni <- c(0, 0, 0, rep(0, lni))
## ## scramble around names
## names(ni) <- c("hereisoneagene",
##                "oreoisasabgene",
##                "nnhsisthecgene",
##                replicate(lni,
##                          paste(sample(letters, 12), collapse = "")))
## ni <- ni[order(names(ni))]
## ni <- sample(ni)
## fe <- allFitnessEffects(noIntGenes = ni)
## mutator1 <- rep(1, lni + 3)
## pg1 <- rep(1e-7, lni + 3)  ## what breaks is using 1e-7 instead of 1e-9
## ## scramble names of mutator and per-gene too
## names(mutator1) <- sample(names(ni))
## names(pg1) <- sample(names(ni))
## mutator1["hereisoneagene"] <- 100
## mutator1["oreoisasabgene"] <- 5
## mutator1["nnhsisthecgene"] <- 0.01
## m1 <- allMutatorEffects(noIntGenes = mutator1)
## pg1["hereisoneagene"] <- 1e-5
## pg1["oreoisasabgene"] <- 1e-7
## pg1["nnhsisthecgene"] <- 1e-10
## m1.pg1.a <- oncoSimulIndiv(fe,
##                            mu = pg1,
##                            muEF = m1,
##                            finalTime = ft,
##                            mutationPropGrowth = FALSE,
##                            initSize = no,
##                            initMutant = "hereisoneagene",
##                            onlyCancer = FALSE,
##                            verbosity = 6)
## m1.pg1.a

## ## Next is handy to compare
## mu.nop <- function(p, mu) {
##     nni <- names(ni)[p]
##     pmu <- which(names(mu) %in%  nni)
##     mu[-pmu]
## }

## ## genes a and b are 1 and 10
## ni[c(1, 10)]

## ##### iteration 928,
## ## parent:
## 100 * sum(mu.nop(c(10), pg1))
## ## child
## 100 * sum(mu.nop(c(10, 3), pg1))

## ### iteration 255
## ## parent
## 100 * 5 * sum(mu.nop(c(10, 1), pg1))
## ## child
## 100 * 5 * sum(mu.nop(c(10, 1, 11), pg1))





## ## ## Some checks on C++ of mutationPropGrowth and mutator
## ni <- rep(0.4, 20)
## names(ni) <- c("a", "b", "c", "d", paste0("n", 1:16))
## fe <- allFitnessEffects(noIntGenes = ni)
## fm6 <- allMutatorEffects(noIntGenes = c("a" = 15,
##                                         "b" = 15,
##                                         "c" = 15,
##                                         "d" = 15))
    
## set.seed(5)  ## if you use seed of 2, it blows up with huge birth
## ## rate and many mutations > 1
## ## But then, 1.4^20 is > 800!

## oncoSimulIndiv(fe, muEF = fm6, finalTime =30,
##                mutationPropGrowth = TRUE,
##                initSize = 1e4,
##                mu = 1e-06,
##                verbosity = 6,
##                onlyCancer = FALSE, seed = NULL)    
## ## ###### Iteration 40.
## ## ## mutation
## ## ## parent
## ## 20 * 1e-06
## ## ## child
## ## 1.4 * 15 * 1e-06 * 19
## ## ###### Iteration 39
## ## ## mutation
## ## ## parent
## ## 1.4 * 15 * 1e-06 * 19
## ## ## chlid
## ## 1.4 * 1.4 * 15 * 1e-06 * 18



## ## If you want to check the internal running of C++, use verbosity >=
## ## 3. For instance, at iteration 3 a new species is created from the
## ## wildtype. The new clone has mutated gene 1. Thus the new mutation rate
## ## should be

## 10 * 3 * 1e-6  ## 10 for the effect of a, 3 for the number of remaining
## ## genes, and 1e-6 for the baseline mutation rate

## That is 3e-5, as shown.

## At iteration 9 we create a mutant at 3 from wildtype and its mutation rate is
## 30 * 3 * 1e-6


## At iteration 254 we create species 1,2,4 from 1,2. The parent mutation is

## 10 * 20 * 2 * 1e-6 = 4e-4

## and the child's is

## 10 * 20 * 40 * 1 * 1e-6 = 0.008

## set.seed(1)
## fe <- allFitnessEffects(noIntGenes = c("a" = 0.12,
##                                        "b" = 0.14,
##                                        "c" = 0.16,
##                                        "d" = 0.11))
## fm6 <- allMutatorEffects(noIntGenes = c("a" = 10,
##                                         "b" = 20,
##                                         "c" = 30,
##                                         "d" = 40))
## oncoSimulIndiv(fe, muEF = fm6, finalTime =100,
##                mutationPropGrowth = FALSE,
##                initSize = 1e5,
##                verbosity = 6,
##                onlyCancer = FALSE,
##                seed = NULL)



cat(paste("\n Finished test.mutator.R test at", date()))













## ### These are too convoluted and trying to pick too tiny diffs. Substitute
## ###  by diffs initMutant and by per-gene-mut rates and mutator

## date()
## test_that("Num clones: Mutator and var mut rates and init and really scrambled names", {
##     pseed <- sample(9999999, 1)
##     set.seed(pseed)
##     cat("\n z2: the seed is", pseed, "\n")
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
##                            onlyCancer = FALSE, mc.cores = 2)
##     m1.pg2.a <- oncoSimulPop(pops,
##                            fe,
##                            mu = pg2,
##                            muEF = m1,
##                            finalTime = ft,
##                            mutationPropGrowth = FALSE,
##                            initSize = no,
##                            initMutant = "hereisoneagene",
##                            onlyCancer = FALSE, mc.cores = 2)
##     m1.pg1.c <- oncoSimulPop(pops,
##                            fe,
##                            mu = pg1,
##                            muEF = m1,
##                            finalTime = ft,
##                            mutationPropGrowth = FALSE,
##                            initSize = no,
##                            initMutant = "nnhsisthecgene",
##                            onlyCancer = FALSE, mc.cores = 2)
##     m1.pg2.c <- oncoSimulPop(pops,
##                            fe,
##                            mu = pg2,
##                            muEF = m1,
##                            finalTime = ft,
##                            mutationPropGrowth = FALSE,
##                            initSize = no,
##                            initMutant = "nnhsisthecgene",
##                            onlyCancer = FALSE, mc.cores = 2)
##     m1.pg1.b <- oncoSimulPop(pops,
##                            fe,
##                            mu = pg1,
##                            muEF = m1,
##                            finalTime = ft,
##                            mutationPropGrowth = FALSE,
##                            initSize = no,
##                            initMutant = "oreoisasabgene",
##                            onlyCancer = FALSE, mc.cores = 2)
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
##                            onlyCancer = FALSE, mc.cores = 2)
##     m2.pg2.a <- oncoSimulPop(pops,
##                            fe,
##                            mu = pg2,
##                            muEF = m2,
##                            finalTime = ft,
##                            mutationPropGrowth = FALSE,
##                            initSize = no,
##                            initMutant = "hereisoneagene",
##                            onlyCancer = FALSE, mc.cores = 2)
##     m2.pg1.c <- oncoSimulPop(pops,
##                            fe,
##                            mu = pg1,
##                            muEF = m2,
##                            finalTime = ft,
##                            mutationPropGrowth = FALSE,
##                            initSize = no,
##                            initMutant = "nnhsisthecgene",
##                            onlyCancer = FALSE, mc.cores = 2)
##     m2.pg2.c <- oncoSimulPop(pops,
##                            fe,
##                            mu = pg2,
##                            muEF = m2,
##                            finalTime = ft,
##                            mutationPropGrowth = FALSE,
##                            initSize = no,
##                            initMutant = "nnhsisthecgene",
##                            onlyCancer = FALSE, mc.cores = 2)
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
##     set.seed(pseed)
##     cat("\n x3: the seed is", pseed, "\n")
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
##                            onlyCancer = FALSE, mc.cores = 2)
##     m1.pg2.a <- oncoSimulPop(pops,
##                            fe,
##                            mu = pg2,
##                            muEF = m1,
##                            finalTime = ft, model = "McFL",
##                            mutationPropGrowth = FALSE,
##                            initSize = no,
##                            initMutant = "hereisoneagene",
##                            onlyCancer = FALSE, mc.cores = 2)
##     m1.pg1.c <- oncoSimulPop(pops,
##                            fe,
##                            mu = pg1,
##                            muEF = m1,
##                            finalTime = ft, model = "McFL",
##                            mutationPropGrowth = FALSE,
##                            initSize = no,
##                            initMutant = "nnhsisthecgene",
##                            onlyCancer = FALSE, mc.cores = 2)
##     m1.pg2.c <- oncoSimulPop(pops,
##                            fe,
##                            mu = pg2,
##                            muEF = m1,
##                            finalTime = ft, model = "McFL",
##                            mutationPropGrowth = FALSE,
##                            initSize = no,
##                            initMutant = "nnhsisthecgene",
##                            onlyCancer = FALSE, mc.cores = 2)
##     m1.pg1.b <- oncoSimulPop(pops,
##                            fe,
##                            mu = pg1,
##                            muEF = m1,
##                            finalTime = ft, model = "McFL",
##                            mutationPropGrowth = FALSE,
##                            initSize = no,
##                            initMutant = "oreoisasabgene",
##                            onlyCancer = FALSE, mc.cores = 2)
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
##                            onlyCancer = FALSE, mc.cores = 2)
##     m2.pg2.a <- oncoSimulPop(pops,
##                            fe,
##                            mu = pg2,
##                            muEF = m2,
##                            finalTime = ft, model = "McFL",
##                            mutationPropGrowth = FALSE,
##                            initSize = no,
##                            initMutant = "hereisoneagene",
##                            onlyCancer = FALSE, mc.cores = 2)
##     m2.pg1.c <- oncoSimulPop(pops,
##                            fe,
##                            mu = pg1,
##                            muEF = m2,
##                            finalTime = ft, model = "McFL",
##                            mutationPropGrowth = FALSE,
##                            initSize = no,
##                            initMutant = "nnhsisthecgene",
##                            onlyCancer = FALSE, mc.cores = 2)
##     m2.pg2.c <- oncoSimulPop(pops,
##                            fe,
##                            mu = pg2,
##                            muEF = m2,
##                            finalTime = ft, model = "McFL",
##                            mutationPropGrowth = FALSE,
##                            initSize = no,
##                            initMutant = "nnhsisthecgene",
##                            onlyCancer = FALSE, mc.cores = 2)
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
##     set.seed(pseed)
##     cat("\n x4: the seed is", pseed, "\n")
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
##                            onlyCancer = FALSE, seed = NULL, mc.cores = 2)
##     m1.pg2.a <- oncoSimulPop(pops,
##                            fe,
##                            mu = pg2,
##                            muEF = m1,
##                            finalTime = ft,
##                            mutationPropGrowth = TRUE,
##                            initSize = no,
##                            initMutant = "hereisoneagene",
##                            onlyCancer = FALSE, seed = NULL, mc.cores = 2)
##     m1.pg1.c <- oncoSimulPop(pops,
##                            fe,
##                            mu = pg1,
##                            muEF = m1,
##                            finalTime = ft,
##                            mutationPropGrowth = TRUE,
##                            initSize = no,
##                            initMutant = "nnhsisthecgene",
##                            onlyCancer = FALSE, seed = NULL, mc.cores = 2)
##     m1.pg2.c <- oncoSimulPop(pops,
##                            fe,
##                            mu = pg2,
##                            muEF = m1,
##                            finalTime = ft,
##                            mutationPropGrowth = TRUE,
##                            initSize = no,
##                            initMutant = "nnhsisthecgene",
##                            onlyCancer = FALSE, seed = NULL, mc.cores = 2)
##     m1.pg1.b <- oncoSimulPop(pops,
##                            fe,
##                            mu = pg1,
##                            muEF = m1,
##                            finalTime = ft,
##                            mutationPropGrowth = TRUE,
##                            initSize = no,
##                            initMutant = "oreoisasabgene",
##                            onlyCancer = FALSE, seed = NULL, mc.cores = 2)
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
##                            onlyCancer = FALSE, seed = NULL, mc.cores = 2)
##     m2.pg2.a <- oncoSimulPop(pops,
##                            fe,
##                            mu = pg2,
##                            muEF = m2,
##                            finalTime = ft,
##                            mutationPropGrowth = TRUE,
##                            initSize = no,
##                            initMutant = "hereisoneagene",
##                            onlyCancer = FALSE, seed = NULL, mc.cores = 2)
##     m2.pg1.c <- oncoSimulPop(pops,
##                            fe,
##                            mu = pg1,
##                            muEF = m2,
##                            finalTime = ft,
##                            mutationPropGrowth = TRUE,
##                            initSize = no,
##                            initMutant = "nnhsisthecgene",
##                            onlyCancer = FALSE, seed = NULL, mc.cores = 2)
##     m2.pg2.c <- oncoSimulPop(pops,
##                            fe,
##                            mu = pg2,
##                            muEF = m2,
##                            finalTime = ft,
##                            mutationPropGrowth = TRUE,
##                            initSize = no,
##                            initMutant = "nnhsisthecgene",
##                            onlyCancer = FALSE, seed = NULL, mc.cores = 2)
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
##     set.seed(pseed)
##     cat("\n x5: the seed is", pseed, "\n")
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
##                            onlyCancer = FALSE, seed = NULL, mc.cores = 2)
##     m1.pg2.a <- oncoSimulPop(pops,
##                            fe,
##                            mu = pg2,
##                            muEF = m1,
##                            finalTime = ft, model = "McFL",
##                            mutationPropGrowth = TRUE,
##                            initSize = no,
##                            initMutant = "hereisoneagene",
##                            onlyCancer = FALSE, seed = NULL, mc.cores = 2)
##     m1.pg1.c <- oncoSimulPop(pops,
##                            fe,
##                            mu = pg1,
##                            muEF = m1,
##                            finalTime = ft, model = "McFL",
##                            mutationPropGrowth = TRUE,
##                            initSize = no,
##                            initMutant = "nnhsisthecgene",
##                            onlyCancer = FALSE, seed = NULL, mc.cores = 2)
##     m1.pg2.c <- oncoSimulPop(pops,
##                            fe,
##                            mu = pg2,
##                            muEF = m1,
##                            finalTime = ft, model = "McFL",
##                            mutationPropGrowth = TRUE,
##                            initSize = no,
##                            initMutant = "nnhsisthecgene",
##                            onlyCancer = FALSE, seed = NULL, mc.cores = 2)
##     m1.pg1.b <- oncoSimulPop(pops,
##                            fe,
##                            mu = pg1,
##                            muEF = m1,
##                            finalTime = ft, model = "McFL",
##                            mutationPropGrowth = TRUE,
##                            initSize = no,
##                            initMutant = "oreoisasabgene",
##                            onlyCancer = FALSE, seed = NULL, mc.cores = 2)
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
##                            onlyCancer = FALSE, seed = NULL, mc.cores = 2)
##     m2.pg2.a <- oncoSimulPop(pops,
##                            fe,
##                            mu = pg2,
##                            muEF = m2,
##                            finalTime = ft, model = "McFL",
##                            mutationPropGrowth = TRUE,
##                            initSize = no,
##                            initMutant = "hereisoneagene",
##                            onlyCancer = FALSE, seed = NULL, mc.cores = 2)
##     m2.pg1.c <- oncoSimulPop(pops,
##                            fe,
##                            mu = pg1,
##                            muEF = m2,
##                            finalTime = ft, model = "McFL",
##                            mutationPropGrowth = TRUE,
##                            initSize = no,
##                            initMutant = "nnhsisthecgene",
##                            onlyCancer = FALSE, seed = NULL, mc.cores = 2)
##     m2.pg2.c <- oncoSimulPop(pops,
##                            fe,
##                            mu = pg2,
##                            muEF = m2,
##                            finalTime = ft, model = "McFL",
##                            mutationPropGrowth = TRUE,
##                            initSize = no,
##                            initMutant = "nnhsisthecgene",
##                            onlyCancer = FALSE, seed = NULL, mc.cores = 2)
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
##     set.seed(pseed)
##     cat("\n u5: the seed is", pseed, "\n")
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
##                            onlyCancer = FALSE, seed = NULL, mc.cores = 2)
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
##                            onlyCancer = FALSE, seed = NULL, mc.cores = 2)
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
##     set.seed(pseed)
##     cat("\n l-mmd2: the seed is", pseed, "\n")
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
##                        onlyCancer = FALSE,
##                        initSize = no,
##                        finalTime = ft,
##                        seed = NULL, keepEvery = ft
##                        )
##     gc()
##     b2 <- oncoSimulPop(reps,
##                        f1,
##                        mu = mu,
##                        muEF = mut2,
##                        onlyCancer = FALSE,
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

