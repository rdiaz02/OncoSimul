### Since some tests are slow and some tests are very fragile, for now I
### leave date() and seed()

### This test takes about 10 seconds
cat(paste("\n Starting test.mutator.R test at", date()))
date()

## RNGkind("L'Ecuyer-CMRG") ## for the mclapplies

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

## The threshold is rather arbitrary. Of course, you would not expect most
## of them to pass if they were false. But note that we cannot set too
## tiny a threshold unless we increase number of repetitions a lot, and
## that makes the test run for very long.

p.value.threshold <- 0.01

date()
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
    bb <- oncoSimulIndiv(fe1, mu = muvar2, onlyCancer = FALSE, detectionProb = NA,
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
                                     onlyCancer = FALSE, detectionProb = NA,
                                     initSize = 1000,
                                     finalTime = 1,
                                     seed = NULL, mc.cores = 2)),
                  "Population of OncoSimul",
                  fixed = TRUE)
    expect_output(print(oncoSimulPop(4,
                                     fe,
                                     mu = moo,
                                     onlyCancer = FALSE, detectionProb = NA,
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
    expect_equal(evalGenotypeFitAndMut("a, b, e", fe, fm),
                     c(1.3 * 1.1, 10))
    expect_equal(evalGenotypeFitAndMut("a, b, c, e", fe, fm),
                     c(1.3 * 1.5 * 1.1, 10 * 5))
})

test_that("expect output oncoSimulIndiv", {
    fe <- allFitnessEffects(noIntGenes = c("a" = 0.2,
                                           "c" = 0.4,
                                           "d" = 0.6,
                                           "e" = 0.1))
    fm <- allMutatorEffects(noIntGenes = c("a" = 10,
                                           "c" = 5))
    expect_output(print(oncoSimulIndiv(fe, muEF = fm, sampleEvery = 0.01,
                                 keepEvery = 5)),
                  "Individual OncoSimul trajectory",
                  fixed = TRUE)
    expect_output(print(oncoSimulIndiv(fe, muEF = fm, sampleEvery = 0.01,
                                 keepEvery = 5)),
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



test_that("eval mut genotypes, echo", {
    fe <- allFitnessEffects(epistasis = c("a : b" = 0.3,
                                          "b : c" = 0.5),
                            noIntGenes = c("e" = 0.1),
                            drvNames = c(letters[1:3]))
    fm <- allMutatorEffects(noIntGenes = c("a" = 10,
                                           "c" = 5))
    expect_output(evalGenotypeMut("a, c", fm, echo = TRUE),
                  "Mutation rate product", fixed = TRUE)
})


test_that("evaluating genotype and mutator, Bozic", {
    fe <- allFitnessEffects(epistasis = c("a : b" = 0.3,
                                          "b : c" = 0.5),
                            drvNames = letters[1:3])
    fm <- allMutatorEffects(noIntGenes = c("a" = 10,
                                           "c" = 5))
    ou <- evalAllGenotypesFitAndMut(fe, fm, order = FALSE,
                                    model = "Bozic")
    expect_equivalent(ou[, 2],
                      c(1, 1, 1, 1 - .3, 1, 1 -.5, (1 -.3) * (1 - .5))
                      )
    expect_equivalent(ou[, 3],
                      c(10, 1, 5, 10, 10 * 5, 5, 10 * 5)
                      )
})


test_that("mut and fitness both needed when needed", {
    fe <- allFitnessEffects(epistasis = c("a : b" = 0.3,
                                          "b : c" = 0.5),
                            noIntGenes = c("e" = 0.1),
                            drvNames = c(letters[1:3]))
    fm <- allMutatorEffects(noIntGenes = c("a" = 10,
                                           "c" = 5))
    expect_error(evalGenotypeFitAndMut("a, b", fe),
                 'argument "mutatorEffects" is missing',
                 fixed = TRUE)
    expect_error(evalGenotypeFitAndMut("a, b", fm),
                 "genotype contains NAs or genes not in fitnessEffects",
                 fixed = TRUE)
    expect_error(evalGenotypeFitAndMut("a, b", mutatorEffects = fm),
                 'argument "fitnessEffects" is missing',
                 fixed = TRUE)
    expect_error(evalAllGenotypesFitAndMut(fe, order = TRUE),
                 'argument "mutatorEffects" is missing',
                 fixed = TRUE)
    expect_error(evalAllGenotypesFitAndMut(fe, order = FALSE),
                 'argument "mutatorEffects" is missing',
                 fixed = TRUE)
    expect_error(evalAllGenotypesFitAndMut(fm),
                 'argument "mutatorEffects" is missing',
                 fixed = TRUE)
    expect_error(evalAllGenotypesFitAndMut(mutatorEffects = fm),
                'argument "fitnessEffects" is missing',
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




test_that("evaluating genotype and mutator", {
    fe <- allFitnessEffects(epistasis = c("a : b" = 0.3,
                                          "b : c" = 0.5),
                            noIntGenes = c("e" = 0.1),
                            drvNames = letters[1:3])
    fm <- allMutatorEffects(noIntGenes = c("a" = 10,
                                           "c" = 5))
    ou <- evalAllGenotypesFitAndMut(fe, fm, order = FALSE)
    expect_equivalent(dplyr::filter(ou, Genotype == "a")[, c(2, 3)],
                    c(1, 10))
    expect_equivalent(dplyr::filter(ou, Genotype == "b")[, c(2, 3)],
                    c(1, 1))
    expect_equivalent(dplyr::filter(ou, Genotype == "c")[, c(2, 3)],
                    c(1, 5))
    expect_equivalent(dplyr::filter(ou, Genotype == "e")[, c(2, 3)],
                    c(1.1, 1))
    expect_equivalent(dplyr::filter(ou, Genotype == "a, b, c")[, c(2, 3)],
                    c(1.3 * 1.5, 10 * 5))
    expect_equivalent(dplyr::filter(ou, Genotype == "b, c, e")[, c(2, 3)],
                    c(1.5 * 1.1, 5))
    expect_equivalent(dplyr::filter(ou, Genotype == "a, b, c, e")[, c(2, 3)],
                    c(1.3 * 1.5 * 1.1, 10 * 5))
    oo <- evalAllGenotypesFitAndMut(fe, fm, order = TRUE)
    expect_equivalent(dplyr::filter(oo, Genotype == "a")[, c(2, 3)],
                    c(1, 10))
    expect_equivalent(dplyr::filter(oo, Genotype == "b")[, c(2, 3)],
                    c(1, 1))
    expect_equivalent(dplyr::filter(oo, Genotype == "c")[, c(2, 3)],
                    c(1, 5))
    expect_equivalent(dplyr::filter(oo, Genotype == "e")[, c(2, 3)],
                    c(1.1, 1))
    expect_equivalent(dplyr::filter(oo, Genotype == "a > b > c")[, c(2, 3)],
                    c(1.3 * 1.5, 10 * 5))
    expect_equivalent(dplyr::filter(oo, Genotype == "b > c > e")[, c(2, 3)],
                    c(1.5 * 1.1, 5))
    expect_equivalent(dplyr::filter(oo, Genotype == "e > b > c")[, c(2, 3)],
                    c(1.5 * 1.1, 5))
    expect_equivalent(dplyr::filter(oo, Genotype == "a > b > c > e")[, c(2, 3)],
                    c(1.3 * 1.5 * 1.1, 10 * 5))
})


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
    e1 <- evalAllGenotypesFitAndMut(f1, mut1, order = FALSE,
                                    addwt = TRUE)
    e2 <- evalAllGenotypesFitAndMut(f1, mut2, order = FALSE,
                                    addwt = TRUE)
    expect_identical(e1$Fitness, rep(1, 64))
    expect_identical(e1$MutatorFactor, c(1, rep(5, 63)))
    expect_identical(e2$Fitness, rep(1, 64))
    expect_identical(
        e2$MutatorFactor,
                    c(1, rep(5, 7),
                      rep(25, 8), 5,
                      rep(25, 4), 5,
                      rep(25, 5),
                      rep(125, 4), 25, 25,
                      rep(125, 4), rep(25, 6),
                      rep(125, 4), 25,
                      rep(125, 8), 25,
                      rep(125, 7))
    )
})
date()





date()
test_that("McFL: Relative ordering of number of clones with mut prop growth and init and scrambled names", {
    max.tries <- 4  
    for(tries in 1:max.tries) {
    ## Can occasionally blow up with pE.f: pE not finite.
    cat("\n x2gh: a runif is", runif(1), "\n")
    pops <- 50
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
    fm1 <- allMutatorEffects(noIntGenes = c("thisistheagene" = 8))
    mpg <- oncoSimulPop(pops, fe, muEF = fm1,
                        finalTime = ft,
                        mutationPropGrowth = TRUE,
                        initSize = no, model = "McFL",
                        initMutant = "thisistheagene",
                        onlyCancer = FALSE, detectionProb = NA, seed = NULL, mc.cores = 2)
    mnpg <- oncoSimulPop(pops, fe, muEF = fm1,
                         finalTime = ft,
                         mutationPropGrowth = FALSE,
                         initSize = no, model = "McFL",
                         initMutant = "thisistheagene",
                         onlyCancer = FALSE, detectionProb = NA, seed = NULL, mc.cores = 2)
    pg <- oncoSimulPop(pops, fe, 
                       finalTime = ft,
                       mutationPropGrowth = TRUE,
                       initSize = no, model = "McFL",
                       initMutant = "thisistheagene",
                       onlyCancer = FALSE, detectionProb = NA, seed = NULL, mc.cores = 2)
    npg <- oncoSimulPop(pops, fe, 
                        finalTime = ft,
                        mutationPropGrowth = FALSE,
                        initSize = no, model = "McFL",
                        initMutant = "thisistheagene",
                        onlyCancer = FALSE, detectionProb = NA, seed = NULL, mc.cores = 2)
    ##   ## I once saw a weird thing
    ## expect_true(var(summary(mpg)$NumClones) > 1e-4)
    ## expect_true(var(summary(mnpg)$NumClones) > 1e-4)
    ## expect_true(var(summary(pg)$NumClones) > 1e-4)
    ## expect_true(var(summary(npg)$NumClones) > 1e-4)
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
        if( T1 && T3 && T4 && T5 && T6 && T7 && T8 ) break;    
    }
    cat(paste("\n done tries", tries, "\n"))
    expect_true(T1 && T3 && T4 && T5 && T6 && T7 && T8)
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
    max.tries <- 4
    for(tries in 1:max.tries) {

    
    
    cat("\n u6: a runif is", runif(1), "\n")
    pops <- 20
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
                           onlyCancer = FALSE, detectionProb = NA, seed = NULL, mc.cores = 2)
    expect_true(sm("oreoisasabgene", m1.pg1.b) == totalind(m1.pg1.b))
    enom("oreoisasabgene", pg1, no, pops)
    snom("oreoisasabgene", m1.pg1.b)
    p.fail <- 1e-3
    T1 <- (chisq.test(snom("oreoisasabgene", m1.pg1.b),
                      p = pnom("oreoisasabgene", pg1, no, pops))$p.value > p.fail)
        if(T1) break;
    }
    cat(paste("\n done tries", tries, "\n"))
    expect_true(T1)
})
date()


date()
test_that("McFL, Expect freq genotypes, mutator and var mut rates", {
    max.tries <- 4
    for(tries in 1:max.tries) {
    
    ## We test that mutator does not affect expected frequencies of
    ## mutated genes: they are given by the mutation rate of each gene.
    
    
    cat("\n mcfu6: a runif is", runif(1), "\n")
    pops <- 20
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
                           onlyCancer = FALSE, detectionProb = NA, seed = NULL, mc.cores = 2)
    expect_true(sm("oreoisasabgene", m1.pg1.b) == totalind(m1.pg1.b))
    enom("oreoisasabgene", pg1, no, pops)
    snom("oreoisasabgene", m1.pg1.b)
    p.fail <- 1e-3
    T1 <- (chisq.test(snom("oreoisasabgene", m1.pg1.b),
                           p = pnom("oreoisasabgene", pg1, no, pops))$p.value > p.fail)
        summary(m1.pg1.b)[, c(1:3, 8:9)]
        if(T1 ) break;
    }

    cat(paste("\n done tries", tries, "\n"))
    expect_true( T1 )
})
date()






date()
test_that("McFL: Same mu vector, different mutator; diffs in number muts, tiny t", {
    max.tries <- 4
    for(tries in 1:max.tries) {


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
    pop10 <- oncoSimulPop(pops,
                        fe,
                        mu = muvector,
                        muEF = m10,
                        model = "McFL",
                        finalTime = ft,
                        mutationPropGrowth = FALSE,
                        initSize = no,
                        initMutant = names(mutator10),
                        onlyCancer = FALSE, detectionProb = NA, seed = NULL, mc.cores = 2)
    pop100 <- oncoSimulPop(pops,
                        fe,
                        mu = muvector,
                        muEF = m100,
                        model = "McFL",                        
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
        if( T1 && T2 && T3 ) break;
    }
    cat(paste("\n done tries", tries, "\n"))
    expect_true(T1 && T2 && T3 )
})
date()


date()
test_that("McFL: Same mu vector, different mutator; diffs in number muts, larger t", {
    ## reproduction, death, and double and possibly triple mutants. We
    ## decrease init pop size to make this fast.
        max.tries <- 4
    for(tries in 1:max.tries) {


    
    
    cat("\n nm3: a runif is", runif(1), "\n")
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
    pop10 <- oncoSimulPop(pops,
                        fe,
                        mu = muvector,
                        muEF = m10,
                        model = "McFL",                        
                        finalTime = ft,
                        mutationPropGrowth = FALSE,
                        initSize = no,
                        initMutant = names(mutator10),
                        onlyCancer = FALSE, detectionProb = NA, seed = NULL, mc.cores = 2)
    pop100 <- oncoSimulPop(pops,
                        fe,
                        mu = muvector,
                        muEF = m100,
                        model = "McFL",                        
                        finalTime = ft,
                        mutationPropGrowth = FALSE,
                        initSize = no,
                        initMutant = names(mutator10),
                        onlyCancer = FALSE, detectionProb = NA, seed = NULL, mc.cores = 2)
    ## number of total mutations
    expect_true(smAnomPi(pop10, names(mutator10)) < smAnomPi(pop100, names(mutator100)))
    ## number of clones
    T1 <- (wilcox.test(NClones(pop100),  NClones(pop10), alternative = "greater")$p.value < p.value.threshold)
    T2 <- (t.test(mutsPerClone(pop100), mutsPerClone(pop10), alternative = "greater")$p.value < p.value.threshold)
        if(T1 && T2) break;
    }
    cat(paste("\n done tries", tries, "\n"))
    expect_true(T1 && T2 )
    })
date()




date()
test_that(" MCFL Init with different mutators", {
        max.tries <- 4
    for(tries in 1:max.tries) {


    
    
    cat("\n mcfl_z2: a runif is", runif(1), "\n")
    pops <- 20
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
                           onlyCancer = FALSE, detectionProb = NA, mc.cores = 2)
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
                             onlyCancer = FALSE, detectionProb = NA, mc.cores = 2)
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
                           onlyCancer = FALSE, detectionProb = NA, mc.cores = 2)
    T1  <- (wilcox.test(NClones(m1.pg1.a),  NClones(m1.pg1.b), alternative = "greater")$p.value < p.value.threshold)
    T2  <- (wilcox.test(NClones(m1.pg1.b),  NClones(m1.pg1.c), alternative = "greater")$p.value < p.value.threshold)
    T3  <- (t.test(mutsPerClone(m1.pg1.a), mutsPerClone(m1.pg1.b), alternative = "greater")$p.value < p.value.threshold)
    T4  <- (t.test(mutsPerClone(m1.pg1.b), mutsPerClone(m1.pg1.c), alternative = "greater")$p.value < p.value.threshold)
    T5  <- (smAnomPi(m1.pg1.a, "hereisoneagene") >
                smAnomPi(m1.pg1.b, "oreoisasabgene"))
    T6  <- (smAnomPi(m1.pg1.b, "oreoisasabgene") >
                smAnomPi(m1.pg1.c, "nnhsisthecgene"))
    summary(m1.pg1.a)[, c(1:3, 8:9)]
    summary(m1.pg1.b)[, c(1:3, 8:9)]
    summary(m1.pg1.c)[, c(1:3, 8:9)]
        if(T1 && T2 && T3 && T4 && T5 && T6) break;
    }
            cat(paste("\n done tries", tries, "\n"))
    expect_true(T1 && T2 && T3 && T4 && T5 && T6)

})
date()



 

date() 
test_that("Mutator, several modules differences", {
    max.tries <- 4 
    for(tries in 1:max.tries) {
    cat("\n mmdSM1: a runif is", runif(1), "\n")
    reps <- 60
    no <- 5e3
    ft <- 50 ## you need it large enough to get enough hits
    mu <- 1e-5
    ln <- 100
    m1 <- 5 ## if this is too large, easy to get it to blow.
    ## Having three renders it too unpredictable in time.
    ## ni <- rep(0, 3 * ln)
    ## gna <- paste0("a", 1:ln)
    ## gnb <- paste0("b", 1:ln)
    ## gnc <- paste0("c", 1:ln)
    ## names(ni) <- c(gna, gnb, gnc)
    ## gn1 <- paste(c(gna, gnb, gnc), collapse = ", ")
    ## gna <- paste(gna, collapse = ", ")
    ## gnb <- paste(gnb, collapse = ", ")
    ## gnc <- paste(gnc, collapse = ", ")
    ## mut1 <- allMutatorEffects(epistasis = c("A" = m1),
    ##                           geneToModule = c("A" = gn1))
    ## mut2 <- allMutatorEffects(epistasis = c("A" = m1,
    ##                                         "B" = m1,
    ##                                         "C" = m1),
    ##                           geneToModule = c("A" = gna,
    ##                                            "B" = gnb,
    ##                                            "C" = gnc))
    ni <- rep(0, 2 * ln)
    gna <- paste0("a", 1:ln)
    gnb <- paste0("b", 1:ln)
    names(ni) <- c(gna, gnb  )
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
    b1 <- oncoSimulPop(reps,
                       f1,
                       mu = mu,
                       muEF = mut1,
                       onlyCancer = FALSE, detectionProb = NA,
                       initSize = no,
                       finalTime = ft,
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
                       seed = NULL, mc.cores = 2
                       )
    gc()
    ## p.value.threshold <- 0.05 ## diffs here small, unless large reps
    ## summary(b2)[, c(1:3, 8:9)]
    ## summary(b1)[, c(1:3, 8:9)]
    ## mean(mutsPerClone(b2))
    ## mean(mutsPerClone(b1))
    ## This is, of course, affected by sampling only at end: we do not see
    ## the many intermediate events.
    T1 <- ( wilcox.test(summary(b2)$NumClones,
                             summary(b1)$NumClones, alternative = "greater")$p.value < p.value.threshold)
    T2 <- (t.test(mutsPerClone(b2), mutsPerClone(b1), alternative = "greater")$p.value < p.value.threshold)
    ## it very rarely fails; what are the p-values?
    print(suppressWarnings(wilcox.test(summary(b2)$NumClones,
                                       summary(b1)$NumClones, alternative = "greater")$p.value))
    print(suppressWarnings(t.test(mutsPerClone(b2), mutsPerClone(b1), alternative = "greater")$p.value))
    if( T1 && T2 ) break;
    }
    cat(paste ("\n done tries ", tries, "\n"))
    expect_true( (T1 && T2) )
    })
date()



date() 
test_that("Mutator and mutPropGrowth, mcfl", {
        max.tries <- 4
    for(tries in 1:max.tries) {


    ## we stop on size
    ## Note that names of modules are different. Just for fun.
    cat("\n mmpg_mcfl: a runif is", runif(1), "\n")
    reps <- 14
    no <- 5e3
    ds <- 1.5e4
    ft <- 1000 
    mu <- 2e-5
    ln <- 50
    m1 <- 1 ## if this is too large, easy to get it to blow.
    m50 <- 50
    gna <- paste0("a", 1:ln)
    gna <- paste(gna, collapse = ", ")
    f1 <- allFitnessEffects(epistasis = c("A" = 1),
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
                          sampleEvery = 0.05,
                          keepEvery = 5,
                          seed = NULL, mc.cores = 2, model = "McFL"
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
                          sampleEvery = 0.05,
                          keepEvery = 5,
                          seed = NULL, mc.cores = 2, model = "McFL"
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
                          sampleEvery = 0.05,
                          keepEvery = 5,
                          seed = NULL, mc.cores = 2, model = "McFL"
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
                          sampleEvery = 0.05,
                          keepEvery = 5,
                          seed = NULL, mc.cores = 2, model = "McFL"
                          )
    ##gc()
    summary(m1.pg)[, c(1:3, 8:9)]
    summary(m50.pg)[, c(1:3, 8:9)]
    summary(m1.npg)[, c(1:3, 8:9)]
    summary(m50.npg)[, c(1:3, 8:9)]
    ## Over mutator, as we have mutPropGrowth, clones, etc, increase
    T1  <- ( wilcox.test(summary(m1.pg)$NumClones,
                             summary(m1.npg)$NumClones,
                             alternative = "greater")$p.value < p.value.threshold)
    T2  <- (t.test(mutsPerClone(m1.pg),
                       mutsPerClone(m1.npg),
                       alternative = "greater")$p.value < p.value.threshold)
    T3  <- ( wilcox.test(summary(m50.pg)$NumClones,
                             summary(m50.npg)$NumClones,
                             alternative = "greater")$p.value < p.value.threshold)
    T4  <- (t.test(mutsPerClone(m50.pg),
                       mutsPerClone(m50.npg),
                       alternative = "greater")$p.value < p.value.threshold)
    ## Over mutPropGrowth, as we increase mutator, clones, etc, increase
    T5  <- ( wilcox.test(summary(m50.pg)$NumClones,
                             summary(m1.pg)$NumClones,
                             alternative = "greater")$p.value < p.value.threshold)
    T6  <- (t.test(mutsPerClone(m50.pg),
                       mutsPerClone(m1.pg),
                       alternative = "greater")$p.value < p.value.threshold)
    T7  <- ( wilcox.test(summary(m50.npg)$NumClones,
                             summary(m1.npg)$NumClones,
                             alternative = "greater")$p.value < p.value.threshold)
    T8  <- (t.test(mutsPerClone(m50.npg),
                       mutsPerClone(m1.npg),
                       alternative = "greater")$p.value < p.value.threshold)
        if( T1 && T2 && T3 && T4 && T5 && T6 && T7 && T8) break;
    }
            cat(paste("\n done tries", tries, "\n"))
    expect_true(T1 && T2 && T3 && T4 && T5 && T6 && T7 && T8)
    })
date()









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



cat(paste("\n Finished test.mutator.R test at", date(), "\n"))



