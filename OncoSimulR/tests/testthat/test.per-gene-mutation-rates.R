cat("\n Starting at", date(), "\n") ## whole file takes about 30 seconds
## When submitting, probably move half of the tests (mcfl?) to the "long" file.

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



mutsPerCloneLast <- function(x, per.pop.mean = TRUE) {
    ## Only clones alive in the last period
    aliveLast <- function(u) {
        pbt <- u$pops.by.time
        which(pbt[nrow(pbt), -1] >= 1)
    }
    perCl <- function(z) {
        this <- aliveLast(z)
        unlist(lapply(z$GenotypesWDistinctOrderEff[this],  length))
    }
    perCl2 <- function(z) {
        this <- aliveLast(z)
        mean(unlist(lapply(z$GenotypesWDistinctOrderEff[this], length)))
    }
    if(per.pop.mean)    
        unlist(lapply(x, function(u) perCl2(u)))
    else
        lapply(x, function(u) perCl(u))
}

## we could have used this below . Oh well
## totalind <- function(out) {
##     ## total num indivs
##   sum(unlist(lapply(out, function(x) x$TotalPopSize)))  
## }


RNGkind("L'Ecuyer-CMRG") ## for the mclapplies
## RNGkind("Mersenne-Twister")

test_that("single named gene in mut. fail", {
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n s01: the seed is", pseed, "\n")
    fea9 <- allFitnessEffects(noIntGenes = c("m" = 0.1, "D" = 0.3))
    muvar <- c("m" = 1e-5)
    expect_error(oncoSimulIndiv(fea9, mu = muvar),
                  "A length 1 mutation, but named",
                  fixed = TRUE)
} )

test_that("Per-gene mutation rates with old poset format, fail", {
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n s02: the seed is", pseed, "\n")
    data(examplePosets)
    p701 <- examplePosets[["p701"]]
    muvar <- c(rep(1e-5, 4), rep(1e-6, 3))
    names(muvar) <- letters[1:7]
    expect_error(oncoSimulIndiv(p701, mu = muvar),
                 "Per-gene mutation rates cannot be used with the old poset format")
} )

test_that("Only no-int, and sorting", {
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n s03: the seed is", pseed, "\n")
    fea9 <- allFitnessEffects(noIntGenes = c("m" = 0.1, "D" = 0.3))
    ## OncoSimulR:::allNamedGenes(fea9)
    muvar <- c("m" = 1e-5, "D" = 1e-7)
    expect_output(oncoSimulIndiv(fea9, mu = muvar,
                                 sampleEvery = 0.03,
                                 keepEvery = 5),
                  "Individual OncoSimul trajectory",
                  fixed = TRUE)
    fea8 <- allFitnessEffects(noIntGenes =
                                  c("m" = 0.1,
                                    "D" = 0.1,
                                    "z" = 0.1,
                                    "e" = 0.1,
                                    "U" = 0.1
                                    ))
    ## OncoSimulR:::allNamedGenes(fea8)
    muvar2 <- c("U" = 1e-5, "z" = 1e-5, "e" = 1e-5, "m" = 1e-5, "D" = 1e-5)
    expect_output(oncoSimulIndiv(fea8, mu = muvar2, sampleEvery = 0.03,
                                 keepEvery = 5),
                  "Individual OncoSimul trajectory",
                  fixed = TRUE)
} )

test_that("Only no-int, unnamed, fail", {
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n s04: the seed is", pseed, "\n")
    fea9 <- allFitnessEffects(noIntGenes = c(0.1, 0.3))
    OncoSimulR:::allNamedGenes(fea9)
    muvar <- c(1e-5, 1e-7)
    expect_error(oncoSimulIndiv(fea9, mu = muvar),
                 "When using per-gene mutation rates the mu vector must be named",
                 fixed = TRUE)
} )

test_that("Only one, named, fail", {
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n s05: the seed is", pseed, "\n")
    fea9 <- allFitnessEffects(noIntGenes = c("m" = 0.1))
    muvar <- c("m" = 1e-5)
    expect_error(oncoSimulIndiv(fea9, mu = muvar),
                 "There must be at least two genes (loci) in the fitness effects",
                 fixed = TRUE)
} )

test_that("Only no-int, different names, fail", {
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n s06: the seed is", pseed, "\n")
    fea9 <- allFitnessEffects(noIntGenes = c("m" = 0.1, "D" = 0.3))
    ## OncoSimulR:::allNamedGenes(fea9)
    muvar <- c("n" = 1e-5, "D" = 1e-7)
    expect_error(oncoSimulIndiv(fea9, mu = muvar),
                 "When using per-gene mutation rates, names of genes must match",
                  fixed = TRUE)
} )

test_that("Only no-int, different numbers, fail", {
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n s07: the seed is", pseed, "\n")
    fea9 <- allFitnessEffects(noIntGenes = c("m" = 0.1, "D" = 0.3))
    muvar <- c("m" = 1e-5, "D" = 1e-7, "E" = 1e-7)
    expect_error(oncoSimulIndiv(fea9, mu = muvar),
                 "When using per-gene mutation rates, there must be the same number of genes",
                 fixed = TRUE)
    fea9 <- allFitnessEffects(noIntGenes = c("m" = 0.1, "D" = 0.3, "E" = 0.1))
    muvar <- c("m" = 1e-5, "D" = 1e-7)
    expect_error(oncoSimulIndiv(fea9, mu = muvar),
                 "When using per-gene mutation rates, there must be the same number of genes",
                 fixed = TRUE)    
} )


date()
test_that("0 or negative mu not allowed", {
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n s15: the seed is", pseed, "\n")
    muvar2 <- c("U" = 0, "z" = 1e-5, "e" = 1e-5, "m" = 1e-5, "D" = 1e-5)
    ni1 <- rep(0.02, 5)
    names(ni1) <- names(muvar2)
    fe1 <- allFitnessEffects(noIntGenes = ni1)
    no <- 1e5
    expect_error(oncoSimulIndiv(fe1, mu = muvar2, onlyCancer = FALSE,
                                initSize = no,
                                finalTime = 1
                                ),
                 "At least one per-gene mutation rate is negative or less",
                 fixed = TRUE)
    muvar2 <- c("U" = 1e-70, "z" = 1e-5, "e" = 1e-5, "m" = 1e-5, "D" = 1e-5)
    ni1 <- rep(0.02, 5)
    names(ni1) <- names(muvar2)
    fe1 <- allFitnessEffects(noIntGenes = ni1)
    no <- 1e5
    expect_error(oncoSimulIndiv(fe1, mu = muvar2, onlyCancer = FALSE,
                                initSize = no,
                                finalTime = 1
                                ),
                 "At least one per-gene mutation rate is negative or less",
                 fixed = TRUE)
    muvar2 <- c("U" = 1e-4, "z" = -0.2, "e" = 1e-5, "m" = 1e-5, "D" = 1e-5)
    ni1 <- rep(0.02, 5)
    names(ni1) <- names(muvar2)
    fe1 <- allFitnessEffects(noIntGenes = ni1)
    no <- 1e5
    expect_error(oncoSimulIndiv(fe1, mu = muvar2, onlyCancer = FALSE,
                                initSize = no,
                                finalTime = 1
                                ),
                 "(at least one) mutation rate (mu) is negative",
                 fixed = TRUE)
})
date()



date()
test_that("Same freqs, chisq, when s=0", {
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n s08: the seed is", pseed, "\n")
    muvar2 <- c("U" = 1e-5, "z" = 1e-5, "e" = 1e-5, "m" = 1e-5, "D" = 1e-5)
    ni1 <- rep(0, 5)
    names(ni1) <- names(muvar2)
    fe1 <- allFitnessEffects(noIntGenes = ni1)
    no <- 1e7
    reps <- 200
    bb <- oncoSimulPop(reps,
                       fe1, mu = muvar2,
                       onlyCancer = FALSE,
                       initSize = no,
                       finalTime = 0.001,
                       seed = NULL, mc.cores = 2
                       )
    (expectedC <- no*reps*muvar2)
    ## colSums(OncoSimulR:::geneCounts(bb))
    ## chisq.test(colSums(OncoSimulR:::geneCounts(bb)))
        ## It will fail with prob ~ p.fail
    p.fail <- 1e-3
    expect_true(chisq.test(colSums(OncoSimulR:::geneCounts(bb)),
                           p = expectedC/sum(expectedC))$p.value > p.fail)
})
date()

date()
test_that("Same freqs, chisq, when s", {
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n s09: the seed is", pseed, "\n")
    muvar2 <- c("U" = 1e-5, "z" = 1e-5, "e" = 1e-5, "m" = 1e-5, "D" = 1e-5)
    ni1 <- rep(0.02, 5)
    names(ni1) <- names(muvar2)
    fe1 <- allFitnessEffects(noIntGenes = ni1)
    no <- 1e7
    reps <- 200
    bb <- oncoSimulPop(reps,
                       fe1, mu = muvar2, onlyCancer = FALSE,
                       initSize = no,
                       finalTime = 0.001,
                       seed = NULL, mc.cores = 2
                       )
    (expectedC <- no*reps*muvar2)
    ## colSums(OncoSimulR:::geneCounts(bb))
    p.fail <- 1e-3
    expect_true(chisq.test(colSums(OncoSimulR:::geneCounts(bb)),
                           p = expectedC/sum(expectedC))$p.value > p.fail)
})
date()

date()
test_that("Different freqs as they should be ordered and chisq, when s=0",{
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n s010: the seed is", pseed, "\n")
    ## muvar2 <- c("U" = 1e-3, "z" = 1e-7, "e" = 1e-6, "m" = 1e-5, "D" = 1e-4)
    ## too small mu: observed often 0 and expected below 1. Make larger
    muvar2 <- c("U" = 1e-3, "z" = 5e-5, "e" = 5e-4, "m" = 5e-3, "D" = 1e-4)
    ni1 <- rep(0, 5)
    names(ni1) <- names(muvar2)
    fe1 <- allFitnessEffects(noIntGenes = ni1)
    no <- 5e7
    reps <- 800
    bb <- oncoSimulPop(reps,
                       fe1, mu = muvar2, onlyCancer = FALSE,
                       initSize = no,
                       finalTime = .00001,
                       seed =NULL, mc.cores = 2
                       )
    (expectedC <- no*reps*muvar2)
    colSums(OncoSimulR:::geneCounts(bb))
    ## It will fail with prob ~ p.fail
    p.fail <- 1e-3
    expect_true(chisq.test(colSums(OncoSimulR:::geneCounts(bb)),
                           p = expectedC/sum(expectedC))$p.value > p.fail)
    ## the ordering criterion could fail even if things are OK
    expect_equal(
        order(colSums(OncoSimulR:::geneCounts(bb))),
        order(expectedC))    
})
date()

date()
test_that("Different freqs as they should be ordered, when s and t > 1", {
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n s011: the seed is", pseed, "\n")
    muvar2 <- c("U" = 1e-3, "z" = 1e-7, "e" = 1e-6, "m" = 1e-5, "D" = 1e-4)
    ni2 <- rep(0.01, 5)
    names(ni2) <- names(muvar2)
    fe1 <- allFitnessEffects(noIntGenes = ni2)
    no <- 5e5
    reps <- 100
    bb <- oncoSimulPop(reps,
                       fe1, mu = muvar2, onlyCancer = FALSE,
                       initSize = no,
                       finalTime = 2,
                       mutationPropGrowth = FALSE, ## cleaner, though no real effect
                       seed = NULL, mc.cores = 2
                       )
    (expectedC <- no*reps*muvar2)
    colSums(OncoSimulR:::geneCounts(bb))
    ## A chisq will not work as we increase finalTime. But ordering of
    ## freqs. should.
    expect_equal(
        order(colSums(OncoSimulR:::geneCounts(bb))),
        order(expectedC))
})
date()

date()
test_that("Different freqs as they should be ordered, when s and t> 1, again", {
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n s012: the seed is", pseed, "\n")
    muvar2 <- c("U" = 1e-3, "z" = 5e-5, "e" = 5e-4, "m" = 5e-3, "D" = 1e-4)
    ## we increase s and finalTime
    ni2 <- rep(0.1, 5)
    names(ni2) <- names(muvar2)
    fe1 <- allFitnessEffects(noIntGenes = ni2)
    no <- 1e5
    reps <- 100
    bb <- oncoSimulPop(reps,
                       fe1, mu = muvar2, onlyCancer = FALSE,
                       initSize = no,
                       finalTime = 4,
                       mutationPropGrowth = FALSE, ## cleaner, tough no real effect
                       seed = NULL, mc.cores = 2
                       )
    (expectedC <- no*reps*muvar2)
    colSums(OncoSimulR:::geneCounts(bb))
    expect_equal(
        order(colSums(OncoSimulR:::geneCounts(bb))),
        order(expectedC))
})
date()

date()
test_that("Complex fitness specification, s diffs, tiny finalTime, systematic mu", {
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n s1: the seed is", pseed, "\n")
    p4 <- data.frame(parent = c(rep("Root", 4), "A", "B", "D", "E", "C", "F"),
                     child = c("A", "B", "D", "E", "C", "C", "F", "F", "G", "G"),
                     s = c(0.1, 0.2, 0.3, 0.4, 0.1, 0.1, 0.2, 0.2, 0.3, 0.3),
                     sh = c(rep(0, 4), c(-.9, -.9), c(-.95, -.95), c(-.99, -.99)),
                     typeDep = c(rep("--", 4), 
                                 "XMPN", "XMPN", "MN", "MN", "SM", "SM"))
    oe <- c("C > F" = -0.1, "H > I" = 0.12)
    sm <- c("I:J"  = -.1)
    sv <- c("-K:M" = -.5, "K:-M" = -.5)
    epist <- c(sm, sv)
    modules <- c("Root" = "Root", "A" = "a1",
                 "B" = "b1, b2", "C" = "c1",
                 "D" = "d1, d2", "E" = "e1",
                 "F" = "f1, f2", "G" = "g1",
                 "H" = "h1, h2", "I" = "i1",
                 "J" = "j1, j2", "K" = "k1, k2", "M" = "m1")
    noint <- runif(5, min = 0.051, max = 0.1)
    names(noint) <- paste0("n", 1:5)
    fea <- allFitnessEffects(rT = p4, epistasis = epist, orderEffects = oe,
                             noIntGenes = noint, geneToModule = modules)
    nfea <- OncoSimulR:::allNamedGenes(fea)$Gene
    ## systematic spacing
    muvar <- sample(seq(from = 5e-6, to = 1e-3, length.out = length(nfea)))
    names(muvar) <- nfea
    no <- 5e7
    reps <- 300
    bb <- oncoSimulPop(reps,
                       fea, mu = muvar,
                       onlyCancer = FALSE,
                       initSize = no,
                       finalTime = 0.0001,
                       seed = NULL, mc.cores = 2
                       )
    (expectedC <- no*reps*muvar)
    colSums(OncoSimulR:::geneCounts(bb))
    ## expectedC - colSums(OncoSimulR:::geneCounts(bb))
    ## chisq.test(colSums(OncoSimulR:::geneCounts(bb)),
    ##                        p = expectedC/sum(expectedC))
    p.fail <- 1e-3
    expect_true(chisq.test(colSums(OncoSimulR:::geneCounts(bb)),
                           p = expectedC/sum(expectedC))$p.value > p.fail)
    ## Order will not necessarily be OK, because of random fluctuations,
    ## when muvar diffs are tiny
    ## expect_equal(
    ##     order(colSums(OncoSimulR:::geneCounts(bb))),
    ##     order(expectedC))
})
date()



date()
test_that("Complex fitness specification, tiny s diffs", {
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n s2: the seed is", pseed, "\n")
    p4 <- data.frame(parent = c(rep("Root", 4), "A", "B", "D", "E", "C", "F"),
                     child = c("A", "B", "D", "E", "C", "C", "F", "F", "G", "G"),
                     s = c(0.00001, 0.00002, 0.00003, 0.00004, 0.00001, 0.00001, 0.00002, 0.00002, 0.00003, 0.00003),
                     sh = c(rep(0, 4), c(-.0000009, -.0000009), c(-.00000095, -.00000095), c(-.00000099, -.00000099)),
                     typeDep = c(rep("--", 4), 
                                 "XMPN", "XMPN", "MN", "MN", "SM", "SM"))
    oe <- c("C > F" = -0.00001, "H > I" = 0.000012)
    sm <- c("I:J"  = -.00001)
    sv <- c("-K:M" = -.000005, "K:-M" = -.000005)
    epist <- c(sm, sv)
    modules <- c("Root" = "Root", "A" = "a1",
                 "B" = "b1, b2", "C" = "c1",
                 "D" = "d1, d2", "E" = "e1",
                 "F" = "f1, f2", "G" = "g1",
                 "H" = "h1, h2", "I" = "i1",
                 "J" = "j1, j2", "K" = "k1, k2", "M" = "m1")
    noint <- runif(5, min = 0.0000051, max = 0.00001)
    names(noint) <- paste0("n", 1:5)
    fea <- allFitnessEffects(rT = p4, epistasis = epist, orderEffects = oe,
                             noIntGenes = noint, geneToModule = modules)
    nfea <- OncoSimulR:::allNamedGenes(fea)$Gene
    ## Now, random muvar
    ## muvar <- sample(seq(from = 5e-6, to = 1e-3, length.out = length(nfea)))
    muvar <- runif(length(nfea), min = 5e-6, max = 1e-3)
    names(muvar) <- nfea
    no <- 5e7
    reps <- 100
    bb <- oncoSimulPop(reps,
                       fea, mu = muvar,
                       onlyCancer = FALSE,
                       initSize = no,
                       finalTime = .0001,
                       seed = NULL, mc.cores = 2
                       )
    (expectedC <- no*reps*muvar)
    colSums(OncoSimulR:::geneCounts(bb))
    ## chisq.test(colSums(OncoSimulR:::geneCounts(bb)),
    ##                        p = expectedC/sum(expectedC))
    p.fail <- 1e-3
    expect_true(chisq.test(colSums(OncoSimulR:::geneCounts(bb)),
                           p = expectedC/sum(expectedC))$p.value > p.fail)
    ## Even with systematic spacing, you need huge reps to even out the
    ## sampling effects on order. And ordering tested above several
    ## times. This is an overkill.
    ## expect_equal(
    ##     order(colSums(OncoSimulR:::geneCounts(bb))),
    ##     order(expectedC))
})
date()

#### Repeating above, but with McFL

test_that("McFL: Per-gene mutation rates with old poset format, fail", {
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n sz01: the seed is", pseed, "\n")
    data(examplePosets)
    p701 <- examplePosets[["p701"]]
    muvar <- c(rep(1e-5, 4), rep(1e-6, 3))
    names(muvar) <- letters[1:7]
    expect_error(oncoSimulIndiv(p701, mu = muvar, model = "McFL"),
                 "Per-gene mutation rates cannot be used with the old poset format")
} )

test_that("McFL: Only no-int, and sorting", {
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n sz02: the seed is", pseed, "\n")
    fea9 <- allFitnessEffects(noIntGenes = c("m" = 0.1, "D" = 0.3))
    ## OncoSimulR:::allNamedGenes(fea9)
    muvar <- c("m" = 1e-5, "D" = 1e-7)
    expect_output(oncoSimulIndiv(fea9, mu = muvar, model = "McFL",
                                 sampleEvery = 0.03,
                                 keepEvery = 5,
                                 finalTime = 20),
                  "Individual OncoSimul trajectory",
                  fixed = TRUE)
    fea8 <- allFitnessEffects(noIntGenes =
                                  c("m" = 0.1,
                                    "D" = 0.1,
                                    "z" = 0.1,
                                    "e" = 0.1,
                                    "U" = 0.1
                                    ))
    ## OncoSimulR:::allNamedGenes(fea8)
    muvar2 <- c("U" = 1e-5, "z" = 1e-5, "e" = 1e-5, "m" = 1e-5, "D" = 1e-5)
    expect_output(oncoSimulIndiv(fea8, mu = muvar2,
                                 model = "McFL",
                                 sampleEvery = 0.03,
                                 keepEvery = 5,
                                 seed = NULL,
                                 finalTime = 20),
                  "Individual OncoSimul trajectory", 
                  fixed = TRUE)
} )

test_that("McFL: Only no-int, unnamed, fail", {
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n sz03: the seed is", pseed, "\n")
    fea9 <- allFitnessEffects(noIntGenes = c(0.1, 0.3))
    OncoSimulR:::allNamedGenes(fea9)
    muvar <- c(1e-5, 1e-7)
    expect_error(oncoSimulIndiv(fea9, mu = muvar, model = "McFL",
                                finalTime = 20),
                 "When using per-gene mutation rates the mu vector must be named",
                 fixed = TRUE)
} )

test_that("McFL: Only one, named, fail", {
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n sz04: the seed is", pseed, "\n")
    fea9 <- allFitnessEffects(noIntGenes = c("m" = 0.1))
    muvar <- c("m" = 1e-5)
    expect_error(oncoSimulIndiv(fea9, mu = muvar, model = "McFL",
                                finalTime = 20),
                 "There must be at least two genes (loci) in the fitness effects",
                 fixed = TRUE)
} )

test_that("McFL: Only no-int, different names, fail", {
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n sz05: the seed is", pseed, "\n")
    fea9 <- allFitnessEffects(noIntGenes = c("m" = 0.1, "D" = 0.3))
    ## OncoSimulR:::allNamedGenes(fea9)
    muvar <- c("n" = 1e-5, "D" = 1e-7)
    expect_error(oncoSimulIndiv(fea9, mu = muvar, model = "McFL",
                                finalTime = 20),
                 "When using per-gene mutation rates, names of genes must match",
                  fixed = TRUE)
} )

test_that("McFL: Only no-int, different numbers, fail", {
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n sz06: the seed is", pseed, "\n")
    fea9 <- allFitnessEffects(noIntGenes = c("m" = 0.1, "D" = 0.3))
    muvar <- c("m" = 1e-5, "D" = 1e-7, "E" = 1e-7)
    expect_error(oncoSimulIndiv(fea9, mu = muvar, model = "McFL"),
                 "When using per-gene mutation rates, there must be the same number of genes",
                 fixed = TRUE)
    fea9 <- allFitnessEffects(noIntGenes = c("m" = 0.1, "D" = 0.3, "E" = 0.1))
    muvar <- c("m" = 1e-5, "D" = 1e-7)
    expect_error(oncoSimulIndiv(fea9, mu = muvar, model = "McFL"),
                 "When using per-gene mutation rates, there must be the same number of genes",
                 fixed = TRUE)    
} )

date()
test_that("McFL: Same freqs, chisq, when s=0", {
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n s3: the seed is", pseed, "\n")
    muvar2 <- c("U" = 1e-5, "z" = 1e-5, "e" = 1e-5, "m" = 1e-5, "D" = 1e-5)
    ni1 <- rep(0, 5)
    names(ni1) <- names(muvar2)
    fe1 <- allFitnessEffects(noIntGenes = ni1)
    no <- 5e7
    reps <- 100
    bb <- oncoSimulPop(reps,
                       fe1, mu = muvar2, onlyCancer = FALSE,
                       model = "McFL",
                       initSize = no,
                       finalTime = 0.001,
                       seed= NULL, mc.cores = 2
                       )
    (expectedC <- no*reps*muvar2)
    colSums(OncoSimulR:::geneCounts(bb))
    ## It will fail with prob ~ p.fail
    p.fail <- 1e-3
    expect_true(chisq.test(colSums(OncoSimulR:::geneCounts(bb)),
                           p = expectedC/sum(expectedC))$p.value > p.fail)
})
date()

date()
test_that("McFL: Same freqs, chisq, when s", {
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n s4: the seed is", pseed, "\n")
    muvar2 <- c("U" = 1e-5, "z" = 1e-5, "e" = 1e-5, "m" = 1e-5, "D" = 1e-5)
    ni1 <- rep(0.02, 5)
    names(ni1) <- names(muvar2)
    fe1 <- allFitnessEffects(noIntGenes = ni1)
    no <- 1e7
    reps <- 200
    bb <- oncoSimulPop(reps,
                       fe1, mu = muvar2,
                       model = "McFL",
                       onlyCancer = FALSE,
                       initSize = no,
                       finalTime = 0.001,
                       seed = NULL, mc.cores = 2
                       )
    (expectedC <- no*reps*muvar2)
    colSums(OncoSimulR:::geneCounts(bb))
    ## It will fail with prob ~ p.fail
    p.fail <- 1e-3
    expect_true(chisq.test(colSums(OncoSimulR:::geneCounts(bb)),
                           p = expectedC/sum(expectedC))$p.value > p.fail)
})
date()

date()
test_that("McFL: Different freqs as they should be ordered and chisq, when s=0", {
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n s5: the seed is", pseed, "\n")
    muvar2 <- c("U" = 1e-3, "z" = 5e-3, "e" = 1e-4, "m" = 5e-5, "D" = 5e-4)
    ni1 <- rep(0, 5)
    names(ni1) <- names(muvar2)
    fe1 <- allFitnessEffects(noIntGenes = ni1)
    no <- 1e7
    reps <- 250
    bb <- oncoSimulPop(reps,
                       fe1, mu = muvar2, onlyCancer = FALSE,
                       model = "McFL",
                       initSize = no,
                       finalTime = 0.0001,
                       seed = NULL, mc.cores = 2                       
                       )
    (expectedC <- no*reps*muvar2)
    colSums(OncoSimulR:::geneCounts(bb))
    ## It will fail with prob ~ p.fail
    p.fail <- 1e-3
    expect_true(chisq.test(colSums(OncoSimulR:::geneCounts(bb)),
                           p = expectedC/sum(expectedC))$p.value > p.fail)
    expect_equal(
        order(colSums(OncoSimulR:::geneCounts(bb))),
        order(expectedC))    
    ## chisq.test(colSums(OncoSimulR:::geneCounts(bb)),
    ##            p = expectedC/sum(expectedC))
})
date()

date()
test_that("McFL: Different freqs as they should be ordered when s and t > 1", {
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n s6: the seed is", pseed, "\n")
    muvar2 <- c("U" = 1e-3, "z" = 5e-3, "e" = 1e-4, "m" = 5e-5, "D" = 5e-4)
    ni2 <- rep(0.01, 5)
    names(ni2) <- names(muvar2)
    fe1 <- allFitnessEffects(noIntGenes = ni2)
    no <- 1e5
    reps <- 50
    bb <- oncoSimulPop(reps,
                       fe1, mu = muvar2, onlyCancer = FALSE,
                       model = "McFL",
                       initSize = no,
                       finalTime = 4,
                       mutationPropGrowth = FALSE, ## cleaner, though no real effect
                       seed = NULL, mc.cores = 2
                       )
    (expectedC <- no*reps*muvar2)
    colSums(OncoSimulR:::geneCounts(bb))
    expect_equal(
        order(colSums(OncoSimulR:::geneCounts(bb))),
        order(expectedC))
})
date()

date()
test_that("McFL: Different freqs as they should be ordered when s and t > 1, again", {
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n s7: the seed is", pseed, "\n")
    ## Increase s and time
    muvar2 <- c("U" = 1e-3, "z" = 5e-3, "e" = 1e-4, "m" = 5e-5, "D" = 5e-4)
    ni2 <- rep(0.2, 5)
    names(ni2) <- names(muvar2)
    fe1 <- allFitnessEffects(noIntGenes = ni2)
    no <- 1e5
    reps <- 20
    bb <- oncoSimulPop(reps,
                       fe1, mu = muvar2, onlyCancer = FALSE,
                       model = "McFL",
                       initSize = no,
                       finalTime = 10,
                       mutationPropGrowth = FALSE, ## cleaner, though no real effect
                       seed = NULL, mc.cores = 2
                       )
    (expectedC <- no*reps*muvar2)
    colSums(OncoSimulR:::geneCounts(bb))
    expect_equal(
        order(colSums(OncoSimulR:::geneCounts(bb))),
        order(expectedC))
})
date()

date()
test_that("McFL: Complex fitness specification, s diffs, tiny finalTime, systematic mu", {
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n s8: the seed is", pseed, "\n")
    p4 <- data.frame(parent = c(rep("Root", 4), "A", "B", "D", "E", "C", "F"),
                     child = c("A", "B", "D", "E", "C", "C", "F", "F", "G", "G"),
                     s = c(0.1, 0.2, 0.3, 0.4, 0.1, 0.1, 0.2, 0.2, 0.3, 0.3),
                     sh = c(rep(0, 4), c(-.9, -.9), c(-.95, -.95), c(-.99, -.99)),
                     typeDep = c(rep("--", 4), 
                                 "XMPN", "XMPN", "MN", "MN", "SM", "SM"))
    oe <- c("C > F" = -0.1, "H > I" = 0.12)
    sm <- c("I:J"  = -.1)
    sv <- c("-K:M" = -.5, "K:-M" = -.5)
    epist <- c(sm, sv)
    modules <- c("Root" = "Root", "A" = "a1",
                 "B" = "b1, b2", "C" = "c1",
                 "D" = "d1, d2", "E" = "e1",
                 "F" = "f1, f2", "G" = "g1",
                 "H" = "h1, h2", "I" = "i1",
                 "J" = "j1, j2", "K" = "k1, k2", "M" = "m1")
    noint <- runif(5, min = 0.051, max = 0.1)
    names(noint) <- paste0("n", 1:5)
    fea <- allFitnessEffects(rT = p4, epistasis = epist, orderEffects = oe,
                             noIntGenes = noint, geneToModule = modules)
    nfea <- OncoSimulR:::allNamedGenes(fea)$Gene
    ## systematic spacing
    muvar <- sample(seq(from = 5e-6, to = 1e-3, length.out = length(nfea)))
    names(muvar) <- nfea
    no <- 5e7
    reps <- 300
    bb <- oncoSimulPop(reps,
                       fea, mu = muvar,
                       onlyCancer = FALSE,
                       initSize = no,
                       finalTime = 0.0001,
                       seed = NULL, mc.cores = 2,
                       model = "McFL"
                       )
    (expectedC <- no*reps*muvar)
    colSums(OncoSimulR:::geneCounts(bb))
    ## expectedC - colSums(OncoSimulR:::geneCounts(bb))
    ## chisq.test(colSums(OncoSimulR:::geneCounts(bb)),
    ##                        p = expectedC/sum(expectedC))
    p.fail <- 1e-3
    expect_true(chisq.test(colSums(OncoSimulR:::geneCounts(bb)),
                           p = expectedC/sum(expectedC))$p.value > p.fail)
})
date()

date()
test_that("McFL:Complex fitness specification, tiny s diffs", {
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n s9: the seed is", pseed, "\n")
    p4 <- data.frame(parent = c(rep("Root", 4), "A", "B", "D", "E", "C", "F"),
                     child = c("A", "B", "D", "E", "C", "C", "F", "F", "G", "G"),
                     s = c(0.00001, 0.00002, 0.00003, 0.00004, 0.00001, 0.00001, 0.00002, 0.00002, 0.00003, 0.00003),
                     sh = c(rep(0, 4), c(-.0000009, -.0000009), c(-.00000095, -.00000095), c(-.00000099, -.00000099)),
                     typeDep = c(rep("--", 4), 
                                 "XMPN", "XMPN", "MN", "MN", "SM", "SM"))
    oe <- c("C > F" = -0.00001, "H > I" = 0.000012)
    sm <- c("I:J"  = -.00001)
    sv <- c("-K:M" = -.000005, "K:-M" = -.000005)
    epist <- c(sm, sv)
    modules <- c("Root" = "Root", "A" = "a1",
                 "B" = "b1, b2", "C" = "c1",
                 "D" = "d1, d2", "E" = "e1",
                 "F" = "f1, f2", "G" = "g1",
                 "H" = "h1, h2", "I" = "i1",
                 "J" = "j1, j2", "K" = "k1, k2", "M" = "m1")
    noint <- runif(5, min = 0.0000051, max = 0.00001)
    names(noint) <- paste0("n", 1:5)
    fea <- allFitnessEffects(rT = p4, epistasis = epist, orderEffects = oe,
                             noIntGenes = noint, geneToModule = modules)
    nfea <- OncoSimulR:::allNamedGenes(fea)$Gene
    ## Now, random muvar
    ## muvar <- sample(seq(from = 5e-6, to = 1e-3, length.out = length(nfea)))
    muvar <- runif(length(nfea), min = 5e-6, max = 1e-3)
    names(muvar) <- nfea
    no <- 5e7
    reps <- 100
    bb <- oncoSimulPop(reps,
                       fea, mu = muvar,
                       onlyCancer = FALSE,
                       initSize = no,
                       finalTime = .0001,
                       model = "McFL",
                       seed = NULL, mc.cores = 2
                       )
    (expectedC <- no*reps*muvar)
    colSums(OncoSimulR:::geneCounts(bb))
    ## chisq.test(colSums(OncoSimulR:::geneCounts(bb)),
    ##                        p = expectedC/sum(expectedC))
    p.fail <- 1e-3
    expect_true(chisq.test(colSums(OncoSimulR:::geneCounts(bb)),
                           p = expectedC/sum(expectedC))$p.value > p.fail)
    ## Even with systematic spacing, you need huge reps to even out the
    ## sampling effects on order. And ordering tested above several
    ## times. This is an overkill.
    ## expect_equal(
    ##     order(colSums(OncoSimulR:::geneCounts(bb))),
    ##     order(expectedC))
})
date()


date()
test_that("get.gene.counts exercising for NA case", {
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n s10: the seed is", pseed, "\n")
    muvar2 <- c("U" = 1e-3, "z" = 1e-7, "e" = 1e-6, "m" = 1e-5, "D" = 1e-4)
    ni1 <- rep(0, 5)
    names(ni1) <- names(muvar2)
    fe1 <- allFitnessEffects(noIntGenes = ni1)
    ou1 <- oncoSimulIndiv(fe1, mu = muvar2,
                          initSize = 100,
                          onlyCancer = FALSE,
                          seed = NULL)
    expect_output(OncoSimulR:::get.gene.counts(ou1),
                  "$counts",
                  fixed = TRUE)
    expect_output(OncoSimulR:::geneCounts(ou1),
                  "0",
                  fixed = TRUE)
})






date()
test_that("Init mutant  with tiny mutation always present", {
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n s11: the seed is", pseed, "\n")
    ## We check the two init mutants have same frequency, and it is the
    ## largest. So we also verify no additional mutation to
    ## either one (since their freq is same) and we verify present in all.
    p4 <- data.frame(parent = c(rep("Root", 4), "A", "B", "D", "E", "C", "F"),
                     child = c("A", "B", "D", "E", "C", "C", "F", "F", "G", "G"),
                     s = c(0.01, 0.02, 0.03, 0.04, 0.1, 0.1, 0.2, 0.2, 0.3, 0.3),
                     sh = c(rep(0, 4), c(-.9, -.9), c(-.95, -.95), c(-.99, -.99)),
                     typeDep = c(rep("--", 4), 
                                 "XMPN", "XMPN", "MN", "MN", "SM", "SM"))
    oe <- c("C > F" = -0.1, "H > I" = 0.12)
    sm <- c("I:J"  = -1)
    sv <- c("-K:M" = -.5, "K:-M" = -.5)
    epist <- c(sm, sv)
    modules <- c("Root" = "Root", "A" = "a1",
                 "B" = "b1, b2", "C" = "c1",
                 "D" = "d1, d2", "E" = "e1",
                 "F" = "f1, f2", "G" = "g1",
                 "H" = "h1, h2", "I" = "i1",
                 "J" = "j1, j2", "K" = "k1, k2", "M" = "m1")
    noint <- runif(5, min = 0.01, max = 0.1)
    names(noint) <- paste0("n", 1:5)
    fea <- allFitnessEffects(rT = p4, epistasis = epist, orderEffects = oe,
                             noIntGenes = noint, geneToModule = modules)
    nfea <- OncoSimulR:::allNamedGenes(fea)$Gene
    ## muvar <- runif(length(nfea), min = 1e-7, max = 1e-3) ## too tiny
    ## diffs sometimes for order comp
    muvar <- sample(seq(from = 1e-7, to = 1e-4, length.out = length(nfea)))
    names(muvar) <- nfea
    muvar["h2"] <- 1e-13
    muvar["i1"] <- 2e-13
    no <- 1e5
    reps <- 20
    bb <- oncoSimulPop(reps,
                       fea, mu = muvar,
                       onlyCancer = FALSE,
                       initSize = no,
                       sampleEvery = 0.03,
                       keepEvery = 2,
                       finalTime = 10000, ## yes, huge; we only get close to 4 or 5.
                       detectionDrivers = 4,
                       mutationPropGrowth = FALSE, ## yes, exclude this possible effect
                       initMutant = "h2 > i1",
                       mc.cores = 2
                       )
    (expectedC <- no*reps*muvar)
    ccs <- colSums(OncoSimulR:::geneCounts(bb))
    totalindivs <- sum(unlist(lapply(bb, function(x) x$TotalPopSize)))
    expect_true(ccs["h2"] == ccs["i1"])
    expect_true(ccs["h2"] == totalindivs)
    expect_true(all(ccs["h2"] > ccs[!(names(ccs) %in% c("h2", "i1"))]))
    ## of course, no agreement with chi-square or ordering
    p.fail <- 1e-6
    expect_true(chisq.test(colSums(OncoSimulR:::geneCounts(bb)),
                           p = expectedC/sum(expectedC))$p.value < p.fail)
    expect_false(
        identical(
            order(colSums(OncoSimulR:::geneCounts(bb))),
            order(expectedC)))
})
date()

date()
test_that("McFL: Init mutant with tiny mutation always present", {
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n s12: the seed is", pseed, "\n")
    p4 <- data.frame(parent = c(rep("Root", 4), "A", "B", "D", "E", "C", "F"),
                     child = c("A", "B", "D", "E", "C", "C", "F", "F", "G", "G"),
                     s = c(0.01, 0.02, 0.03, 0.04, 0.1, 0.1, 0.2, 0.2, 0.3, 0.3),
                     sh = c(rep(0, 4), c(-.9, -.9), c(-.95, -.95), c(-.99, -.99)),
                     typeDep = c(rep("--", 4), 
                                 "XMPN", "XMPN", "MN", "MN", "SM", "SM"))
    oe <- c("C > F" = -0.1, "H > I" = 0.12)
    sm <- c("I:J"  = -1)
    sv <- c("-K:M" = -.5, "K:-M" = -.5)
    epist <- c(sm, sv)
    modules <- c("Root" = "Root", "A" = "a1",
                 "B" = "b1, b2", "C" = "c1",
                 "D" = "d1, d2", "E" = "e1",
                 "F" = "f1, f2", "G" = "g1",
                 "H" = "h1, h2", "I" = "i1",
                 "J" = "j1, j2", "K" = "k1, k2", "M" = "m1")
    noint <- runif(5, min = 0.01, max = 0.1)
    names(noint) <- paste0("n", 1:5)
    fea <- allFitnessEffects(rT = p4, epistasis = epist, orderEffects = oe,
                             noIntGenes = noint, geneToModule = modules)
    nfea <- OncoSimulR:::allNamedGenes(fea)$Gene
    ## muvar <- runif(length(nfea), min = 1e-7, max = 1e-3) ## too tiny
    ## diffs sometimes for order comp
    muvar <- sample(seq(from = 1e-7, to = 1e-5, length.out = length(nfea)))
    names(muvar) <- nfea
    muvar["h2"] <- 3e-13
    muvar["i1"] <- 1e-13
    no <- 1e4
    reps <- 10
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n s12b: the seed is", pseed, "\n")
    bb <- oncoSimulPop(5, ##reps,
                       fea, mu = muvar,
                       onlyCancer = FALSE,
                       initSize = no,
                       model = "McFL",
                       sampleEvery = 0.02,
                       keepEvery = 2,
                       finalTime = 300,
                       mutationPropGrowth = FALSE, ## yes, exclude this possible effect
                       initMutant = "h2 > i1",
                       mc.cores = 2
                       )
    (expectedC <- no*reps*muvar)
    ccs <- colSums(OncoSimulR:::geneCounts(bb))
    totalindivs <- sum(unlist(lapply(bb, function(x) x$TotalPopSize)))
    expect_true(ccs["h2"] == ccs["i1"])
    expect_true(ccs["h2"] == totalindivs)
    expect_true(all(ccs["h2"] > ccs[!(names(ccs) %in% c("h2", "i1"))]))
    ## this will occasionally fail
    p.fail <- 1e-6
    expect_true(chisq.test(colSums(OncoSimulR:::geneCounts(bb)),
                           p = expectedC/sum(expectedC))$p.value < p.fail)
    expect_false(
        identical(
            order(colSums(OncoSimulR:::geneCounts(bb))),
            order(expectedC)))
})
date()

date()
test_that("Different freqs as they should be ordered and chisq, when s  and a tiny mu", {
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n s13: the seed is", pseed, "\n")
    muvar2 <- c("U" = 1e-13, "z" = 1e-7, "e" = 1e-6, "m" = 1e-5, "D" = 1e-4)
    ni2 <- rep(0.01, 5)
    names(ni2) <- names(muvar2)
    ni2["U"] <- 0.5
    fe1 <- allFitnessEffects(noIntGenes = ni2)
    no <- 1e6
    reps <- 100
    bb <- oncoSimulPop(reps,
                       fe1, mu = muvar2, onlyCancer = FALSE,
                       initSize = no,
                       finalTime = 5,
                       mutationPropGrowth = FALSE, ## cleaner, though no real effect
                       seed = NULL, mc.cores = 2
                       )
    (expectedC <- no*reps*muvar2)
    colSums(OncoSimulR:::geneCounts(bb))
    expect_true(colSums(OncoSimulR:::geneCounts(bb))[1] == 0)
    expect_equal(
        order(colSums(OncoSimulR:::geneCounts(bb))),
        order(expectedC))
    ## A chisq will not work as we increase finalTime.
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n s13b: the seed is", pseed, "\n")
    bb <- oncoSimulPop(reps,
                       fe1, mu = muvar2, onlyCancer = FALSE,
                       initSize = no,
                       finalTime = .001,
                       mutationPropGrowth = FALSE, ## cleaner, though no real effect
                       seed = NULL, mc.cores = 2
                       )
    (expectedC <- no*reps*muvar2)
    colSums(OncoSimulR:::geneCounts(bb))
    expect_true(colSums(OncoSimulR:::geneCounts(bb))[1] == 0)
    ## This will fail sometimes
    p.fail <- 1e-3
    expect_true(chisq.test(colSums(OncoSimulR:::geneCounts(bb))[-1],
                           p = expectedC[-1]/sum(expectedC))$p.value > p.fail)
    ## could fail because very small freqs
    ## expect_equal(
    ##     order(colSums(OncoSimulR:::geneCounts(bb))),
    ##     order(expectedC))
}
)
date()

date()
test_that("McFL: Different freqs as they should be ordered and chisq, when s  and a tiny mu", {
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n s14: the seed is", pseed, "\n")
    muvar2 <- c("U" = 1e-13, "z" = 5e-7, "e" = 1e-6, "m" = 1e-5, "D" = 1e-4)
    ni2 <- rep(0.01, 5)
    names(ni2) <- names(muvar2)
    ni2["U"] <- 0.5
    fe1 <- allFitnessEffects(noIntGenes = ni2)
    no <- 1e5
    reps <- 400
    ## Beware: with McFL final time has to be small or we will see
    ## mutations in U as it is the only one left. 
    bb <- oncoSimulPop(reps,
                       fe1, mu = muvar2, onlyCancer = FALSE,
                       initSize = no,
                       model = "McFL",
                       finalTime = 50,
                       mutationPropGrowth = FALSE, ## cleaner, though no real effect
                       seed = NULL, mc.cores = 2
                       )
    (expectedC <- no*reps*muvar2)
    colSums(OncoSimulR:::geneCounts(bb))
    expect_true(colSums(OncoSimulR:::geneCounts(bb))[1] == 0)
    expect_equal(
        order(colSums(OncoSimulR:::geneCounts(bb))),
        order(expectedC))
    ## A chisq will not work as we increase finalTime.
    no <- 1e7
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n s14b: the seed is", pseed, "\n")
    bb <- oncoSimulPop(reps,
                       fe1, mu = muvar2, onlyCancer = FALSE,
                       initSize = no,
                       model = "McFL",
                       finalTime = .001,
                       mutationPropGrowth = FALSE, ## cleaner, though no real effect
                       seed = NULL, mc.cores = 2
                       )
    (expectedC <- no*reps*muvar2)
    colSums(OncoSimulR:::geneCounts(bb))
    expect_true(colSums(OncoSimulR:::geneCounts(bb))[1] == 0)
    ## This will fail sometimes
    p.fail <- 1e-3
    expect_true(chisq.test(colSums(OncoSimulR:::geneCounts(bb))[-1],
                           p = expectedC[-1]/sum(expectedC))$p.value > p.fail)
    ## expect_equal(
    ##     order(colSums(OncoSimulR:::geneCounts(bb))),
    ##     order(expectedC))
}
)
date()


date()
test_that("Different freqs as they should be ordered and chisq, when s=0, and initMutant",{
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n s17: the seed is", pseed, "\n")
    muvar2 <- c("U" = 1e-3, "z" = 5e-5, "e" = 1e-6, "m" = 1e-5, "D" = 1e-4)
    ni1 <- rep(0, 5)
    names(ni1) <- names(muvar2)
    fe1 <- allFitnessEffects(noIntGenes = ni1)
    no <- 1e7
    reps <- 200
    bb <- oncoSimulPop(reps,
                       fe1, mu = muvar2, onlyCancer = FALSE,
                       initSize = no,
                       finalTime = .001,
                       seed =NULL,
                       initMutant = "e",
                       mc.cores = 2
                       )
    (expectedC <- no*reps*muvar2)
    ccs <- colSums(OncoSimulR:::geneCounts(bb))
    totalindivs <- sum(unlist(lapply(bb, function(x) x$TotalPopSize)))
    expect_true(ccs["e"] == totalindivs)
    ## It will fail with prob ~ p.fail
    p.fail <- 1e-3
    expect_true(chisq.test(colSums(OncoSimulR:::geneCounts(bb))[-3],
                           p = expectedC[-3]/sum(expectedC[-3]))$p.value > p.fail)
    expect_equal(
        order(colSums(OncoSimulR:::geneCounts(bb))[-3]),
        order(expectedC[-3]))
})
date()

date()
test_that("McFL: Different freqs as they should be ordered and chisq, when s=0, and initMutant",{
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n s18: the seed is", pseed, "\n")
    muvar2 <- c("U" = 1e-3, "z" = 5e-5, "e" = 1e-6, "m" = 1e-5, "D" = 1e-4)
    ni1 <- rep(0, 5)
    names(ni1) <- names(muvar2)
    fe1 <- allFitnessEffects(noIntGenes = ni1)
    no <- 1e7
    reps <- 200
    bb <- oncoSimulPop(reps,
                       fe1, mu = muvar2, onlyCancer = FALSE,
                       initSize = no,
                       finalTime = 0.001,
                       seed =NULL,
                       model = "McFL",
                       initMutant = "m",
                       mc.cores = 2
                       )
    (expectedC <- no*reps*muvar2)
    ccs <- colSums(OncoSimulR:::geneCounts(bb))
    totalindivs <- sum(unlist(lapply(bb, function(x) x$TotalPopSize)))
    expect_true(ccs["m"] == totalindivs)
    ## It will fail with prob ~ p.fail
    p.fail <- 1e-3
    expect_true(chisq.test(colSums(OncoSimulR:::geneCounts(bb))[-4],
                           p = expectedC[-4]/sum(expectedC[-4]))$p.value > p.fail)
    expect_equal(
        order(colSums(OncoSimulR:::geneCounts(bb))[-4]),
        order(expectedC[-4]))
    ## chisq.test(colSums(OncoSimulR:::geneCounts(bb)),
    ##            p = expectedC/sum(expectedC))
})
date()

date()
test_that("Different freqs as they are expected with chisq, when s=0 and initMutant, many genotypes",{
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n s19: the seed is", pseed, "\n")
    ft <- 0.001 ## yes, small
    lni <- 100  ## 16
    muvar2 <- runif(lni, min = 1e-4, max = 1e-3)
    names(muvar2) <- c(replicate(lni,
                                 paste(sample(letters, 12), collapse = "")))
    names(muvar2)[3] <- "e"
    muvar2[3] <- 1e-9
    ni1 <- rep(0, lni)
    names(ni1) <- names(muvar2)
    fe1 <- allFitnessEffects(noIntGenes = ni1)
    no <- 1e7
    reps <- 100
    bb <- oncoSimulPop(reps,
                       fe1, mu = muvar2, onlyCancer = FALSE,
                       initSize = no,
                       finalTime = ft,
                       seed =NULL,
                       initMutant = "e",
                       mc.cores = 2
                       )
    (expectedC <- no*reps*muvar2)
    (ccs <- colSums(OncoSimulR:::geneCounts(bb)))
    totalindivs <- sum(unlist(lapply(bb, function(x) x$TotalPopSize)))
    expect_true(ccs["e"] == totalindivs)
    ## It will fail with prob ~ p.fail
    p.fail <- 1e-3
    expect_true(chisq.test(colSums(OncoSimulR:::geneCounts(bb))[-3],
                           p = expectedC[-3]/sum(expectedC[-3]))$p.value > p.fail)
})
date()


date()
test_that("McFL: Different freqs as they are expected with chisq, when s=0 and initMutant, many genotypes",{
    
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n mcs19: the seed is", pseed, "\n")
    ft <- 0.001 ## yes, small
    lni <- 100  ## 16
    muvar2 <- runif(lni, min = 1e-4, max = 1e-3)
    names(muvar2) <- c(replicate(lni,
                                 paste(sample(letters, 12), collapse = "")))
    names(muvar2)[3] <- "e"
    muvar2[3] <- 1e-9
    ni1 <- rep(0, lni)
    names(ni1) <- names(muvar2)
    fe1 <- allFitnessEffects(noIntGenes = ni1)
    no <- 1e7
    reps <- 100
    bb <- oncoSimulPop(reps,
                       fe1, mu = muvar2, onlyCancer = FALSE,
                       initSize = no,
                       finalTime = ft,
                       seed =NULL,
                       initMutant = "e",
                       model = "McFL",
                       mc.cores = 2
                       )
    (expectedC <- no*reps*muvar2)
    (ccs <- colSums(OncoSimulR:::geneCounts(bb)))
    totalindivs <- sum(unlist(lapply(bb, function(x) x$TotalPopSize)))
    expect_true(ccs["e"] == totalindivs)
    ## It will fail with prob ~ p.fail
    p.fail <- 1e-3
    expect_true(chisq.test(colSums(OncoSimulR:::geneCounts(bb))[-3],
                           p = expectedC[-3]/sum(expectedC[-3]))$p.value > p.fail)
    
})
date()


date()
test_that("Different freqs as they should be ordered and chisq, when s=0, and initMutant",{
    ## More on the above, with less variation. But yet another set of tests.
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n s20: the seed is", pseed, "\n")
    muvar2 <- c("U" = 1e-3, "z" = 1e-7, "e" = 1e-6, "m" = 1e-5, "D" = 1e-4)
    ## moderately small mu
    muvar2[] <- 1e-5
    muvar2["e"] <- 1e-3
    ni1 <- rep(0, 5)
    names(ni1) <- names(muvar2)
    fe1 <- allFitnessEffects(noIntGenes = ni1)
    no <- 1e7
    reps <- 200
    bb <- oncoSimulPop(reps,
                       fe1, mu = muvar2, onlyCancer = FALSE,
                       initSize = no,
                       finalTime = 0.002,
                       mutationPropGrowth = FALSE, ## cleaner, though no real effect
                       seed =NULL,
                       initMutant = "e",
                       mc.cores = 2
                       )
    (expectedC <- no*reps*muvar2)
    ccs <- colSums(OncoSimulR:::geneCounts(bb))
    totalindivs <- sum(unlist(lapply(bb, function(x) x$TotalPopSize)))
    expect_true(ccs["e"] == totalindivs)
    ## It will fail with prob ~ p.fail
    p.fail <- 1e-3
    expect_true(chisq.test(colSums(OncoSimulR:::geneCounts(bb))[-3],
                           p = expectedC[-3]/sum(expectedC[-3]))$p.value > p.fail)
    muvar2 <- c("U" = 1e-3, "z" = 1e-7, "e" = 1e-6, "m" = 1e-5, "D" = 1e-4)
    ## relatively large mu
    muvar2[] <- 1e-3
    muvar2["e"] <- 1e-6
    ni1 <- rep(0, 5)
    names(ni1) <- names(muvar2)
    fe1 <- allFitnessEffects(noIntGenes = ni1)
    no <- 1e7
    reps <- 200
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n s20b: the seed is", pseed, "\n")
    bb <- oncoSimulPop(reps,
                       fe1, mu = muvar2, onlyCancer = FALSE,
                       initSize = no,
                       finalTime = 0.002,
                       mutationPropGrowth = FALSE, ## cleaner, though no real effect
                       seed =NULL,
                       initMutant = "e",
                       mc.cores = 2
                       )
    (expectedC <- no*reps*muvar2)
    ccs <- colSums(OncoSimulR:::geneCounts(bb))
    totalindivs <- sum(unlist(lapply(bb, function(x) x$TotalPopSize)))
    expect_true(ccs["e"] == totalindivs)
    ## It will fail with prob ~ p.fail
    p.fail <- 1e-3
    expect_true(chisq.test(colSums(OncoSimulR:::geneCounts(bb))[-3],
                           p = expectedC[-3]/sum(expectedC[-3]))$p.value > p.fail)
    ## nope, as many are equal
    ## expect_equal(
    ##     order(colSums(OncoSimulR:::geneCounts(bb))[-3]),
    ##     order(expectedC[-3]))
    ## chisq.test(colSums(OncoSimulR:::geneCounts(bb)),
    ##            p = expectedC/sum(expectedC))
    ## some moderate, one very large
    muvar2 <- c("U" = 1e-3, "z" = 1e-7, "e" = 1e-6, "m" = 1e-5, "D" = 1e-4)
    muvar2[] <- 1e-4
    muvar2["e"] <- 1e-6
    muvar2[4] <- 1e-2
    ni1 <- rep(0, 5)
    names(ni1) <- names(muvar2)
    fe1 <- allFitnessEffects(noIntGenes = ni1)
    no <- 1e7
    reps <- 200
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n s20c: the seed is", pseed, "\n")
    bb <- oncoSimulPop(reps,
                       fe1, mu = muvar2, onlyCancer = FALSE,
                       initSize = no,
                       finalTime = 0.002,
                       mutationPropGrowth = FALSE, ## cleaner, though no real effect
                       seed =NULL,
                       initMutant = "e",
                       mc.cores = 2
                       )
    (expectedC <- no*reps*muvar2)
    (ccs <- colSums(OncoSimulR:::geneCounts(bb)))
    totalindivs <- sum(unlist(lapply(bb, function(x) x$TotalPopSize)))
    expect_true(ccs["e"] == totalindivs)
    ## It will fail with prob ~ p.fail
    p.fail <- 1e-3
    expect_true(chisq.test(colSums(OncoSimulR:::geneCounts(bb))[-3],
                           p = expectedC[-3]/sum(expectedC[-3]))$p.value > p.fail)
})
date()




date()
test_that("McFL: Different freqs as they should be ordered and chisq, when s=0, and initMutant",{
    ## More on the above, with less variation. But yet another set of tests.
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n mcs20: the seed is", pseed, "\n")
    muvar2 <- c("U" = 1e-3, "z" = 1e-7, "e" = 1e-6, "m" = 1e-5, "D" = 1e-4)
    ## moderately small mu
    muvar2[] <- 1e-5
    muvar2["e"] <- 1e-3
    ni1 <- rep(0, 5)
    names(ni1) <- names(muvar2)
    fe1 <- allFitnessEffects(noIntGenes = ni1)
    no <- 1e7
    reps <- 200
    bb <- oncoSimulPop(reps,
                       fe1, mu = muvar2, onlyCancer = FALSE,
                       initSize = no,
                       finalTime = 0.002,
                       mutationPropGrowth = FALSE, ## cleaner, though no real effect
                       seed =NULL,
                       initMutant = "e",
                       model = "McFL",
                       mc.cores = 2
                       )
    (expectedC <- no*reps*muvar2)
    ccs <- colSums(OncoSimulR:::geneCounts(bb))
    totalindivs <- sum(unlist(lapply(bb, function(x) x$TotalPopSize)))
    expect_true(ccs["e"] == totalindivs)
    ## It will fail with prob ~ p.fail
    p.fail <- 1e-3
    expect_true(chisq.test(colSums(OncoSimulR:::geneCounts(bb))[-3],
                           p = expectedC[-3]/sum(expectedC[-3]))$p.value > p.fail)
    muvar2 <- c("U" = 1e-3, "z" = 1e-7, "e" = 1e-6, "m" = 1e-5, "D" = 1e-4)
    ## relatively large mu
    muvar2[] <- 1e-3
    muvar2["e"] <- 1e-6
    ni1 <- rep(0, 5)
    names(ni1) <- names(muvar2)
    fe1 <- allFitnessEffects(noIntGenes = ni1)
    no <- 1e7
    reps <- 200
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n mcs20b: the seed is", pseed, "\n")
    bb <- oncoSimulPop(reps,
                       fe1, mu = muvar2, onlyCancer = FALSE,
                       initSize = no,
                       finalTime = 0.002,
                       mutationPropGrowth = FALSE, ## cleaner, though no real effect
                       seed =NULL,
                       initMutant = "e",
                       model = "McFL",
                       mc.cores = 2
                       )
    (expectedC <- no*reps*muvar2)
    ccs <- colSums(OncoSimulR:::geneCounts(bb))
    totalindivs <- sum(unlist(lapply(bb, function(x) x$TotalPopSize)))
    expect_true(ccs["e"] == totalindivs)
    ## It will fail with prob ~ p.fail
    p.fail <- 1e-3
    expect_true(chisq.test(colSums(OncoSimulR:::geneCounts(bb))[-3],
                           p = expectedC[-3]/sum(expectedC[-3]))$p.value > p.fail)
    ## nope, as many are equal
    ## expect_equal(
    ##     order(colSums(OncoSimulR:::geneCounts(bb))[-3]),
    ##     order(expectedC[-3]))
    ## chisq.test(colSums(OncoSimulR:::geneCounts(bb)),
    ##            p = expectedC/sum(expectedC))
    ## some moderate, one very large
    muvar2 <- c("U" = 1e-3, "z" = 1e-7, "e" = 1e-6, "m" = 1e-5, "D" = 1e-4)
    muvar2[] <- 1e-4
    muvar2["e"] <- 1e-6
    muvar2[4] <- 1e-2
    ni1 <- rep(0, 5)
    names(ni1) <- names(muvar2)
    fe1 <- allFitnessEffects(noIntGenes = ni1)
    no <- 1e7
    reps <- 200
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n mcs20c: the seed is", pseed, "\n")
    bb <- oncoSimulPop(reps,
                       fe1, mu = muvar2, onlyCancer = FALSE,
                       initSize = no,
                       finalTime = 0.002,
                       mutationPropGrowth = FALSE, ## cleaner, though no real effect
                       seed =NULL,
                       initMutant = "e",
                       model = "McFL",
                       mc.cores = 2
                       )
    (expectedC <- no*reps*muvar2)
    (ccs <- colSums(OncoSimulR:::geneCounts(bb)))
    totalindivs <- sum(unlist(lapply(bb, function(x) x$TotalPopSize)))
    expect_true(ccs["e"] == totalindivs)
    ## It will fail with prob ~ p.fail
    p.fail <- 1e-3
    expect_true(chisq.test(colSums(OncoSimulR:::geneCounts(bb))[-3],
                           p = expectedC[-3]/sum(expectedC[-3]))$p.value > p.fail)
})
date()



date()
test_that("Expect freqs, num clones, muts per clone for different per-gene-mut",{
## More on the above, with less variation. But yet another set of tests.
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n df1: the seed is", pseed, "\n")
    ng <- 10
    ni <- rep(0, ng)
    m1 <- runif(ng, min = 1e-6, max = 1e-5)
    m2 <- runif(ng, min = 1e-4, max = 1e-3)
    names(ni) <- names(m1) <- names(m2) <- c(replicate(ng,
                                 paste(sample(letters, 12), collapse = "")))
    fe1 <- allFitnessEffects(noIntGenes = ni)
    ft <- 0.01
    no <- 1e7
    reps <- 80
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n df1a: the seed is", pseed, "\n")
    b1 <- oncoSimulPop(reps,
                       fe1,
                       mu = m1,
                       onlyCancer = FALSE,
                       initSize = no,
                       finalTime = ft,
                       mutationPropGrowth = FALSE, ## cleaner, though no real effect
                       seed =NULL,
                       mc.cores = 2
                       )
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n df1b: the seed is", pseed, "\n")
    b2 <- oncoSimulPop(reps,
                       fe1,
                       mu = m2,
                       onlyCancer = FALSE,
                       initSize = no,
                       finalTime = ft,
                       mutationPropGrowth = FALSE, ## cleaner, though no real effect
                       seed =NULL,
                       mc.cores = 2
                       )
    (expected1 <- no*reps*m1)
    (expected2 <- no*reps*m2)
    (cc1 <- colSums(OncoSimulR:::geneCounts(b1)))
    (cc2 <- colSums(OncoSimulR:::geneCounts(b2)))    
    ## It will fail with prob ~ p.fail
    p.fail <- 1e-3
    expect_true(chisq.test(colSums(OncoSimulR:::geneCounts(b1)),
                           p = expected1/sum(expected1))$p.value > p.fail)
    expect_true(chisq.test(colSums(OncoSimulR:::geneCounts(b2)),
                           p = expected2/sum(expected2))$p.value > p.fail)
    expect_true( median(summary(b2)$NumClones) >
                 median(summary(b1)$NumClones))
    ## Note the short time, so this is not always very different as few
    ## have double or triple mutants
    expect_true( mean(mutsPerClone(b2)) >
                 mean(mutsPerClone(b1)))
})

date()
test_that("Num clones, muts per clone for different per-gene-mut",{
    ## Like previous, but larger finalTime, so no longer chi-square test
    ## here.
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n df2: the seed is", pseed, "\n")
    ng <- 40
    ni <- rep(0, ng)
    m1 <- runif(ng, min = 1e-6, max = 1e-5)
    m2 <- runif(ng, min = 1e-4, max = 1e-3)
    names(ni) <- names(m1) <- names(m2) <- c(replicate(ng,
                                 paste(sample(letters, 12), collapse = "")))
    fe1 <- allFitnessEffects(noIntGenes = ni)
    ft <- 2
    no <- 1e5
    reps <- 20
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n df2a: the seed is", pseed, "\n")
    b1 <- oncoSimulPop(reps,
                       fe1,
                       mu = m1,
                       onlyCancer = FALSE,
                       initSize = no,
                       finalTime = ft,
                       mutationPropGrowth = FALSE, ## cleaner, though no real effect
                       seed =NULL,
                       mc.cores = 2
                       )
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n df2b: the seed is", pseed, "\n")
    b2 <- oncoSimulPop(reps,
                       fe1,
                       mu = m2,
                       onlyCancer = FALSE,
                       initSize = no,
                       finalTime = ft,
                       mutationPropGrowth = FALSE, ## cleaner, though no real effect
                       seed =NULL,
                       mc.cores = 2
                       )
    expect_true( median(summary(b2)$NumClones) >
                 median(summary(b1)$NumClones))
    ## Note the short time, so this is not always very different as few
    ## have double or triple mutants
    expect_true( mean(mutsPerClone(b2)) >
                 mean(mutsPerClone(b1)))
})
date()








date()
test_that("McFL: Expect freqs, num clones, muts per clone for different per-gene-mut",{
## More on the above, with less variation. But yet another set of tests.
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n mcdf1: the seed is", pseed, "\n")
    ng <- 10
    ni <- rep(0, ng)
    m1 <- runif(ng, min = 1e-6, max = 1e-5)
    m2 <- runif(ng, min = 1e-4, max = 1e-3)
    names(ni) <- names(m1) <- names(m2) <- c(replicate(ng,
                                 paste(sample(letters, 12), collapse = "")))
    fe1 <- allFitnessEffects(noIntGenes = ni)
    ft <- 0.01
    no <- 1e7
    reps <- 80
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n mcdf1a: the seed is", pseed, "\n")
    b1 <- oncoSimulPop(reps,
                       fe1,
                       mu = m1,
                       onlyCancer = FALSE,
                       initSize = no,
                       finalTime = ft,
                       mutationPropGrowth = FALSE, ## cleaner, though no real effect
                       seed =NULL,
                       model = "McFL",
                       mc.cores = 2
                       )
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n mcdf1b: the seed is", pseed, "\n")
    b2 <- oncoSimulPop(reps,
                       fe1,
                       mu = m2,
                       onlyCancer = FALSE,
                       initSize = no,
                       finalTime = ft,
                       mutationPropGrowth = FALSE, ## cleaner, though no real effect
                       seed =NULL,
                       model = "McFL",
                       mc.cores = 2
                       )
    (expected1 <- no*reps*m1)
    (expected2 <- no*reps*m2)
    (cc1 <- colSums(OncoSimulR:::geneCounts(b1)))
    (cc2 <- colSums(OncoSimulR:::geneCounts(b2)))    
    ## It will fail with prob ~ p.fail
    p.fail <- 1e-3
    expect_true(chisq.test(colSums(OncoSimulR:::geneCounts(b1)),
                           p = expected1/sum(expected1))$p.value > p.fail)
    expect_true(chisq.test(colSums(OncoSimulR:::geneCounts(b2)),
                           p = expected2/sum(expected2))$p.value > p.fail)
    expect_true( median(summary(b2)$NumClones) >
                 median(summary(b1)$NumClones))
    ## Note the short time, so this is not always very different as few
    ## have double or triple mutants
    expect_true( mean(mutsPerClone(b2)) >
                 mean(mutsPerClone(b1)))
})

date()
test_that("MCFL: Num clones, muts per clone for different per-gene-mut",{
    ## Like previous, but larger finalTime, so no longer chi-square test
    ## here.
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n mcdf2: the seed is", pseed, "\n")
    ng <- 40
    ni <- rep(0, ng)
    m1 <- runif(ng, min = 1e-6, max = 1e-5)
    m2 <- runif(ng, min = 1e-4, max = 1e-3)
    names(ni) <- names(m1) <- names(m2) <- c(replicate(ng,
                                 paste(sample(letters, 12), collapse = "")))
    fe1 <- allFitnessEffects(noIntGenes = ni)
    ft <- 2
    no <- 1e5
    reps <- 20
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n mcdf2a: the seed is", pseed, "\n")
    b1 <- oncoSimulPop(reps,
                       fe1,
                       mu = m1,
                       onlyCancer = FALSE,
                       initSize = no,
                       finalTime = ft,
                       mutationPropGrowth = FALSE, ## cleaner, though no real effect
                       seed =NULL,
                       model = "McFL",
                       mc.cores = 2
                       )
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n mcdf2b: the seed is", pseed, "\n")
    b2 <- oncoSimulPop(reps,
                       fe1,
                       mu = m2,
                       onlyCancer = FALSE,
                       initSize = no,
                       finalTime = ft,
                       mutationPropGrowth = FALSE, ## cleaner, though no real effect
                       seed =NULL,
                       model = "McFL",
                       mc.cores = 2
                       )
    expect_true( median(summary(b2)$NumClones) >
                 median(summary(b1)$NumClones))
    ## Note the short time, so this is not always very different as few
    ## have double or triple mutants
    expect_true( mean(mutsPerClone(b2)) >
                 mean(mutsPerClone(b1)))
})
date()


## Most tests above with t >> 0.01 or so have mutationPropGrowth =
## FALSE. Why? mutationPropGrowth will not have any noticeable effect
## unless we let it run for some time and unless there are sizeable
## differences in birth rates between clones. So in most cases above
## setting it to FALSE makes little difference, but just to be cleaner.






## Much nicer: we stop on population size. With mutPropGrowth = TRUE, we
## actually stop slightly earlier. But we have much larger numbers of
## clones, etc. So we are not affected by issues of differences in
## populationSize.

date()
test_that(" And mutPropGrowth, 3",{
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n sz033: the seed is", pseed, "\n")
    muvar2 <- c("U" = 1e-4, "z" = 5e-5, "e" = 5e-4, "m" = 5e-3, "D" = 1e-4)
    ## muvar2 <- c("U" = 5e-5, "z" = 5e-5, "e" = 5e-5, "m" = 5e-5, "D" = 5e-5)
    ni1 <- rep(1.9, 5)
    names(ni1) <- names(muvar2)
    fe1 <- allFitnessEffects(noIntGenes = ni1)
    no <- 1e5
    reps <- 10 
    ft <- 16
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n sz033a: the seed is", pseed, "\n")
    b1 <- oncoSimulPop(reps,
                       fe1, mu = muvar2,
                       mutationPropGrowth = FALSE,
                       onlyCancer = FALSE,
                       initSize = no,
                       finalTime = ft,
                       detectionSize = 5e6,
                       sampleEvery = 0.01,
                       seed =NULL,
                       mc.cores = 2
                       )
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n sz033b: the seed is", pseed, "\n")
    b2 <- oncoSimulPop(reps,
                       fe1, mu = muvar2,
                       mutationPropGrowth = TRUE,
                       onlyCancer = FALSE,
                       initSize = no,
                       finalTime = ft,
                       detectionSize = 5e6,
                       sampleEvery = 0.01,                       
                       seed =NULL,
                       mc.cores = 2
                       )
    summary(b2)[, c(1:3, 8:9)]
    mean(mutsPerClone(b1));mean(mutsPerClone(b2))
    median(summary(b1)$NumClones)
    median(summary(b2)$NumClones)
    ## More mutations in mutationPropGrowth
    expect_true( mean(mutsPerClone(b2)) >
                 mean(mutsPerClone(b1)))
    expect_true( median(summary(b2)$NumClones) >
                 median(summary(b1)$NumClones))
    ## But frequency of mutations about the same? Nope: since very fast
    ## growth and thus non-indep, huge variation in geneCounts in each
    ## run, etc. so hard to compare geneCounts
    ## Just for reference, here
    ## First, look at run to run variation
    ## OncoSimulR:::geneCounts(b1)
    ## OncoSimulR:::geneCounts(b2)
    ## ## The next makes sense
    ## fb1 <- colSums(OncoSimulR:::geneCounts(b1))
    ## fb2 <- colSums(OncoSimulR:::geneCounts(b2))
    ## fb1
    ## fb2
    ## fb2/fb1
    ## fb1/sum(fb1)
    ## fb2/sum(fb2)
    ## (fb2/sum(fb2))/(fb1/sum(fb1))
})
date()




date()
test_that(" McFL: And mutPropGrowth, 3",{
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n mcsz033: the seed is", pseed, "\n")
    muvar2 <- c("U" = 1e-4, "z" = 5e-5, "e" = 5e-4, "m" = 5e-3, "D" = 1e-4)
    ## muvar2 <- c("U" = 5e-5, "z" = 5e-5, "e" = 5e-5, "m" = 5e-5, "D" = 5e-5)
    ni1 <- rep(1.9, 5)
    names(ni1) <- names(muvar2)
    fe1 <- allFitnessEffects(noIntGenes = ni1)
    no <- 1e5
    reps <- 10 
    ft <- 16
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n mcsz033a: the seed is", pseed, "\n")
    b1 <- oncoSimulPop(reps,
                       fe1, mu = muvar2,
                       mutationPropGrowth = FALSE,
                       onlyCancer = FALSE,
                       initSize = no,
                       finalTime = ft,
                       detectionSize = 5e6,
                       sampleEvery = 0.01,
                       model = "McFL",
                       seed =NULL,
                       mc.cores = 2
                       )
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n mcsz033b: the seed is", pseed, "\n")
    b2 <- oncoSimulPop(reps,
                       fe1, mu = muvar2,
                       mutationPropGrowth = TRUE,
                       onlyCancer = FALSE,
                       initSize = no,
                       finalTime = ft,
                       detectionSize = 5e6,
                       sampleEvery = 0.01,
                       model = "McFL",
                       seed =NULL,
                       mc.cores = 2
                       )
    summary(b1)[, c(1:3, 8:9)]
    summary(b2)[, c(1:3, 8:9)]
    mean(mutsPerClone(b1));mean(mutsPerClone(b2))
    median(summary(b1)$NumClones)
    median(summary(b2)$NumClones)
    ## More mutations in mutationPropGrowth
    expect_true( mean(mutsPerClone(b2)) >
                 mean(mutsPerClone(b1)))
    expect_true( median(summary(b2)$NumClones) >
                 median(summary(b1)$NumClones))
    ## But frequency of mutations about the same? Nope: since very fast
    ## growth and thus non-indep, huge variation in geneCounts in each
    ## run, etc. so hard to compare geneCounts
    ## Just for reference, here
    ## First, look at run to run variation
    ## OncoSimulR:::geneCounts(b1)
    ## OncoSimulR:::geneCounts(b2)
    ## ## The next makes sense
    ## fb1 <- colSums(OncoSimulR:::geneCounts(b1))
    ## fb2 <- colSums(OncoSimulR:::geneCounts(b2))
    ## fb1
    ## fb2
    ## fb2/fb1
    ## fb1/sum(fb1)
    ## fb2/sum(fb2)
    ## (fb2/sum(fb2))/(fb1/sum(fb1))
})
date()









date()
test_that("More mutpropgrowth, in modules of s", {
    ## From a similar test in mutPropGrowth, but we have a vector mu
    ## And here, we fix detectionSize, so effects are not due
    ## to larger population sizes.

    ## As previously, stop on population size
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n mpgs3: the seed is", pseed, "\n")
    pops <- 20
    lni <- 1 ## no fitness effects genes
    fni <- 50 ## fitness effects genes
    no <- 1e3
    ft <- 10 ## 5
    s3 <- 3.0
    ## noInt have no fitness effects, but can accumulate mutations
    ni <- rep(0, lni)
    names(ni) <- paste0("ni", 1:lni)
    ## Those with fitness effects in one module, so
    ## neither fitness nor mut. rate blow up
    gn <- paste(paste0("a", 1:fni), collapse = ", ")
    f3 <- allFitnessEffects(epistasis = c("A" = s3),
                            geneToModule = c("A" = gn),
                            noIntGenes = ni)
    mu <- runif(fni + lni, min = 1e-7, max = 1e-4)
    names(mu) <- c(paste0("a", 1:fni), names(ni))
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n mpgs3a: the seed is", pseed, "\n")
    s3.ng <- oncoSimulPop(pops,
                          f3,
                          mu = mu,
                          mutationPropGrowth = FALSE,
                          finalTime =ft,
                          sampleEvery = 0.01,
                          detectionSize = 1e6,
                          detectionDrivers = 9999,
                          initSize = no,
                          onlyCancer = FALSE,
                          seed = NULL, mc.cores = 2)
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n mpgs3b: the seed is", pseed, "\n")
    s3.g <- oncoSimulPop(pops,
                         f3,
                         mu = mu,
                         mutationPropGrowth = TRUE,
                         finalTime =ft,
                         sampleEvery = 0.01,
                         detectionSize = 1e6,
                         detectionDrivers = 9999,
                         initSize = no,
                         onlyCancer = FALSE,
                         seed = NULL, mc.cores = 2)
    summary(s3.g)[, c(1, 2, 3, 8, 9)]
    summary(s3.ng)[, c(1, 2, 3, 8, 9)]
    summary(summary(s3.ng)[, 2])
    summary(summary(s3.g)[, 2])
    expect_true( mean(mutsPerClone(s3.g)) >
                 mean(mutsPerClone(s3.ng)))
    expect_true( median(summary(s3.g)$NumClones) >
                 median(summary(s3.ng)$NumClones))
})
date()





date()
test_that("McFL: More mutpropgrowth, in modules of s", {
    ## From a similar test in mutPropGrowth, but we have a vector mu

    ## With McFL, since total size is bounded from above, the effects of
    ## mutPropGrowth are rarely those of popSize. But we stop on popSize
    ## anyway here.
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n mcmpgs3: the seed is", pseed, "\n")
    pops <- 20
    lni <- 1 ## no fitness effects genes
    fni <- 50 ## fitness effects genes
    no <- 1e3
    ft <-  15 ## 5
    s3 <- 3.0
    ## noInt have no fitness effects, but can accumulate mutations
    ni <- rep(0, lni)
    names(ni) <- paste0("ni", 1:lni)
    ## Those with fitness effects in one module, so
    ## neither fitness nor mut. rate blow up
    gn <- paste(paste0("a", 1:fni), collapse = ", ")
    f3 <- allFitnessEffects(epistasis = c("A" = s3),
                            geneToModule = c("A" = gn),
                            noIntGenes = ni)
    mu <- runif(fni + lni, min = 1e-7, max = 1e-4)
    names(mu) <- c(paste0("a", 1:fni), names(ni))
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n mcmpgs3a: the seed is", pseed, "\n")
    s3.ng <- oncoSimulPop(pops,
                          f3,
                          mu = mu,
                          mutationPropGrowth = FALSE,
                          finalTime =ft,
                          sampleEvery = 0.01,
                          detectionSize = 2.5e4,
                          detectionDrivers = 9999,
                          initSize = no,
                          onlyCancer = FALSE,
                          model = "McFL",
                          seed = NULL, mc.cores = 2)
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n mcmpgs3b: the seed is", pseed, "\n")
    s3.g <- oncoSimulPop(pops,
                         f3,
                         mu = mu,
                         mutationPropGrowth = TRUE,
                         finalTime =ft,
                         sampleEvery = 0.01,
                         detectionSize = 2.5e4,
                         detectionDrivers = 9999,
                         initSize = no,
                         onlyCancer = FALSE,
                         model = "McFL",
                         seed = NULL, mc.cores = 2)
    summary(s3.g)[, c(1, 2, 3, 8, 9)]
    summary(s3.ng)[, c(1, 2, 3, 8, 9)]
    summary(summary(s3.ng)[, 2])
    summary(summary(s3.g)[, 2])
    expect_true( mean(mutsPerClone(s3.g)) >
                 mean(mutsPerClone(s3.ng)))
    expect_true( median(summary(s3.g)$NumClones) >
                 median(summary(s3.ng)$NumClones))
})
date()




date()
test_that("oncoSimulSample: expected vs. observed for different per-gene-mut",{
    ## Here, we test that freqs as they should, but so that the test is
    ## not eternal, we use different settings of reps and no

    ## We want to get about 2 mutations in each population, but not more
    ## than 2, as that would mean we are not picking diffrences between
    ## muts with different mut rate, because we have
    ## wholePopulationSample.
    
    ## We probably want about a mean or median number of clones of about 2
    ## or so. Though if fewer, better but then to have power in the
    ## chi-square we need much larger reps (as usual, if ft increase, etc,
    ## we increase the reproduction/death events, which then screws up
    ## simple expectations for chi-square)
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n oss11: the seed is", pseed, "\n")
    ng <- 10
    ni <- rep(0, ng)
    m1 <- seq(from = 1e-7, to = 1e-4, length.out = ng)
    m2 <- runif(ng, min = 1e-6, max = 1e-3)
    names(ni) <- names(m1) <- names(m2) <- c(replicate(ng,
                                 paste(sample(letters, 12), collapse = "")))
    fe1 <- allFitnessEffects(noIntGenes = ni)
    ft <- 5e-3
    no <- 5e5 # delicate as if this is huge, we get the cc1 or cc2 below
              # to be equal to reps in many genes, because they are
              # present in at least one cell in all populations
    reps <- 5000 ## large because few events with small mut
                 ## freqs. o.w. chi has many small cells.
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    x <- 1e-20
    cat("\n oss1a: the seed is", pseed, "\n")
    b1 <- oncoSimulSample(reps,
                          fe1,
                          mu = m1,
                          onlyCancer = FALSE,
                          initSize = no,
                          finalTime = ft,
                          mutationPropGrowth = FALSE, ## cleaner, though no real effect
                          seed =NULL,
                          sampleEvery = 0.01, thresholdWhole = 1e-20,
                          detectionSize = 1e9,
                          detectionDrivers = 9999,
                          )
    b1$popSummary[, c(1:3, 8:9)]
    summary(b1$popSummary[, "NumClones"])
    (expected1 <- no*reps*m1)
    (cc1 <- colSums(b1$popSample))
    if( (any(cc1 == reps)) )
        warning("The test is likely to fail because reps == cc1 or cc2")
    ## It will fail with prob ~ p.fail
    p.fail <- 1e-3
    expect_true(chisq.test(cc1,
                           p = expected1/sum(expected1))$p.value > p.fail)

    reps <- 500
    no <- 1e4
    ft <- 0.03
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n oss1b: the seed is", pseed, "\n")
    b2 <- oncoSimulSample(reps,
                       fe1,
                       mu = m2,
                       onlyCancer = FALSE,
                       initSize = no,
                       finalTime = ft,
                       mutationPropGrowth = FALSE, ## cleaner, though no real effect
                       seed =NULL,
                       sampleEvery = 0.01, thresholdWhole = 1e-20,
                       detectionSize = 1e9,
                       detectionDrivers = 9999,
                       )
    summary(b2$popSummary[, "NumClones"])
    ## we detect anything that is present in at least one case.
    ## Not exactly the same as what we did in oncoSimulPop
    (expected2 <- no*reps*m2)
    (cc2 <- colSums(b2$popSample))
    if( (any(cc2 == reps)))
        warning("The test is likely to fail because reps == cc1 or cc2")
    ## It will fail with prob ~ p.fail
    p.fail <- 1e-3
    expect_true(chisq.test(cc2,
                           p = expected2/sum(expected2))$p.value > p.fail)


})
date()





date()
test_that("oncoSimulSample comparing different per-gene-mut",{
    ## No attempt to compare against expected (other tests do that). We
    ## just verify that larger mutations rates lead to more total
    ## mutations and clones.
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n oss2: the seed is", pseed, "\n")
    ng <- 10
    ni <- rep(0, ng)
    m1 <- runif(ng, min = 1e-6, max = 1e-5)
    m2 <- runif(ng, min = 1e-4, max = 1e-3)
    names(ni) <- names(m1) <- names(m2) <- c(replicate(ng,
                                 paste(sample(letters, 12), collapse = "")))
    fe1 <- allFitnessEffects(noIntGenes = ni)
    ft <- .05 ## if you make it too large, then all pops will have at
             ## least one cell with one of the genes mutated. You can see
             ## this when cc1 or cc2 have most/all entries equal to reps.
    no <- 1e5 # delicate as if this is huge, we get the cc1 or cc2 below
              # to be equal to reps in many genes, because they are
              # present in at least one cell in all populations
    reps <- 500
    x <- 1e-20
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n oss2a: the seed is", pseed, "\n")
    b1 <- oncoSimulSample(reps,
                          fe1,
                          mu = m1,
                          onlyCancer = FALSE,
                          initSize = no,
                          finalTime = ft,
                          mutationPropGrowth = FALSE, ## cleaner, though no real effect
                          seed =NULL,
                          thresholdWhole = x
                          )
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n oss2b: the seed is", pseed, "\n")
    b2 <- oncoSimulSample(reps,
                       fe1,
                       mu = m2,
                       onlyCancer = FALSE,
                       initSize = no,
                       finalTime = ft,
                       mutationPropGrowth = FALSE, ## cleaner, though no real effect
                       seed =NULL,
                       thresholdWhole = x
                       )
    ## we detect anything that is present in at least one case.
    ## Not exactly the same as what we did in oncoSimulPop
    (cc1 <- colSums(b1$popSample))
    (cc2 <- colSums(b2$popSample))
    expect_true(sum(cc2) > sum(cc1))
    ## This is very similar to above, like assimilating a pop to a clone
    mutsPerClone1 <- rowSums(b1$popSample)
    mutsPerClone2 <- rowSums(b2$popSample)
    expect_true( mean(mutsPerClone2) >
                 mean(mutsPerClone1))
    expect_true( median(b2$popSummary[, "NumClones"]) >
                 median(b1$popSummary[, "NumClones"]))
    ## Note the short time, so this is not always very different as few
    ## have double or triple mutants
})
date()

date()
test_that("McFL: oncoSimulSample: expected vs. observed for different per-gene-mut",{
    ## Here, we test that freqs as they should, but so that the test is
    ## not eternal, we use different settings of reps and no

    ## We probably want about a mean or median number of clones of about 2
    ## or so. Though if fewer, better but then to have power in the
    ## chi-square we need much larger reps (as usual, if ft increase, etc,
    ## we increase the reproduction/death events, which then screws up
    ## simple expectations for chi-square)
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n mcoss11: the seed is", pseed, "\n")
    ng <- 10
    ni <- rep(0, ng)
    m1 <- runif(ng, min = 5e-8, max = 1e-5)
    m2 <- runif(ng, min = 1e-4, max = 1e-3)
    names(ni) <- names(m1) <- names(m2) <- c(replicate(ng,
                                 paste(sample(letters, 12), collapse = "")))
    fe1 <- allFitnessEffects(noIntGenes = ni)
    ft <- 0.03
    no <- 5e5 # delicate as if this is huge, we get the cc1 or cc2 below
              # to be equal to reps in many genes, because they are
              # present in at least one cell in all populations
    reps <- 600
    x <- 1e-20
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n mcoss1a: the seed is", pseed, "\n")
    b1 <- oncoSimulSample(reps,
                          fe1,
                          mu = m1,
                          onlyCancer = FALSE,
                          initSize = no,
                          finalTime = ft,
                          mutationPropGrowth = FALSE, ## cleaner, though no real effect
                          seed =NULL,
                          model = "McFL",
                          thresholdWhole = x
                          )
    summary(b1$popSummary[, "NumClones"])
    b1$popSummary[, c(1:3, 8:9)]
    (expected1 <- no*reps*m1)
    (cc1 <- colSums(b1$popSample))
    if( (any(cc1 == reps)) )
        warning("The test is likely to fail because reps == cc1 or cc2")
    ## It will fail with prob ~ p.fail
    p.fail <- 1e-3
    expect_true(chisq.test(cc1,
                           p = expected1/sum(expected1))$p.value > p.fail)
    reps <- 500
    no <- 1e4
    ft <- 0.03
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n mcoss1b: the seed is", pseed, "\n")
    b2 <- oncoSimulSample(reps,
                       fe1,
                       mu = m2,
                       onlyCancer = FALSE,
                       initSize = no,
                       finalTime = ft,
                       mutationPropGrowth = FALSE, ## cleaner, though no real effect
                       seed =NULL,
                       model = "McFL",
                       thresholdWhole = x
                       )
    summary(b2$popSummary[, "NumClones"])
    ## we detect anything that is present in at least one case.
    ## Not exactly the same as what we did in oncoSimulPop
    (expected2 <- no*reps*m2)
    (cc2 <- colSums(b2$popSample))
    if( (any(cc2 == reps)))
        warning("The test is likely to fail because reps == cc1 or cc2")
    ## It will fail with prob ~ p.fail
    p.fail <- 1e-3
    expect_true(chisq.test(cc2,
                           p = expected2/sum(expected2))$p.value > p.fail)

    
})
date()





date()
test_that("McFL: oncoSimulSample comparing different per-gene-mut",{
    ## No attempt to compare against expected (other tests do that). We
    ## just verify that larger mutations rates lead to more total
    ## mutations and clones.
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n mcoss2: the seed is", pseed, "\n")
    ng <- 10
    ni <- rep(0, ng)
    m1 <- runif(ng, min = 1e-6, max = 1e-5)
    m2 <- runif(ng, min = 1e-4, max = 1e-3)
    names(ni) <- names(m1) <- names(m2) <- c(replicate(ng,
                                 paste(sample(letters, 12), collapse = "")))
    fe1 <- allFitnessEffects(noIntGenes = ni)
    ft <- .05 ## if you make it too large, then all pops will have at
             ## least one cell with one of the genes mutated. You can see
             ## this when cc1 or cc2 have most/all entries equal to reps.
    no <- 1e5 # delicate as if this is huge, we get the cc1 or cc2 below
              # to be equal to reps in many genes, because they are
              # present in at least one cell in all populations
    reps <- 500
    x <- 1e-20
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n mcoss2a: the seed is", pseed, "\n")
    b1 <- oncoSimulSample(reps,
                          fe1,
                          mu = m1,
                          onlyCancer = FALSE,
                          initSize = no,
                          finalTime = ft,
                          mutationPropGrowth = FALSE, ## cleaner, though no real effect
                          seed =NULL,
                          model = "McFL",
                          thresholdWhole = x
                          )
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n mcoss2b: the seed is", pseed, "\n")
    b2 <- oncoSimulSample(reps,
                       fe1,
                       mu = m2,
                       onlyCancer = FALSE,
                       initSize = no,
                       finalTime = ft,
                       mutationPropGrowth = FALSE, ## cleaner, though no real effect
                       seed =NULL,
                       model = "McFL",
                       thresholdWhole = x
                       )
    ## we detect anything that is present in at least one case.
    ## Not exactly the same as what we did in oncoSimulPop
    (cc1 <- colSums(b1$popSample))
    (cc2 <- colSums(b2$popSample))
    expect_true(sum(cc2) > sum(cc1))
    ## This is very similar to above, like assimilating a pop to a clone
    mutsPerClone1 <- rowSums(b1$popSample)
    mutsPerClone2 <- rowSums(b2$popSample)
    expect_true( mean(mutsPerClone2) >
                 mean(mutsPerClone1))
    expect_true( median(b2$popSummary[, "NumClones"]) >
                 median(b1$popSummary[, "NumClones"]))
    ## Note the short time, so this is not always very different as few
    ## have double or triple mutants
})
date()





## From similar tests in mutPropGrwoth-long, but here mu is a vector

cat("\n", date(), "\n")
test_that("oncoSimulSample Without initmutant and modules, fixed size", {
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n osSFPS: the seed is", pseed, "\n")
    pops <- 60
    lni <- 1 ## no fitness effects genes
    fni <- 50 ## fitness effects genes
    no <- 1e4 
    ft <- 9  #4 
    s3 <- 2.5 
    mu <- 1e-5 
    ## noInt have no fitness effects, but can accumulate mutations
    ni <- rep(0, lni)
    names(ni) <- paste0("ni", 1:lni)
    ## Those with fitness effects in one module, so
    ## neither fitness nor mut. rate blow up
    gn <- paste0("a", 1:fni)
    mu <- runif(lni + fni, min = 1e-7, max = 1e-4)
    names(mu) <- c(gn, names(ni))
    gn <- paste(gn, collapse = ", ")
    f3 <- allFitnessEffects(epistasis = c("A" = s3),
                            geneToModule = c("A" = gn),
                            noIntGenes = ni)
   
    x <- 1e-9 ## so basically anything that appears once
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n osSFPSa: the seed is", pseed, "\n")
    b1 <- oncoSimulSample(pops,
                          f3,
                          mu = mu,
                          mutationPropGrowth = FALSE,
                          finalTime =ft,
                          initSize = no,
                          onlyCancer = FALSE,
                          sampleEvery = 0.01,
                          detectionSize = 6e4,
                          detectionDrivers = 99,
                          seed =NULL,
                          thresholdWhole = x)
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n osSFPSb: the seed is", pseed, "\n")
    b2 <- oncoSimulSample(pops,
                         f3,
                         mu = mu,
                         mutationPropGrowth = TRUE,
                         finalTime =ft,
                         initSize = no,
                         onlyCancer = FALSE,
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
    expect_true( mean(mutsPerClone2) >
                 mean(mutsPerClone1))
    expect_true( median(b2$popSummary[, "NumClones"]) >
                 median(b1$popSummary[, "NumClones"]))
})
cat("\n", date(), "\n")


cat("\n", date(), "\n")
test_that("oncoSimulSample Without initmutant and modules, McFL, fixed size", {
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n osSFPSMcFL: the seed is", pseed, "\n")
    pops <- 60
    lni <- 1 ## no fitness effects genes
    fni <- 50 ## fitness effects genes
    no <- 1e4 ## note we use only 10 in the other example below
    ft <- 10  ##4 
    s3 <- 2.5 
    mu <- 1e-5 
    ## noInt have no fitness effects, but can accumulate mutations
    ni <- rep(0, lni)
    names(ni) <- paste0("ni", 1:lni)
    ## Those with fitness effects in one module, so
    ## neither fitness nor mut. rate blow up
    gn <- paste0("a", 1:fni)
    mu <- runif(lni + fni, min = 1e-7, max = 1e-4)
    names(mu) <- c(gn, names(ni))
    gn <- paste(gn, collapse = ", ")
    ## Those with fitness effects in one module, so
    ## neither fitness nor mut. rate blow up
    gn <- paste(paste0("a", 1:fni), collapse = ", ")
    f3 <- allFitnessEffects(epistasis = c("A" = s3),
                            geneToModule = c("A" = gn),
                            noIntGenes = ni)
    x <- 1e-9 ## 1/no
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n osSFPSMcFLa: the seed is", pseed, "\n")
    b1 <- oncoSimulSample(pops,
                          f3,
                          mu = mu,
                          mutationPropGrowth = FALSE,
                          finalTime =ft,
                          initSize = no,
                          onlyCancer = FALSE,
                          sampleEvery = 0.01,
                          detectionSize = 2.5e4,
                          detectionDrivers = 99,
                          seed =NULL,
                          model = "McFL",
                          thresholdWhole = x)
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n osSFPSMcFLb: the seed is", pseed, "\n")
    b2 <- oncoSimulSample(pops,
                         f3,
                         mu = mu,
                         mutationPropGrowth = TRUE,
                         finalTime =ft,
                         initSize = no,
                         onlyCancer = FALSE,
                         sampleEvery = 0.01,
                          detectionSize = 2.5e4,
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



## Repeat some above, allowing for mutPropGrowth

date()
test_that("Different freqs as they should be ordered and chisq, when s  and a tiny mu", {
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n mpg s13: the seed is", pseed, "\n")
    muvar2 <- c("U" = 1e-13, "z" = 1e-7, "e" = 1e-6, "m" = 1e-5, "D" = 1e-4)
    ni2 <- rep(0.01, 5)
    names(ni2) <- names(muvar2)
    ni2["U"] <- 0.5
    fe1 <- allFitnessEffects(noIntGenes = ni2)
    no <- 1e6
    reps <- 600
    bb <- oncoSimulPop(reps,
                       fe1, mu = muvar2, onlyCancer = FALSE,
                       initSize = no,
                       finalTime = 5,
                       mutationPropGrowth = TRUE,
                       seed = NULL, mc.cores = 2
                       )
    (expectedC <- no*reps*muvar2)
    colSums(OncoSimulR:::geneCounts(bb))
    expect_true(colSums(OncoSimulR:::geneCounts(bb))[1] == 0)
    expect_equal(
        order(colSums(OncoSimulR:::geneCounts(bb))),
        order(expectedC))
    ## A chisq will not work as we increase finalTime.
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n mpg s13b: the seed is", pseed, "\n")
    bb <- oncoSimulPop(reps,
                       fe1, mu = muvar2, onlyCancer = FALSE,
                       initSize = no,
                       finalTime = .001,
                       mutationPropGrowth = TRUE,
                       seed = NULL, mc.cores = 2
                       )
    (expectedC <- no*reps*muvar2)
    colSums(OncoSimulR:::geneCounts(bb))
    expect_true(colSums(OncoSimulR:::geneCounts(bb))[1] == 0)
    ## This will fail sometimes
    p.fail <- 1e-3
    expect_true(chisq.test(colSums(OncoSimulR:::geneCounts(bb))[-1],
                           p = expectedC[-1]/sum(expectedC))$p.value > p.fail)
    ## expect_equal(
    ##     order(colSums(OncoSimulR:::geneCounts(bb))),
    ##     order(expectedC))
}
)
date()

date()
test_that("McFL: Different freqs as they should be ordered and chisq, when s  and a tiny mu", {
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n mpg s14: the seed is", pseed, "\n")
    muvar2 <- c("U" = 1e-13, "z" = 5e-7, "e" = 1e-6, "m" = 1e-5, "D" = 1e-4)
    ni2 <- rep(0.01, 5)
    names(ni2) <- names(muvar2)
    ni2["U"] <- 0.5
    fe1 <- allFitnessEffects(noIntGenes = ni2)
    no <- 1e5
    reps <- 400
    ## Beware: with McFL final time has to be small or we will see
    ## mutations in U as it is the only one left. 
    bb <- oncoSimulPop(reps,
                       fe1, mu = muvar2, onlyCancer = FALSE,
                       initSize = no,
                       model = "McFL",
                       finalTime = 50,
                       mutationPropGrowth = TRUE,
                       seed = NULL, mc.cores = 2
                       )
    (expectedC <- no*reps*muvar2)
    colSums(OncoSimulR:::geneCounts(bb))
    expect_true(colSums(OncoSimulR:::geneCounts(bb))[1] == 0)
    expect_equal(
        order(colSums(OncoSimulR:::geneCounts(bb))),
        order(expectedC))
    ## A chisq will not work as we increase finalTime.
    no <- 1e7
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n mpg s14b: the seed is", pseed, "\n")
    bb <- oncoSimulPop(reps,
                       fe1, mu = muvar2, onlyCancer = FALSE,
                       initSize = no,
                       model = "McFL",
                       finalTime = .001,
                       mutationPropGrowth = TRUE,
                       seed = NULL, mc.cores = 2
                       )
    (expectedC <- no*reps*muvar2)
    colSums(OncoSimulR:::geneCounts(bb))
    expect_true(colSums(OncoSimulR:::geneCounts(bb))[1] == 0)
    ## This will fail sometimes
    p.fail <- 1e-3
    expect_true(chisq.test(colSums(OncoSimulR:::geneCounts(bb))[-1],
                           p = expectedC[-1]/sum(expectedC))$p.value > p.fail)
    ## expect_equal(
    ##     order(colSums(OncoSimulR:::geneCounts(bb))),
    ##     order(expectedC))
}
)
date()





date()
test_that("Num clones, muts per clone for different per-gene-mut",{
    ## Like previous, but larger finalTime, so no longer chi-square test
    ## here.
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n mpg df2: the seed is", pseed, "\n")
    ng <- 40
    ni <- rep(0, ng)
    m1 <- runif(ng, min = 1e-6, max = 1e-5)
    m2 <- runif(ng, min = 1e-4, max = 1e-3)
    names(ni) <- names(m1) <- names(m2) <- c(replicate(ng,
                                 paste(sample(letters, 12), collapse = "")))
    fe1 <- allFitnessEffects(noIntGenes = ni)
    ft <- 2
    no <- 1e5
    reps <- 20
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n mpg df2a: the seed is", pseed, "\n")
    b1 <- oncoSimulPop(reps,
                       fe1,
                       mu = m1,
                       onlyCancer = FALSE,
                       initSize = no,
                       finalTime = ft,
                       mutationPropGrowth = TRUE, 
                       seed =NULL,
                       mc.cores = 2
                       )
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n mpg df2b: the seed is", pseed, "\n")
    b2 <- oncoSimulPop(reps,
                       fe1,
                       mu = m2,
                       onlyCancer = FALSE,
                       initSize = no,
                       finalTime = ft,
                       mutationPropGrowth = TRUE, 
                       seed =NULL,
                       mc.cores = 2
                       )
    expect_true( median(summary(b2)$NumClones) >
                 median(summary(b1)$NumClones))
    ## Note the short time, so this is not always very different as few
    ## have double or triple mutants
    expect_true( mean(mutsPerClone(b2)) >
                 mean(mutsPerClone(b1)))
})
date()

date()
test_that("MCFL: Num clones, muts per clone for different per-gene-mut",{
    ## Like previous, but larger finalTime, so no longer chi-square test
    ## here.
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n mpg mcdf2: the seed is", pseed, "\n")
    ng <- 40
    ni <- rep(0, ng)
    m1 <- runif(ng, min = 1e-6, max = 1e-5)
    m2 <- runif(ng, min = 1e-4, max = 1e-3)
    names(ni) <- names(m1) <- names(m2) <- c(replicate(ng,
                                 paste(sample(letters, 12), collapse = "")))
    fe1 <- allFitnessEffects(noIntGenes = ni)
    ft <- 2
    no <- 1e5
    reps <- 20
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n mpg mcdf2a: the seed is", pseed, "\n")
    b1 <- oncoSimulPop(reps,
                       fe1,
                       mu = m1,
                       onlyCancer = FALSE,
                       initSize = no,
                       finalTime = ft,
                       mutationPropGrowth = TRUE,
                       seed =NULL,
                       model = "McFL",
                       mc.cores = 2
                       )
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n mpg mcdf2b: the seed is", pseed, "\n")
    b2 <- oncoSimulPop(reps,
                       fe1,
                       mu = m2,
                       onlyCancer = FALSE,
                       initSize = no,
                       finalTime = ft,
                       mutationPropGrowth = TRUE, 
                       seed =NULL,
                       model = "McFL",
                       mc.cores = 2
                       )
    expect_true( median(summary(b2)$NumClones) >
                 median(summary(b1)$NumClones))
    ## Note the short time, so this is not always very different as few
    ## have double or triple mutants
    expect_true( mean(mutsPerClone(b2)) >
                 mean(mutsPerClone(b1)))
})
date()



##################### If you want to verify step by step that the C++ does
##################### what it should you can, for instance, run this R code
## library(OncoSimulR)
## RNGkind("L'Ecuyer-CMRG")
## set.seed(13)
## muvar2 <- c("U" = 1e-3, "z" = 3e-7, "e" = 5e-6, "m" = 5e-5, "D" = 5e-4)
## ni2 <- rep(0.01, 5)
## names(ni2) <- names(muvar2)
## fe1 <- allFitnessEffects(noIntGenes = ni2)
## no <- 1e5
## bb <- oncoSimulIndiv(fe1, mu = muvar2, onlyCancer = FALSE,
##                      mutationPropGrowth = FALSE,
##                      initSize = no,
##                      finalTime = 560
##                    )

## bb
## with the following C++ in BNB_nr.cpp, right after the line
## tmpParam.mutation = mutationFromParent(mu, tmpParam, popParams[nextMutant],
## 				   newMutations, mutationPropGrowth);
## Add this C++ code, recompile.
## DP1("at mutation");
## 	    Rcpp::Rcout << "\n New Genotype :";
## 	    print_Genotype(newGenotype);
## 	    Rcpp::Rcout << "\n Parent Genotype :";
## 	    print_Genotype(Genotypes[nextMutant]);
## 	    DP2(tmpParam.mutation);
## 	    DP2( popParams[nextMutant].mutation);
## 	    DP2(mutationFromScratch(mu, tmpParam, newGenotype,
## 					       fitnessEffects,
## 				    mutationPropGrowth));
## 	    DP2(mutationFromParent(mu, tmpParam, popParams[nextMutant],
## 				   newMutations, mutationPropGrowth));
## 	    DP1("tmpParam");
## 	    print_spP(tmpParam);
## 	    DP1("nextmutatn")
## 	    print_spP(popParams[nextMutant]);
	    
## 	    DP1("end at mutation");
cat("\n Done at", date(), "\n") ## whole file takes about 6 seconds





## But in the tests below, as we increase finalTime, those with
## mutPropGrowth grow to larger population sizes and thus that in itself
## could explain differences in mutations, etc.

## date()
## test_that(" And mutProGrowth, 1",{
##     pseed <- sample(9999999, 1)
##     set.seed(pseed)
##     cat("\n sz0331: the seed is", pseed, "\n")
##     muvar2 <- c("U" = 1e-5, "z" = 1e-5, "e" = 1e-5, "m" = 1e-5, "D" = 1e-5)
##     ni1 <- rep(0.9, 5)
##     names(ni1) <- names(muvar2)
##     fe1 <- allFitnessEffects(noIntGenes = ni1)
##     no <- 1e4
##     reps <- 50
##     ft <- 20
##     pseed <- sample(9999999, 1)
##     set.seed(pseed)
##     cat("\n sz0331a: the seed is", pseed, "\n")
##     b1 <- oncoSimulPop(reps,
##                        fe1, mu = muvar2,
##                        mutationPropGrowth = FALSE,
##                        onlyCancer = FALSE,
##                        initSize = no,
##                        finalTime = ft,
##                        sampleEvery = 0.01,
##                        seed =NULL,
##                        mc.cores = 2
##                        )
##     pseed <- sample(9999999, 1)
##     set.seed(pseed)
##     cat("\n sz0331b: the seed is", pseed, "\n")
##     b2 <- oncoSimulPop(reps,
##                        fe1, mu = muvar2,
##                        mutationPropGrowth = TRUE,
##                        onlyCancer = FALSE,
##                        initSize = no,
##                        finalTime = ft,
##                        sampleEvery = 0.01,                       
##                        seed =NULL,
##                        mc.cores = 2
##                        )
##     ## summary(b2)[, c(1:3, 8:9)]
##     ## mean(mutsPerClone(b1));mean(mutsPerClone(b2))
##     ## median(summary(b1)$NumClones)
##     ## median(summary(b2)$NumClones)
##     ## More mutations in mutationPropGrowth
##     expect_true( mean(mutsPerClone(b2)) >
##                  mean(mutsPerClone(b1)))
##     expect_true( median(summary(b2)$NumClones) >
##                  median(summary(b1)$NumClones))
##     ## But frequency of mutations about the same? Nope: since very fast
##     ## growth and thus non-indep, huge variation in geneCounts in each
##     ## run, etc. so hard to compare geneCounts
##     ## Just for reference, here
##     ## First, look at run to run variation
##     ## OncoSimulR:::geneCounts(b1)
##     ## OncoSimulR:::geneCounts(b2)
##     ## ## The next makes sense
##     ## fb1 <- colSums(OncoSimulR:::geneCounts(b1))
##     ## fb2 <- colSums(OncoSimulR:::geneCounts(b2))
##     ## fb1
##     ## fb2
##     ## fb2/fb1
##     ## fb1/sum(fb1)
##     ## fb2/sum(fb2)
##     ## (fb2/sum(fb2))/(fb1/sum(fb1))
##     ## summary(b2)[, c(1:3, 8:9)]
##     ## mean(mutsPerClone(b1));mean(mutsPerClone(b2))
##     ## median(summary(b1)$NumClones)
##     ## median(summary(b2)$NumClones)
## })
## date()





## date()
## test_that(" And mutProGrowth, 2",{
##     pseed <- sample(9999999, 1)
##     set.seed(pseed)
##     cat("\n sz033: the seed is", pseed, "\n")
##     muvar2 <- c("U" = 1e-4, "z" = 5e-5, "e" = 5e-4, "m" = 5e-3, "D" = 1e-4)
##     ## muvar2 <- c("U" = 5e-5, "z" = 5e-5, "e" = 5e-5, "m" = 5e-5, "D" = 5e-5)
##     ni1 <- rep(1.9, 5)
##     names(ni1) <- names(muvar2)
##     fe1 <- allFitnessEffects(noIntGenes = ni1)
##     no <- 1e5
##     reps <- 10 
##     ft <- 6
##     pseed <- sample(9999999, 1)
##     set.seed(pseed)
##     cat("\n sz033a: the seed is", pseed, "\n")
##     b1 <- oncoSimulPop(reps,
##                        fe1, mu = muvar2,
##                        mutationPropGrowth = FALSE,
##                        onlyCancer = FALSE,
##                        initSize = no,
##                        finalTime = ft,
##                        sampleEvery = 0.01,
##                        seed =NULL,
##                        mc.cores = 2
##                        )
##     pseed <- sample(9999999, 1)
##     set.seed(pseed)
##     cat("\n sz033b: the seed is", pseed, "\n")
##     b2 <- oncoSimulPop(reps,
##                        fe1, mu = muvar2,
##                        mutationPropGrowth = TRUE,
##                        onlyCancer = FALSE,
##                        initSize = no,
##                        finalTime = ft,
##                        sampleEvery = 0.01,                       
##                        seed =NULL,
##                        mc.cores = 2
##                        )
##     summary(b2)[, c(1:3, 8:9)]
##     mean(mutsPerClone(b1));mean(mutsPerClone(b2))
##     median(summary(b1)$NumClones)
##     median(summary(b2)$NumClones)
##     ## More mutations in mutationPropGrowth
##     expect_true( mean(mutsPerClone(b2)) >
##                  mean(mutsPerClone(b1)))
##     expect_true( median(summary(b2)$NumClones) >
##                  median(summary(b1)$NumClones))
##     ## But frequency of mutations about the same? Nope: since very fast
##     ## growth and thus non-indep, huge variation in geneCounts in each
##     ## run, etc. so hard to compare geneCounts
##     ## Just for reference, here
##     ## First, look at run to run variation
##     ## OncoSimulR:::geneCounts(b1)
##     ## OncoSimulR:::geneCounts(b2)
##     ## ## The next makes sense
##     ## fb1 <- colSums(OncoSimulR:::geneCounts(b1))
##     ## fb2 <- colSums(OncoSimulR:::geneCounts(b2))
##     ## fb1
##     ## fb2
##     ## fb2/fb1
##     ## fb1/sum(fb1)
##     ## fb2/sum(fb2)
##     ## (fb2/sum(fb2))/(fb1/sum(fb1))
## })
## date()


## date()
## test_that(" McFL: And mutProGrowth, 1",{
##     pseed <- sample(9999999, 1)
##     set.seed(pseed)
##     cat("\n Mcsz0331: the seed is", pseed, "\n")
##     muvar2 <- c("U" = 1e-5, "z" = 1e-5, "e" = 1e-5, "m" = 1e-5, "D" = 1e-5)
##     ni1 <- rep(0.9, 5)
##     names(ni1) <- names(muvar2)
##     fe1 <- allFitnessEffects(noIntGenes = ni1)
##     no <- 1e4
##     reps <- 50
##     ft <- 20
##     pseed <- sample(9999999, 1)
##     set.seed(pseed)
##     cat("\n Mcsz0331a: the seed is", pseed, "\n")
##     b1 <- oncoSimulPop(reps,
##                        fe1, mu = muvar2,
##                        mutationPropGrowth = FALSE,
##                        onlyCancer = FALSE,
##                        initSize = no,
##                        finalTime = ft,
##                        sampleEvery = 0.01,
##                        seed =NULL,
##                        model = "McFL",
##                        mc.cores = 2
##                        )
##     pseed <- sample(9999999, 1)
##     set.seed(pseed)
##     cat("\n Mcsz0331b: the seed is", pseed, "\n")
##     b2 <- oncoSimulPop(reps,
##                        fe1, mu = muvar2,
##                        mutationPropGrowth = TRUE,
##                        onlyCancer = FALSE,
##                        initSize = no,
##                        finalTime = ft,
##                        sampleEvery = 0.01,                       
##                        seed =NULL,
##                        model = "McFL",
##                        mc.cores = 2
##                        )
##     ## summary(b2)[, c(1:3, 8:9)]
##     ## mean(mutsPerClone(b1));mean(mutsPerClone(b2))
##     ## median(summary(b1)$NumClones)
##     ## median(summary(b2)$NumClones)
##     ## More mutations in mutationPropGrowth
##     expect_true( mean(mutsPerClone(b2)) >
##                  mean(mutsPerClone(b1)))
##     expect_true( median(summary(b2)$NumClones) >
##                  median(summary(b1)$NumClones))
##     ## But frequency of mutations about the same? Nope: since very fast
##     ## growth and thus non-indep, huge variation in geneCounts in each
##     ## run, etc. so hard to compare geneCounts
##     ## Just for reference, here
##     ## First, look at run to run variation
##     ## OncoSimulR:::geneCounts(b1)
##     ## OncoSimulR:::geneCounts(b2)
##     ## ## The next makes sense
##     ## fb1 <- colSums(OncoSimulR:::geneCounts(b1))
##     ## fb2 <- colSums(OncoSimulR:::geneCounts(b2))
##     ## fb1
##     ## fb2
##     ## fb2/fb1
##     ## fb1/sum(fb1)
##     ## fb2/sum(fb2)
##     ## (fb2/sum(fb2))/(fb1/sum(fb1))
##     ## summary(b2)[, c(1:3, 8:9)]
##     ## mean(mutsPerClone(b1));mean(mutsPerClone(b2))
##     ## median(summary(b1)$NumClones)
##     ## median(summary(b2)$NumClones)
## })
## date()





## date()
## test_that(" McFL: And mutProGrowth, 2",{
##     pseed <- sample(9999999, 1)
##     set.seed(pseed)
##     cat("\n mcsz033: the seed is", pseed, "\n")
##     muvar2 <- c("U" = 1e-4, "z" = 5e-5, "e" = 5e-4, "m" = 5e-3, "D" = 1e-4)
##     ## muvar2 <- c("U" = 5e-5, "z" = 5e-5, "e" = 5e-5, "m" = 5e-5, "D" = 5e-5)
##     ni1 <- rep(1.9, 5)
##     names(ni1) <- names(muvar2)
##     fe1 <- allFitnessEffects(noIntGenes = ni1)
##     no <- 1e5
##     reps <- 10 
##     ft <- 6
##     pseed <- sample(9999999, 1)
##     set.seed(pseed)
##     cat("\n mcsz033a: the seed is", pseed, "\n")
##     b1 <- oncoSimulPop(reps,
##                        fe1, mu = muvar2,
##                        mutationPropGrowth = FALSE,
##                        onlyCancer = FALSE,
##                        initSize = no,
##                        finalTime = ft,
##                        sampleEvery = 0.01,
##                        model = "McFL",
##                        seed =NULL,
##                        mc.cores = 2
##                        )
##     pseed <- sample(9999999, 1)
##     set.seed(pseed)
##     cat("\n mcsz033b: the seed is", pseed, "\n")
##     b2 <- oncoSimulPop(reps,
##                        fe1, mu = muvar2,
##                        mutationPropGrowth = TRUE,
##                        onlyCancer = FALSE,
##                        initSize = no,
##                        finalTime = ft,
##                        sampleEvery = 0.01,                       
##                        seed =NULL,
##                        model = "McFL",
##                        mc.cores = 2
##                        )
##     summary(b2)[, c(1:3, 8:9)]
##     mean(mutsPerClone(b1));mean(mutsPerClone(b2))
##     median(summary(b1)$NumClones)
##     median(summary(b2)$NumClones)
##     ## More mutations in mutationPropGrowth
##     expect_true( mean(mutsPerClone(b2)) >
##                  mean(mutsPerClone(b1)))
##     expect_true( median(summary(b2)$NumClones) >
##                  median(summary(b1)$NumClones))
##     ## But frequency of mutations about the same? Nope: since very fast
##     ## growth and thus non-indep, huge variation in geneCounts in each
##     ## run, etc. so hard to compare geneCounts
##     ## Just for reference, here
##     ## First, look at run to run variation
##     ## OncoSimulR:::geneCounts(b1)
##     ## OncoSimulR:::geneCounts(b2)
##     ## ## The next makes sense
##     ## fb1 <- colSums(OncoSimulR:::geneCounts(b1))
##     ## fb2 <- colSums(OncoSimulR:::geneCounts(b2))
##     ## fb1
##     ## fb2
##     ## fb2/fb1
##     ## fb1/sum(fb1)
##     ## fb2/sum(fb2)
##     ## (fb2/sum(fb2))/(fb1/sum(fb1))
## })
## date()


## The problem below is that there are also large differences in
## population size, so the differences in number of clones, etc,
## attributable to that and not just mutationPropGrwoth = TRUE. Which, of
## course, must be behind the differences in popSize, but that is not what
## we are testing here.

## date()
## test_that(" oncoSimuSample and mutPropGrowth",{
##     pseed <- sample(9999999, 1)
##     set.seed(pseed)
##     cat("\n sz0331: the seed is", pseed, "\n")
##     muvar2 <- c("U" = 5e-5, "z" = 5e-6, "e" = 1e-4, "m" = 1e-5, "D" = 5e-4)
##     ni1 <- rep(0.9, 5)
##     names(ni1) <- names(muvar2)
##     fe1 <- allFitnessEffects(noIntGenes = ni1)
##     no <- 1e3  ## very small, and in all
##     reps <- 100
##     ft <- 12
##     x <- 1e-20 # 0.5 * no ## for detection threshold
##     pseed <- sample(9999999, 1)
##     set.seed(pseed)
##     cat("\n sz0331a: the seed is", pseed, "\n")
##     b1 <- oncoSimulSample(reps,
##                           fe1, mu = muvar2,
##                           mutationPropGrowth = FALSE,
##                           onlyCancer = FALSE,
##                           initSize = no,
##                           finalTime = ft,
##                           sampleEvery = 0.01,
##                           detectionSize = 1e9,
##                           detectionDrivers = 99,
##                           seed =NULL,
##                           thresholdWhole = x
##                           )
##     pseed <- sample(9999999, 1)
##     set.seed(pseed)
##     cat("\n sz0331b: the seed is", pseed, "\n")
##     b2 <- oncoSimulSample(reps,
##                        fe1, mu = muvar2,
##                        mutationPropGrowth = TRUE,
##                        onlyCancer = FALSE,
##                        initSize = no,
##                        finalTime = ft,
##                        sampleEvery = 0.01,
##                        detectionSize = 1e9,
##                        detectionDrivers = 99,
##                        seed =NULL,
##                        thresholdWhole = x
##                        )
##     summary(b2$popSummary[, "NumClones"])
##     summary(b1$popSummary[, "NumClones"])
##     summary(b2$popSummary[, "TotalPopSize"])
##     summary(b1$popSummary[, "TotalPopSize"])
##     b2$popSummary[1:5, c(1:3, 8:9)]
##     b1$popSummary[1:5, c(1:3, 8:9)]
    
##     (cc1 <- colSums(b1$popSample))
##     (cc2 <- colSums(b2$popSample))
##     (mutsPerClone1 <- rowSums(b1$popSample))
##     (mutsPerClone2 <- rowSums(b2$popSample))

##   ## I stop about here; diffs in numclones are obvious. But what about the rest?  
    
##     expect_true(sum(cc2) > sum(cc1))
##     expect_true( mean(mutsPerClone2) >
##                  mean(mutsPerClone1))
##     expect_true( median(b2$popSummary[, "NumClones"]) >
##                  median(b1$popSummary[, "NumClones"]))
## })
## date()


## date()
## test_that(" McFL: oncoSimuSample and mutPropGrowth",{

##     pseed <- sample(9999999, 1)
##     set.seed(pseed)
##     cat("\n sz0331: the seed is", pseed, "\n")
##     muvar2 <- c("U" = 5e-5, "z" = 5e-6, "e" = 1e-4, "m" = 1e-5, "D" = 5e-4)
##     ni1 <- rep(0.9, 5)
##     names(ni1) <- names(muvar2)
##     fe1 <- allFitnessEffects(noIntGenes = ni1)
##     no <- 1e3  ## very small, and in all
##     reps <- 100
##     ft <- 12
##     x <- 0.5 * no ## for detection threshold
##     pseed <- sample(9999999, 1)
##     set.seed(pseed)
##     cat("\n sz0331a: the seed is", pseed, "\n")
##     b1 <- oncoSimulSample(reps,
##                           fe1, mu = muvar2,
##                           mutationPropGrowth = FALSE,
##                           onlyCancer = FALSE,
##                           initSize = no,
##                           finalTime = ft,
##                           sampleEvery = 0.01,
##                           detectionSize = 1e9,
##                           detectionDrivers = 99,
##                           seed =NULL,
##                           thresholdWhole = x
##                           )
##     pseed <- sample(9999999, 1)
##     set.seed(pseed)
##     cat("\n sz0331b: the seed is", pseed, "\n")
##     b2 <- oncoSimulSample(reps,
##                        fe1, mu = muvar2,
##                        mutationPropGrowth = TRUE,
##                        onlyCancer = FALSE,
##                        initSize = no,
##                        finalTime = ft,
##                        sampleEvery = 0.01,
##                        detectionSize = 1e9,
##                        detectionDrivers = 99,
##                        seed =NULL,
##                        thresholdWhole = x
##                        )
##     summary(b2$popSummary[, "NumClones"])
##     summary(b1$popSummary[, "NumClones"])
##     summary(b2$popSummary[, "TotalPopSize"])
##     summary(b1$popSummary[, "TotalPopSize"])
##     b2$popSummary[1:5, c(1:3, 8:9)]
##     b1$popSummary[1:5, c(1:3, 8:9)]
    
##     (cc1 <- colSums(b1$popSample))
##     (cc2 <- colSums(b2$popSample))
##     (mutsPerClone1 <- rowSums(b1$popSample))
##     (mutsPerClone2 <- rowSums(b2$popSample))

##   ## I stop about here; diffs in numclones are obvious. But what about the rest?  
    
##     expect_true(sum(cc2) > sum(cc1))
##     expect_true( mean(mutsPerClone2) >
##                  mean(mutsPerClone1))
##     expect_true( median(b2$popSummary[, "NumClones"]) >
##                  median(b1$popSummary[, "NumClones"]))
    

## })
## date()

