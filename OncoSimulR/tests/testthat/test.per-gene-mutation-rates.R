RNGkind("L'Ecuyer-CMRG") ## for the mclapplies

### FIXME: recheck the time for the chisq tests. And beware with the
### seeds

### FIXME: check failures in Lacerta
## I fix seeds because otherwise with continuous testing I get lots of
## false positives.

test_that("single named gene in mut. fail", {
    fea9 <- allFitnessEffects(noIntGenes = c("m" = 0.1, "D" = 0.3))
    muvar <- c("m" = 1e-5)
    expect_error(oncoSimulIndiv(fea9, mu = muvar),
                  "A length 1 mutation, but named",
                  fixed = TRUE)
} )

test_that("Per-gene mutation rates with old poset format, fail", {
    data(examplePosets)
    p701 <- examplePosets[["p701"]]
    muvar <- c(rep(1e-5, 4), rep(1e-6, 3))
    names(muvar) <- letters[1:7]
    expect_error(oncoSimulIndiv(p701, mu = muvar),
                 "Per-gene mutation rates cannot be used with the old poset format")
} )


test_that("Only no-int, and sorting", {
    fea9 <- allFitnessEffects(noIntGenes = c("m" = 0.1, "D" = 0.3))
    ## OncoSimulR:::allNamedGenes(fea9)
    muvar <- c("m" = 1e-5, "D" = 1e-7)
    expect_output(oncoSimulIndiv(fea9, mu = muvar),
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
    expect_output(oncoSimulIndiv(fea8, mu = muvar2),
                  "Individual OncoSimul trajectory",
                  fixed = TRUE)
} )

test_that("Only no-int, unnamed, fail", {
    fea9 <- allFitnessEffects(noIntGenes = c(0.1, 0.3))
    OncoSimulR:::allNamedGenes(fea9)
    muvar <- c(1e-5, 1e-7)
    expect_error(oncoSimulIndiv(fea9, mu = muvar),
                 "When using per-gene mutation rates the mu vector must be named",
                 fixed = TRUE)
} )

test_that("Only one, named, fail", {
    fea9 <- allFitnessEffects(noIntGenes = c("m" = 0.1))
    muvar <- c("m" = 1e-5)
    expect_error(oncoSimulIndiv(fea9, mu = muvar),
                 "There must be at least two genes (loci) in the fitness effects",
                 fixed = TRUE)
} )

test_that("Only no-int, different names, fail", {
    fea9 <- allFitnessEffects(noIntGenes = c("m" = 0.1, "D" = 0.3))
    ## OncoSimulR:::allNamedGenes(fea9)
    muvar <- c("n" = 1e-5, "D" = 1e-7)
    expect_error(oncoSimulIndiv(fea9, mu = muvar),
                 "When using per-gene mutation rates, names of genes must match",
                  fixed = TRUE)
} )

test_that("Only no-int, different numbers, fail", {
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


test_that("Same freqs, chisq, when s=0 and t = 1", {
  
    
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
    colSums(OncoSimulR:::geneCounts(bb))
    chisq.test(colSums(OncoSimulR:::geneCounts(bb)))
        ## It will fail with prob ~ p.fail
    p.fail <- 1e-3
    expect_true(chisq.test(colSums(OncoSimulR:::geneCounts(bb)),
                           p = expectedC/sum(expectedC))$p.value > p.fail)


})


test_that("Same freqs, chisq, when s and t = 1", {
    muvar2 <- c("U" = 1e-5, "z" = 1e-5, "e" = 1e-5, "m" = 1e-5, "D" = 1e-5)
    ni1 <- rep(0.02, 5)
    names(ni1) <- names(muvar2)
    fe1 <- allFitnessEffects(noIntGenes = ni1)
    no <- 1e5
    reps <- 10000
    set.seed(134)
    bb <- oncoSimulPop(reps,
                       fe1, mu = muvar2, onlyCancer = FALSE,
                       initSize = no,
                       finalTime = 1,
                       seed = NULL, mc.cores = 2
                       )
    (expectedC <- no*reps*muvar2)
    colSums(OncoSimulR:::geneCounts(bb))
    chisq.test(colSums(OncoSimulR:::geneCounts(bb)))
        ## It will fail with prob ~ p.fail
    p.fail <- 1e-3
    expect_true(chisq.test(colSums(OncoSimulR:::geneCounts(bb)),
                           p = expectedC/sum(expectedC))$p.value > p.fail)
})





test_that("Different freqs as they should be ordered and chisq, when s=0 and t = 1",
{
    muvar2 <- c("U" = 1e-3, "z" = 1e-7, "e" = 1e-6, "m" = 1e-5, "D" = 1e-4)
    ni1 <- rep(0, 5)
    names(ni1) <- names(muvar2)
    fe1 <- allFitnessEffects(noIntGenes = ni1)
    no <- 5e5
    reps <- 1000
    set.seed(703)
    bb <- oncoSimulPop(reps,
                       fe1, mu = muvar2, onlyCancer = FALSE,
                       initSize = no,
                       finalTime = 1,
                       seed =NULL
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


test_that("Different freqs as they should be ordered and chisq, when s and t = 1", {
    muvar2 <- c("U" = 1e-3, "z" = 1e-7, "e" = 1e-6, "m" = 1e-5, "D" = 1e-4)
    ni2 <- rep(0.01, 5)
    names(ni2) <- names(muvar2)
    fe1 <- allFitnessEffects(noIntGenes = ni2)
    no <- 5e5
    reps <- 1000
    set.seed(987)
    bb <- oncoSimulPop(reps,
                       fe1, mu = muvar2, onlyCancer = FALSE,
                       initSize = no,
                       finalTime = 1,
                       seed = NULL, mc.cores = 2
                       )
    (expectedC <- no*reps*muvar2)
    colSums(OncoSimulR:::geneCounts(bb))
    ## A chisq will not work as we increase finalTime. But ordering of
    ## freqs. should.
    ## This will fail sometimes
    p.fail <- 1e-3
    expect_true(chisq.test(colSums(OncoSimulR:::geneCounts(bb)),
                           p = expectedC/sum(expectedC))$p.value > p.fail)
    expect_equal(
        order(colSums(OncoSimulR:::geneCounts(bb))),
        order(expectedC))
})


test_that("Different freqs as they should be ordered, when s and t> 1", {
    muvar2 <- c("U" = 1e-3, "z" = 1e-7, "e" = 1e-6, "m" = 1e-5, "D" = 1e-4)
    ni2 <- rep(0.01, 5)
    names(ni2) <- names(muvar2)
    fe1 <- allFitnessEffects(noIntGenes = ni2)
    no <- 5e5
    reps <- 1000
    set.seed(123456)
    bb <- oncoSimulPop(reps,
                       fe1, mu = muvar2, onlyCancer = FALSE,
                       initSize = no,
                       finalTime = 4,
                       seed = NULL, mc.cores = 2
                       )
    (expectedC <- no*reps*muvar2)
    colSums(OncoSimulR:::geneCounts(bb))
    ## A chisq will not work as we increase finalTime. But ordering of
    ## freqs. should.
    ## chisq.test(colSums(OncoSimulR:::geneCounts(bb)),
    ##            p = expectedC/sum(expectedC))
    expect_equal(
        order(colSums(OncoSimulR:::geneCounts(bb))),
        order(expectedC))
})



test_that("Complex fitness specification, tiny s diffs", {
    set.seed(1)
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
    ## muvar <- runif(length(nfea), min = 1e-7, max = 1e-3) ## too tiny
    ## diffs sometimes for order comp
    muvar <- sample(seq(from = 1e-7, to = 1e-3, length.out = length(nfea)))
    names(muvar) <- nfea
    no <- 1e4
    reps <- 5000
    bb <- oncoSimulPop(reps,
                       fea, mu = muvar,
                       onlyCancer = FALSE,
                       initSize = no,
                       finalTime = 1,
                       seed = NULL, mc.cores = 2
                       )
    (expectedC <- no*reps*muvar)
    colSums(OncoSimulR:::geneCounts(bb))
    expectedC - colSums(OncoSimulR:::geneCounts(bb))
    ## chisq.test(colSums(OncoSimulR:::geneCounts(bb)),
    ##                        p = expectedC/sum(expectedC))
    p.fail <- 1e-3
    expect_true(chisq.test(colSums(OncoSimulR:::geneCounts(bb)),
                           p = expectedC/sum(expectedC))$p.value > p.fail)
    expect_equal(
        order(colSums(OncoSimulR:::geneCounts(bb))),
        order(expectedC))
})


test_that("Complex fitness specification, large diffs", {
    set.seed(89)
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
    muvar <- sample(seq(from = 1e-7, to = 1e-3, length.out = length(nfea)))
    names(muvar) <- nfea
    no <- 1e4
    reps <- 2000
    bb <- oncoSimulPop(reps,
                       fea, mu = muvar,
                       onlyCancer = FALSE,
                       initSize = no,
                       finalTime = 1,
                       seed = NULL, mc.cores = 2
                       )
    (expectedC <- no*reps*muvar)
    colSums(OncoSimulR:::geneCounts(bb))
    expectedC - colSums(OncoSimulR:::geneCounts(bb))
    ## this will occasionally fail
    p.fail <- 1e-6
    expect_true(chisq.test(colSums(OncoSimulR:::geneCounts(bb)),
                           p = expectedC/sum(expectedC))$p.value < p.fail)
    expect_false(
        identical(
            order(colSums(OncoSimulR:::geneCounts(bb))),
            order(expectedC)))
})


#### Repeating above, but with McFL


test_that("McFL: Per-gene mutation rates with old poset format, fail", {
    data(examplePosets)
    p701 <- examplePosets[["p701"]]
    muvar <- c(rep(1e-5, 4), rep(1e-6, 3))
    names(muvar) <- letters[1:7]
    expect_error(oncoSimulIndiv(p701, mu = muvar, model = "McFL"),
                 "Per-gene mutation rates cannot be used with the old poset format")
} )


test_that("McFL: Only no-int, and sorting", {
    fea9 <- allFitnessEffects(noIntGenes = c("m" = 0.1, "D" = 0.3))
    ## OncoSimulR:::allNamedGenes(fea9)
    muvar <- c("m" = 1e-5, "D" = 1e-7)
    expect_output(oncoSimulIndiv(fea9, mu = muvar, model = "McFL", finalTime = 20),
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
    expect_output(oncoSimulIndiv(fea8, mu = muvar2, model = "McFL", finalTime = 20),
                  "Individual OncoSimul trajectory",
                  fixed = TRUE)
} )

test_that("McFL: Only no-int, unnamed, fail", {
    fea9 <- allFitnessEffects(noIntGenes = c(0.1, 0.3))
    OncoSimulR:::allNamedGenes(fea9)
    muvar <- c(1e-5, 1e-7)
    expect_error(oncoSimulIndiv(fea9, mu = muvar, model = "McFL",
                                finalTime = 20),
                 "When using per-gene mutation rates the mu vector must be named",
                 fixed = TRUE)
} )

test_that("McFL: Only one, named, fail", {
    fea9 <- allFitnessEffects(noIntGenes = c("m" = 0.1))
    muvar <- c("m" = 1e-5)
    expect_error(oncoSimulIndiv(fea9, mu = muvar, model = "McFL",
                                finalTime = 20),
                 "There must be at least two genes (loci) in the fitness effects",
                 fixed = TRUE)
} )

test_that("McFL: Only no-int, different names, fail", {
    fea9 <- allFitnessEffects(noIntGenes = c("m" = 0.1, "D" = 0.3))
    ## OncoSimulR:::allNamedGenes(fea9)
    muvar <- c("n" = 1e-5, "D" = 1e-7)
    expect_error(oncoSimulIndiv(fea9, mu = muvar, model = "McFL",
                                finalTime = 20),
                 "When using per-gene mutation rates, names of genes must match",
                  fixed = TRUE)
} )

test_that("McFL: Only no-int, different numbers, fail", {
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

test_that("McFL: Same freqs, chisq, when s=0 and t = 1", {
    muvar2 <- c("U" = 1e-5, "z" = 1e-5, "e" = 1e-5, "m" = 1e-5, "D" = 1e-5)
    ni1 <- rep(0, 5)
    names(ni1) <- names(muvar2)
    fe1 <- allFitnessEffects(noIntGenes = ni1)
    no <- 1e5
    reps <- 5000
    set.seed(65098)
    bb <- oncoSimulPop(reps,
                       fe1, mu = muvar2, onlyCancer = FALSE,
                       model = "McFL",
                       initSize = no,
                       finalTime = 1,
                       seed= NULL
                       )
    (expectedC <- no*reps*muvar2)
    colSums(OncoSimulR:::geneCounts(bb))
    chisq.test(colSums(OncoSimulR:::geneCounts(bb)))
    ## It will fail with prob ~ p.fail
    p.fail <- 1e-3
    expect_true(chisq.test(colSums(OncoSimulR:::geneCounts(bb)),
                           p = expectedC/sum(expectedC))$p.value > p.fail)
})


test_that("McFL: Same freqs, chisq, when s and t = 1", {
    set.seed(3018)
    muvar2 <- c("U" = 1e-5, "z" = 1e-5, "e" = 1e-5, "m" = 1e-5, "D" = 1e-5)
    ni1 <- rep(0.02, 5)
    names(ni1) <- names(muvar2)
    fe1 <- allFitnessEffects(noIntGenes = ni1)
    no <- 1e5
    reps <- 2000
    bb <- oncoSimulPop(reps,
                       fe1, mu = muvar2,
                       model = "McFL",
                       onlyCancer = FALSE,
                       initSize = no,
                       finalTime = 1,
                       seed = NULL, mc.cores = 2
                       )
    (expectedC <- no*reps*muvar2)
    colSums(OncoSimulR:::geneCounts(bb))
    chisq.test(colSums(OncoSimulR:::geneCounts(bb)))
    ## It will fail with prob ~ p.fail
    p.fail <- 1e-3
    expect_true(chisq.test(colSums(OncoSimulR:::geneCounts(bb)),
                           p = expectedC/sum(expectedC))$p.value > p.fail)
})

test_that("McFL: Different freqs as they should be ordered and chisq, when s=0 and t = 1", {
    muvar2 <- c("U" = 1e-3, "z" = 1e-7, "e" = 1e-6, "m" = 1e-5, "D" = 1e-4)
    ni1 <- rep(0, 5)
    names(ni1) <- names(muvar2)
    fe1 <- allFitnessEffects(noIntGenes = ni1)
    no <- 2e5
    reps <- 5000
    set.seed(8)
    bb <- oncoSimulPop(reps,
                       fe1, mu = muvar2, onlyCancer = FALSE,
                       model = "McFL",
                       initSize = no,
                       finalTime = 1,
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


test_that("McFL: Different freqs as they should be ordered and chisq, when s and t = 1", {
    muvar2 <- c("U" = 1e-3, "z" = 1e-7, "e" = 1e-6, "m" = 1e-5, "D" = 1e-4)
    ni2 <- rep(0.01, 5)
    names(ni2) <- names(muvar2)
    fe1 <- allFitnessEffects(noIntGenes = ni2)
    no <- 2e5
    reps <- 2000
    set.seed(65)
    bb <- oncoSimulPop(reps,
                       fe1, mu = muvar2, onlyCancer = FALSE,
                       model = "McFL",
                       initSize = no,
                       finalTime = 1,
                       seed = NULL, mc.cores = 2
                       )
    (expectedC <- no*reps*muvar2)
    colSums(OncoSimulR:::geneCounts(bb))
    ## A chisq will not work as we increase finalTime. But ordering of
    ## freqs. should.
    ## This will fail sometimes
    p.fail <- 1e-3
    expect_true(chisq.test(colSums(OncoSimulR:::geneCounts(bb)),
                           p = expectedC/sum(expectedC))$p.value > p.fail)
    expect_equal(
        order(colSums(OncoSimulR:::geneCounts(bb))),
        order(expectedC))
})


test_that("McFL: Different freqs as they should be ordered, when s and t> 1", {
    set.seed(1)
    muvar2 <- c("U" = 1e-3, "z" = 1e-7, "e" = 1e-6, "m" = 1e-5, "D" = 1e-4)
    ni2 <- rep(0.01, 5)
    names(ni2) <- names(muvar2)
    fe1 <- allFitnessEffects(noIntGenes = ni2)
    no <- 1e5
    reps <- 1000
    bb <- oncoSimulPop(reps,
                       fe1, mu = muvar2, onlyCancer = FALSE,
                       model = "McFL",
                       initSize = no,
                       finalTime = 4,
                       seed = NULL, mc.cores = 2
                       )
    (expectedC <- no*reps*muvar2)
    colSums(OncoSimulR:::geneCounts(bb))
    ## A chisq will not work as we increase finalTime. But ordering of
    ## freqs. should.
    ## chisq.test(colSums(OncoSimulR:::geneCounts(bb)),
    ##            p = expectedC/sum(expectedC))
    expect_equal(
        order(colSums(OncoSimulR:::geneCounts(bb))),
        order(expectedC))
})



test_that("McFL: Complex fitness specification, tiny s diffs", {
    set.seed(1) 
    p4 <- data.frame(parent = c(rep("Root", 4), "A", "B", "D", "E", "C", "F"),
                     child = c("A", "B", "D", "E", "C", "C", "F", "F", "G", "G"),
                     s = c(0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001),
                     sh = c(rep(0, 4), c(-.00000009, -.00000009), c(-.000000095, -.000000095), c(-.0000099, -.0000099)),
                     typeDep = c(rep("--", 4), 
                                 "XMPN", "XMPN", "MN", "MN", "SM", "SM"))
    oe <- c("C > F" = -0.0000001, "H > I" = 0.00000012)
    sm <- c("I:J"  = -.0000001)
    sv <- c("-K:M" = -.00000145, "K:-M" = -.0000015)
    epist <- c(sm, sv)
    modules <- c("Root" = "Root", "A" = "a1",
                 "B" = "b1, b2", "C" = "c1",
                 "D" = "d1, d2", "E" = "e1",
                 "F" = "f1, f2", "G" = "g1",
                 "H" = "h1, h2", "I" = "i1",
                 "J" = "j1, j2", "K" = "k1, k2", "M" = "m1")
    noint <- runif(5, min = 0.0000051, max = 0.0000081)
    names(noint) <- paste0("n", 1:5)
    fea <- allFitnessEffects(rT = p4, epistasis = epist, orderEffects = oe,
                             noIntGenes = noint, geneToModule = modules)
    nfea <- OncoSimulR:::allNamedGenes(fea)$Gene
    ## muvar <- runif(length(nfea), min = 1e-7, max = 1e-3) ## too tiny
    ## diffs sometimes for order comp
    muvar <- sample(seq(from = 1e-7, to = 1e-3, length.out = length(nfea)))
    names(muvar) <- nfea
    no <- 1e3
    reps <- 3500
    bb <- oncoSimulPop(reps,
                       fea, mu = muvar,
                       model = "McFL",
                       onlyCancer = FALSE,
                       initSize = no,
                       finalTime = 1,
                       seed = NULL, mc.cores = 2
                       )
    (expectedC <- no*reps*muvar)
    colSums(OncoSimulR:::geneCounts(bb))
    expectedC - colSums(OncoSimulR:::geneCounts(bb))
    ## chisq.test(colSums(OncoSimulR:::geneCounts(bb)),
    ##                        p = expectedC/sum(expectedC))
    p.fail <- 1e-3
    expect_true(chisq.test(colSums(OncoSimulR:::geneCounts(bb)),
                           p = expectedC/sum(expectedC))$p.value > p.fail)
    expect_equal(
        order(colSums(OncoSimulR:::geneCounts(bb))),
        order(expectedC))
})


test_that("McFL: Complex fitness specification, large diffs", {
    set.seed(4)
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
    muvar <- sample(seq(from = 1e-7, to = 1e-3, length.out = length(nfea)))
    names(muvar) <- nfea
    no <- 1e3
    reps <- 1000
    bb <- oncoSimulPop(reps,
                       fea, mu = muvar,
                       model = "McFL",
                       onlyCancer = FALSE,
                       initSize = no,
                       finalTime = 1,
                       seed = NULL, mc.cores = 2
                       )
    (expectedC <- no*reps*muvar)
    colSums(OncoSimulR:::geneCounts(bb))
    expectedC - colSums(OncoSimulR:::geneCounts(bb))
    ## this will occasionally fail
    p.fail <- 1e-6
    expect_true(chisq.test(colSums(OncoSimulR:::geneCounts(bb)),
                           p = expectedC/sum(expectedC))$p.value < p.fail)
    expect_false(
        identical(
            order(colSums(OncoSimulR:::geneCounts(bb))),
            order(expectedC)))
})



test_that("get.gene.counts exercising for NA case", {
    set.seed(1)
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


test_that("Init mutant with tiny mutation", {
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
    no <- 1e2
    reps <- 500
    bb <- oncoSimulPop(reps,
                       fea, mu = muvar,
                       onlyCancer = FALSE,
                       initSize = no,
                       finalTime = 10000,
                       detectionDrivers = 4,
                       mutationPropGrowth = FALSE,
                       initMutant = "h2 > i1"
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


test_that("McFL: Init mutant with tiny mutation", {
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
    no <- 1e3
    reps <- 500
    bb <- oncoSimulPop(5, ##reps,
                       fea, mu = muvar,
                       onlyCancer = FALSE,
                       initSize = no,
                       model = "McFL",
                       finalTime = 600,
                       initMutant = "h2 > i1"
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

test_that("Different freqs as they should be ordered and chisq, when s  and a tiny mu", {
    muvar2 <- c("U" = 1e-13, "z" = 1e-7, "e" = 1e-6, "m" = 1e-5, "D" = 1e-4)
    ni2 <- rep(0.01, 5)
    names(ni2) <- names(muvar2)
    ni2["U"] <- 0.5
    fe1 <- allFitnessEffects(noIntGenes = ni2)
    no <- 5e5
    reps <- 1000
    set.seed(1)
    bb <- oncoSimulPop(reps,
                       fe1, mu = muvar2, onlyCancer = FALSE,
                       initSize = no,
                       finalTime = 5,
                       seed = NULL, mc.cores = 2
                       )
    (expectedC <- no*reps*muvar2)
    colSums(OncoSimulR:::geneCounts(bb))
    expect_true(colSums(OncoSimulR:::geneCounts(bb))[1] == 0)
    expect_equal(
        order(colSums(OncoSimulR:::geneCounts(bb))),
        order(expectedC))
    ## A chisq will not work as we increase finalTime.
    set.seed(2)
    bb <- oncoSimulPop(reps,
                       fe1, mu = muvar2, onlyCancer = FALSE,
                       initSize = no,
                       finalTime = 1,
                       seed = NULL, mc.cores = 2
                       )
    (expectedC <- no*reps*muvar2)
    colSums(OncoSimulR:::geneCounts(bb))
    expect_true(colSums(OncoSimulR:::geneCounts(bb))[1] == 0)
    ## This will fail sometimes
    p.fail <- 1e-3
    expect_true(chisq.test(colSums(OncoSimulR:::geneCounts(bb))[-1],
                           p = expectedC[-1]/sum(expectedC))$p.value > p.fail)
    expect_equal(
        order(colSums(OncoSimulR:::geneCounts(bb))),
        order(expectedC))
}
)

test_that("McFL: Different freqs as they should be ordered and chisq, when s  and a tiny mu", {
    muvar2 <- c("U" = 1e-13, "z" = 5e-7, "e" = 1e-6, "m" = 1e-5, "D" = 1e-4)
    ni2 <- rep(0.01, 5)
    names(ni2) <- names(muvar2)
    ni2["U"] <- 0.5
    fe1 <- allFitnessEffects(noIntGenes = ni2)
    no <- 1e5
    reps <- 100
    ## Beware: with McFL final time has to be small or we will see
    ## mutations in U as it is the only one left. 
    set.seed(26)
    bb <- oncoSimulPop(reps,
                       fe1, mu = muvar2, onlyCancer = FALSE,
                       initSize = no,
                       model = "McFL",
                       finalTime = 200,
                       seed = NULL, mc.cores = 2
                       )
    (expectedC <- no*reps*muvar2)
    colSums(OncoSimulR:::geneCounts(bb))
    expect_true(colSums(OncoSimulR:::geneCounts(bb))[1] == 0)
    expect_equal(
        order(colSums(OncoSimulR:::geneCounts(bb))),
        order(expectedC))
    ## A chisq will not work as we increase finalTime.
    set.seed(35)
    no <- 1e5
    bb <- oncoSimulPop(reps,
                       fe1, mu = muvar2, onlyCancer = FALSE,
                       initSize = no,
                       model = "McFL",
                       finalTime = 1,
                       seed = NULL, mc.cores = 2
                       )
    (expectedC <- no*reps*muvar2)
    colSums(OncoSimulR:::geneCounts(bb))
    expect_true(colSums(OncoSimulR:::geneCounts(bb))[1] == 0)
    ## This will fail sometimes
    p.fail <- 1e-3
    expect_true(chisq.test(colSums(OncoSimulR:::geneCounts(bb))[-1],
                           p = expectedC[-1]/sum(expectedC))$p.value > p.fail)
    expect_equal(
        order(colSums(OncoSimulR:::geneCounts(bb))),
        order(expectedC))
}
)




test_that("0 or negative mu not allowed", {
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





test_that("Different freqs as they should be ordered and chisq, when s=0 and t = 1, and initMutant",
{
    muvar2 <- c("U" = 1e-3, "z" = 1e-7, "e" = 1e-6, "m" = 1e-5, "D" = 1e-4)
    ni1 <- rep(0, 5)
    names(ni1) <- names(muvar2)
    fe1 <- allFitnessEffects(noIntGenes = ni1)
    no <- 5e5
    reps <- 1000
    set.seed(6305)
    bb <- oncoSimulPop(reps,
                       fe1, mu = muvar2, onlyCancer = FALSE,
                       initSize = no,
                       finalTime = 1,
                       seed =NULL,
                       initMutant = "e"
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
    ## chisq.test(colSums(OncoSimulR:::geneCounts(bb)),
    ##            p = expectedC/sum(expectedC))
})


test_that("McFL: Different freqs as they should be ordered and chisq, when s=0 and t = 1, and initMutant",
{
    muvar2 <- c("U" = 1e-3, "z" = 1e-7, "e" = 1e-6, "m" = 1e-5, "D" = 1e-4)
    ni1 <- rep(0, 5)
    names(ni1) <- names(muvar2)
    fe1 <- allFitnessEffects(noIntGenes = ni1)
    no <- 1e5
    reps <- 1000
    set.seed(465)
    bb <- oncoSimulPop(reps,
                       fe1, mu = muvar2, onlyCancer = FALSE,
                       initSize = no,
                       finalTime = 1,
                       seed =NULL,
                       model = "McFL",
                       initMutant = "m"
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





test_that("Different freqs as they are expected with chisq, when s=0 and initMutant, many genotypes",{
    ## Occasionally, the c++ code can blow up as a large mu, etc
    set.seed(6305)
    ft <- 0.2 ## yes, small
    lni <- 100  ## 16
    muvar2 <- runif(lni, min = 1e-9, max = 1e-3)
    names(muvar2) <- c(replicate(lni,
                                 paste(sample(letters, 12), collapse = "")))
    muvar2[2] <- 2e-5 
    names(muvar2)[3] <- "e"
    muvar2[2] <- 1e-3
    ni1 <- rep(0, lni)
    names(ni1) <- names(muvar2)
    fe1 <- allFitnessEffects(noIntGenes = ni1)
    no <- 1e5
    reps <- 500
    bb <- oncoSimulPop(reps,
                       fe1, mu = muvar2, onlyCancer = FALSE,
                       initSize = no,
                       finalTime = ft,
                       seed =NULL,
                       initMutant = "e"
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





## More on the above, with less variation. But yet other set of tests.
test_that("Different freqs as they should be ordered and chisq, when s=0 and t = 1, and initMutant",
{
    muvar2 <- c("U" = 1e-3, "z" = 1e-7, "e" = 1e-6, "m" = 1e-5, "D" = 1e-4)
    ## moderately small mu
    muvar2[] <- 1e-5
    muvar2["e"] <- 1e-3
    ni1 <- rep(0, 5)
    names(ni1) <- names(muvar2)
    fe1 <- allFitnessEffects(noIntGenes = ni1)
    no <- 5e5
    reps <- 1000
    ## set.seed(6305)
    bb <- oncoSimulPop(reps,
                       fe1, mu = muvar2, onlyCancer = FALSE,
                       initSize = no,
                       finalTime = 0.2,
                       seed =NULL,
                       initMutant = "e"
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
    no <- 5e5
    reps <- 500
    ## set.seed(6305)
    bb <- oncoSimulPop(reps,
                       fe1, mu = muvar2, onlyCancer = FALSE,
                       initSize = no,
                       finalTime = 0.2,
                       seed =NULL,
                       initMutant = "e"
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
    no <- 5e5
    reps <- 500
    ## set.seed(6305)
    bb <- oncoSimulPop(reps,
                       fe1, mu = muvar2, onlyCancer = FALSE,
                       initSize = no,
                       finalTime = 0.2,
                       seed =NULL,
                       initMutant = "e"
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

totalind <- function(out) {
    ## total num indivs
  sum(unlist(lapply(out, function(x) x$TotalPopSize)))  
}


## FIXME: new and done

test_that("Expect freqs, num clones, muts per clone for different per-gene-mut",{
    pseed <- sample(1:9999999, 1)
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
    b1 <- oncoSimulPop(reps,
                       fe1,
                       mu = m1,
                       onlyCancer = FALSE,
                       initSize = no,
                       finalTime = ft,
                       seed =NULL
                       )
    b2 <- oncoSimulPop(reps,
                       fe1,
                       mu = m2,
                       onlyCancer = FALSE,
                       initSize = no,
                       finalTime = ft,
                       seed =NULL
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


test_that("Num clones, muts per clone for different per-gene-mut",{
    ## Like previous, but larger finalTime, so no longer chi-square test
    ## here.
    pseed <- sample(1:9999999, 1)
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
    b1 <- oncoSimulPop(reps,
                       fe1,
                       mu = m1,
                       onlyCancer = FALSE,
                       initSize = no,
                       finalTime = ft,
                       seed =NULL
                       )
    b2 <- oncoSimulPop(reps,
                       fe1,
                       mu = m2,
                       onlyCancer = FALSE,
                       initSize = no,
                       finalTime = ft,
                       seed =NULL
                       )
    expect_true( median(summary(b2)$NumClones) >
                 median(summary(b1)$NumClones))
    ## Note the short time, so this is not always very different as few
    ## have double or triple mutants
    expect_true( mean(mutsPerClone(b2)) >
                 mean(mutsPerClone(b1)))
})

## FIXME: new and done end block



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

