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
    no <- 5e5
    reps <- 10000
    bb <- oncoSimulPop(reps,
                       fe1, mu = muvar2, onlyCancer = FALSE,
                       initSize = no,
                       finalTime = 1
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
    bb <- oncoSimulPop(reps,
                       fe1, mu = muvar2, onlyCancer = FALSE,
                       initSize = no,
                       finalTime = 1
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
    muvar2 <- c("U" = 1e-6, "z" = 1e-7, "e" = 1e-3, "m" = 1e-5, "D" = 1e-4)
    ni1 <- rep(0, 5)
    names(ni1) <- names(muvar2)
    fe1 <- allFitnessEffects(noIntGenes = ni1)
    no <- 5e5
    reps <- 1000
    bb <- oncoSimulPop(reps,
                       fe1, mu = muvar2, onlyCancer = FALSE,
                       initSize = no,
                       finalTime = 1
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
    muvar2 <- c("U" = 1e-6, "z" = 1e-7, "e" = 1e-3, "m" = 1e-5, "D" = 1e-4)
    ni2 <- rep(0.01, 5)
    names(ni2) <- names(muvar2)
    fe1 <- allFitnessEffects(noIntGenes = ni2)
    no <- 5e5
    reps <- 1000
    bb <- oncoSimulPop(reps,
                       fe1, mu = muvar2, onlyCancer = FALSE,
                       initSize = no,
                       finalTime = 1
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
    muvar2 <- c("U" = 1e-6, "z" = 1e-7, "e" = 1e-3, "m" = 1e-5, "D" = 1e-4)
    ni2 <- rep(0.01, 5)
    names(ni2) <- names(muvar2)
    fe1 <- allFitnessEffects(noIntGenes = ni2)
    no <- 5e5
    reps <- 1000
    bb <- oncoSimulPop(reps,
                       fe1, mu = muvar2, onlyCancer = FALSE,
                       initSize = no,
                       finalTime = 4
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

    p4 <- data.frame(parent = c(rep("Root", 4), "A", "B", "D", "E", "C", "F"),
                     child = c("A", "B", "D", "E", "C", "C", "F", "F", "G", "G"),
                     s = c(0.001, 0.002, 0.003, 0.004, 0.001, 0.001, 0.002, 0.002, 0.003, 0.003),
                     sh = c(rep(0, 4), c(-.00009, -.00009), c(-.000095, -.000095), c(-.000099, -.000099)),
                     typeDep = c(rep("--", 4), 
                                 "XMPN", "XMPN", "MN", "MN", "SM", "SM"))
    
    oe <- c("C > F" = -0.001, "H > I" = 0.0012)
    sm <- c("I:J"  = -.001)
    sv <- c("-K:M" = -.0005, "K:-M" = -.0005)
    epist <- c(sm, sv)
    modules <- c("Root" = "Root", "A" = "a1",
                 "B" = "b1, b2", "C" = "c1",
                 "D" = "d1, d2", "E" = "e1",
                 "F" = "f1, f2", "G" = "g1",
                 "H" = "h1, h2", "I" = "i1",
                 "J" = "j1, j2", "K" = "k1, k2", "M" = "m1")
    noint <- runif(5, min = 0.00051, max = 0.001)
    names(noint) <- paste0("n", 1:5)
    fea <- allFitnessEffects(rT = p4, epistasis = epist, orderEffects = oe,
                             noIntGenes = noint, geneToModule = modules)
    nfea <- OncoSimulR:::allNamedGenes(fea)$Gene
    ## muvar <- runif(length(nfea), min = 1e-7, max = 1e-3) ## too tiny
    ## diffs sometimes for order comp
    muvar <- sample(seq(from = 1e-7, to = 1e-3, length.out = length(nfea)))
    names(muvar) <- nfea
    no <- 1e4
    reps <- 10000
    bb <- oncoSimulPop(reps,
                       fea, mu = muvar,
                       onlyCancer = FALSE,
                       initSize = no,
                       finalTime = 1
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
                       finalTime = 1
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


test_that("Per-gene mutation rates with old poset format, fail", {
    data(examplePosets)
    p701 <- examplePosets[["p701"]]
    muvar <- c(rep(1e-5, 4), rep(1e-6, 3))
    names(muvar) <- letters[1:7]
    expect_error(oncoSimulIndiv(p701, mu = muvar, model = "McFL"),
                 "Per-gene mutation rates cannot be used with the old poset format")
} )


test_that("Only no-int, and sorting", {
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

test_that("Only no-int, unnamed, fail", {
    fea9 <- allFitnessEffects(noIntGenes = c(0.1, 0.3))
    OncoSimulR:::allNamedGenes(fea9)
    muvar <- c(1e-5, 1e-7)
    expect_error(oncoSimulIndiv(fea9, mu = muvar, model = "McFL"),
                 "When using per-gene mutation rates the mu vector must be named",
                 fixed = TRUE)
} )

test_that("Only one, named, fail", {
    fea9 <- allFitnessEffects(noIntGenes = c("m" = 0.1))
    muvar <- c("m" = 1e-5)
    expect_error(oncoSimulIndiv(fea9, mu = muvar, model = "McFL"),
                 "There must be at least two genes (loci) in the fitness effects",
                 fixed = TRUE)
} )

test_that("Only no-int, different names, fail", {
    fea9 <- allFitnessEffects(noIntGenes = c("m" = 0.1, "D" = 0.3))
    ## OncoSimulR:::allNamedGenes(fea9)
    muvar <- c("n" = 1e-5, "D" = 1e-7)
    expect_error(oncoSimulIndiv(fea9, mu = muvar, model = "McFL"),
                 "When using per-gene mutation rates, names of genes must match",
                  fixed = TRUE)
} )

test_that("Only no-int, different numbers, fail", {
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




##FIXME continue here FIXME
test_that("Same freqs, chisq, when s=0 and t = 1", {
    muvar2 <- c("U" = 1e-5, "z" = 1e-5, "e" = 1e-5, "m" = 1e-5, "D" = 1e-5)
    ni1 <- rep(0, 5)
    names(ni1) <- names(muvar2)
    fe1 <- allFitnessEffects(noIntGenes = ni1)
    no <- 5e5
    reps <- 10000
    bb <- oncoSimulPop(reps,
                       fe1, mu = muvar2, onlyCancer = FALSE,
                       initSize = no,
                       finalTime = 1
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
    bb <- oncoSimulPop(reps,
                       fe1, mu = muvar2, onlyCancer = FALSE,
                       initSize = no,
                       finalTime = 1
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
    muvar2 <- c("U" = 1e-6, "z" = 1e-7, "e" = 1e-3, "m" = 1e-5, "D" = 1e-4)
    ni1 <- rep(0, 5)
    names(ni1) <- names(muvar2)
    fe1 <- allFitnessEffects(noIntGenes = ni1)
    no <- 5e5
    reps <- 1000
    bb <- oncoSimulPop(reps,
                       fe1, mu = muvar2, onlyCancer = FALSE,
                       initSize = no,
                       finalTime = 1
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
    muvar2 <- c("U" = 1e-6, "z" = 1e-7, "e" = 1e-3, "m" = 1e-5, "D" = 1e-4)
    ni2 <- rep(0.01, 5)
    names(ni2) <- names(muvar2)
    fe1 <- allFitnessEffects(noIntGenes = ni2)
    no <- 5e5
    reps <- 1000
    bb <- oncoSimulPop(reps,
                       fe1, mu = muvar2, onlyCancer = FALSE,
                       initSize = no,
                       finalTime = 1
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
    muvar2 <- c("U" = 1e-6, "z" = 1e-7, "e" = 1e-3, "m" = 1e-5, "D" = 1e-4)
    ni2 <- rep(0.01, 5)
    names(ni2) <- names(muvar2)
    fe1 <- allFitnessEffects(noIntGenes = ni2)
    no <- 5e5
    reps <- 1000
    bb <- oncoSimulPop(reps,
                       fe1, mu = muvar2, onlyCancer = FALSE,
                       initSize = no,
                       finalTime = 4
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



## Use other models above

## For complex:
##  a) have it not be equal with large fitness diffs
##  b) have it be equal with tiny fitness diffs



test_that("Complex fitness specification, tiny s diffs", {

    p4 <- data.frame(parent = c(rep("Root", 4), "A", "B", "D", "E", "C", "F"),
                     child = c("A", "B", "D", "E", "C", "C", "F", "F", "G", "G"),
                     s = c(0.001, 0.002, 0.003, 0.004, 0.001, 0.001, 0.002, 0.002, 0.003, 0.003),
                     sh = c(rep(0, 4), c(-.00009, -.00009), c(-.000095, -.000095), c(-.000099, -.000099)),
                     typeDep = c(rep("--", 4), 
                                 "XMPN", "XMPN", "MN", "MN", "SM", "SM"))
    
    oe <- c("C > F" = -0.001, "H > I" = 0.0012)
    sm <- c("I:J"  = -.001)
    sv <- c("-K:M" = -.0005, "K:-M" = -.0005)
    epist <- c(sm, sv)
    modules <- c("Root" = "Root", "A" = "a1",
                 "B" = "b1, b2", "C" = "c1",
                 "D" = "d1, d2", "E" = "e1",
                 "F" = "f1, f2", "G" = "g1",
                 "H" = "h1, h2", "I" = "i1",
                 "J" = "j1, j2", "K" = "k1, k2", "M" = "m1")
    noint <- runif(5, min = 0.00051, max = 0.001)
    names(noint) <- paste0("n", 1:5)
    fea <- allFitnessEffects(rT = p4, epistasis = epist, orderEffects = oe,
                             noIntGenes = noint, geneToModule = modules)
    nfea <- OncoSimulR:::allNamedGenes(fea)$Gene
    ## muvar <- runif(length(nfea), min = 1e-7, max = 1e-3) ## too tiny
    ## diffs sometimes for order comp
    muvar <- sample(seq(from = 1e-7, to = 1e-3, length.out = length(nfea)))
    names(muvar) <- nfea
    no <- 1e4
    reps <- 10000
    bb <- oncoSimulPop(reps,
                       fea, mu = muvar,
                       onlyCancer = FALSE,
                       initSize = no,
                       finalTime = 1
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
                       finalTime = 1
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
