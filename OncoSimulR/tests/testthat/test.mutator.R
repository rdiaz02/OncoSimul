medianNClones <- function(x) {
    median(summary(x)$NumClones)
}


test_that("eval fitness and mut OK", {
    fe <- allFitnessEffects(epistasis = c("a : b" = 0.3,
                                          "b : c" = 0.5),
                            noIntGenes = c("e" = 0.1))
    fm <- OncoSimulR:::allMutatorEffects(noIntGenes = c("a" = 10,
                                                        "c" = 5))
    expect_output(ou <- evalGenotypeFitAndMut("a", fe, fm),
                  "Genotype", fixed = TRUE)
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
    expect_output(oncoSimulIndiv(fe, muEF = fm),
                  "Individual OncoSimul trajectory",
                  fixed = TRUE)
    expect_output(oncoSimulIndiv(fe),
                  "Individual OncoSimul trajectory",
                  fixed = TRUE)
})





test_that("eval mut genotypes", {
    fe <- allFitnessEffects(epistasis = c("a : b" = 0.3,
                                          "b : c" = 0.5),
                            noIntGenes = c("e" = 0.1))
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
                 "genotype contains NAs or genes not in fitnessEffects",
                 fixed = TRUE)
})

test_that("we evaluate the WT", {
    ## Is fitness of wildtype always 0? Really? Evaluate it.
    ## It is: see evalGenotypeFitness
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
                            noIntGenes = c("e" = 0.1))
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


test_that("Relative ordering of number of clones with mutator effects", {
    pops <- 5
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
                        initSize = 1e6,
                        onlyCancer = FALSE)
    fe <- allFitnessEffects(noIntGenes = c("a" = 0.12,
                                           "b" = 0.14,
                                           "c" = 0.16,
                                           "d" = 0.11))
    fm8 <- allMutatorEffects(noIntGenes = c("a" = 1,
                                            "b" = 1,
                                            "c" = 1,
                                            "d" = 1))
    nc2 <- oncoSimulPop(pops, fe, muEF = fm8, finalTime =250,
                        mutationPropGrowth = FALSE,
                        initSize = 1e6,
                        onlyCancer = FALSE)
    fe <- allFitnessEffects(noIntGenes = c("a" = 0.12,
                                           "b" = 0.14,
                                           "c" = 0.16,
                                           "d" = 0.11))
    fm7 <- allMutatorEffects(noIntGenes = c("a" = 1e-6,
                                            "b" = 1e-6,
                                            "c" = 1e-6,
                                            "d" = 1e-6))
    nc3 <- oncoSimulPop(pops, fe, muEF = fm7, finalTime =250,
                        mutationPropGrowth = FALSE,
                        initSize = 1e6,
                        onlyCancer = FALSE)
    expect_true(median(summary(nc1)$NumClones) > median(summary(nc2)$NumClones))
    expect_true(median(summary(nc2)$NumClones) > median(summary(nc3)$NumClones))
})


test_that("McFL: Relative ordering of number of clones with mutator effects", {
    pops <- 5
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
                        initSize  = 1e6, model = "McFL",
                        onlyCancer = FALSE)
    fe <- allFitnessEffects(noIntGenes = c("a" = 0.12,
                                           "b" = 0.14,
                                           "c" = 0.16,
                                           "d" = 0.11))
    fm8 <- allMutatorEffects(noIntGenes = c("a" = 1,
                                            "b" = 1,
                                            "c" = 1,
                                            "d" = 1))
    nc2 <- oncoSimulPop(pops, fe, muEF = fm8, finalTime =250,
                        mutationPropGrowth = FALSE,
                        initSize  = 1e6, model = "McFL",
                        onlyCancer = FALSE)
    fe <- allFitnessEffects(noIntGenes = c("a" = 0.12,
                                           "b" = 0.14,
                                           "c" = 0.16,
                                           "d" = 0.11))
    fm7 <- allMutatorEffects(noIntGenes = c("a" = 1e-6,
                                            "b" = 1e-6,
                                            "c" = 1e-6,
                                            "d" = 1e-6))
    nc3 <- oncoSimulPop(pops, fe, muEF = fm7, finalTime =250,
                        mutationPropGrowth = FALSE,
                        initSize  = 1e6, model = "McFL",
                        onlyCancer = FALSE)
    expect_true(median(summary(nc1)$NumClones) > median(summary(nc2)$NumClones))
    expect_true(median(summary(nc2)$NumClones) > median(summary(nc3)$NumClones))
})



test_that("Relative ordering of number of clones with init mutant of mutator effects", {
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
                        onlyCancer = FALSE)
    ncb <- oncoSimulPop(pops, fe, muEF = fm6, finalTime =50,
                        mutationPropGrowth = FALSE,
                        initSize = 1e4,
                        initMutant = "b",
                        onlyCancer = FALSE)
    ncc <- oncoSimulPop(pops, fe, muEF = fm6, finalTime =50,
                        mutationPropGrowth = FALSE,
                        initSize = 1e4,
                        initMutant = "c",
                        onlyCancer = FALSE)
    ncd <- oncoSimulPop(pops, fe, muEF = fm6, finalTime =50,
                        mutationPropGrowth = FALSE,
                        initSize = 1e4,
                        initMutant = "d",
                        onlyCancer = FALSE)
    expect_true( median(summary(nca)$NumClones) <
                 median(summary(ncb)$NumClones))
    expect_true(median(summary(ncb)$NumClones) <
                median(summary(ncc)$NumClones) )
    expect_true( median(summary(ncc)$NumClones) <
                 median(summary(ncd)$NumClones) )
})



test_that("Relative ordering of number of clones with init mutant of mutator effects and s = 0", {
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
                        onlyCancer = FALSE)
    ncb <- oncoSimulPop(pops, fe, muEF = fm6, finalTime =50,
                        mutationPropGrowth = FALSE,
                        initSize = 1e4,
                        initMutant = "b",
                        onlyCancer = FALSE)
    ncc <- oncoSimulPop(pops, fe, muEF = fm6, finalTime =50,
                        mutationPropGrowth = FALSE,
                        initSize = 1e4,
                        initMutant = "c",
                        onlyCancer = FALSE)
    ncd <- oncoSimulPop(pops, fe, muEF = fm6, finalTime =50,
                        mutationPropGrowth = FALSE,
                        initSize = 1e4,
                        initMutant = "d",
                        onlyCancer = FALSE)
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
})



test_that("Relative ordering of number of clones with mut prop growth and init and scrambled names", {
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
                        onlyCancer = FALSE)
    mnpg <- oncoSimulPop(pops, fe, muEF = fm1,
                         finalTime = ft,
                         mutationPropGrowth = FALSE,
                         initSize = no,
                         initMutant = "thisistheagene",
                         onlyCancer = FALSE)
    pg <- oncoSimulPop(pops, fe, 
                       finalTime = ft,
                       mutationPropGrowth = TRUE,
                       initSize = no,
                       initMutant = "thisistheagene",
                       onlyCancer = FALSE)
    npg <- oncoSimulPop(pops, fe, 
                        finalTime = ft,
                        mutationPropGrowth = FALSE,
                        initSize = no,
                        initMutant = "thisistheagene",
                        onlyCancer = FALSE)
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
})


test_that("McFL: Relative ordering of number of clones with mut prop growth and init and scrambled names", {
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
                        onlyCancer = FALSE)
    mnpg <- oncoSimulPop(pops, fe, muEF = fm1,
                         finalTime = ft,
                         mutationPropGrowth = FALSE,
                         initSize = no, model = "McFL",
                         initMutant = "thisistheagene",
                         onlyCancer = FALSE)
    pg <- oncoSimulPop(pops, fe, 
                       finalTime = ft,
                       mutationPropGrowth = TRUE,
                       initSize = no, model = "McFL",
                       initMutant = "thisistheagene",
                       onlyCancer = FALSE)
    npg <- oncoSimulPop(pops, fe, 
                        finalTime = ft,
                        mutationPropGrowth = FALSE,
                        initSize = no, model = "McFL",
                        initMutant = "thisistheagene",
                        onlyCancer = FALSE)
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
})


test_that("Num clones: Mutator and var mut rates and init and really scrambled names", {
    pops <- 10
    ft <- 1
    lni <- 200
    no <- 5e5
    ni <- c(0, 0, 0, rep(0, lni))
    ## scramble around names
    names(ni) <- c("hereisoneagene",
                   "bereisisabgene",
                   "nnhsisthecgene",
                   replicate(lni,
                             paste(sample(letters, 12), collapse = "")))
    ni <- ni[order(names(ni))]
    fe <- allFitnessEffects(noIntGenes = ni)
    mutator1 <- mutator2 <- rep(1, lni + 3)
    pg1 <- pg2 <- rep(1e-7, lni + 3)
    ## scramble names of mutator and per-gene too
    names(mutator1) <- sample(names(ni))
    names(mutator2) <- sample(names(ni))
    names(pg1) <- sample(names(ni))
    names(pg2) <- sample(names(ni))
    mutator1["hereisoneagene"] <- 100
    mutator2["hereisoneagene"] <- 1
    mutator1["bereisisabgene"] <- 1
    mutator2["bereisisabgene"] <- 0.01
    mutator1["nnhsisthecgene"] <- 0.01
    mutator2["nnhsisthecgene"] <- 100
    m1 <- allMutatorEffects(noIntGenes = mutator1)
    m2 <- allMutatorEffects(noIntGenes = mutator2)
    pg1["hereisoneagene"] <- 1e-3
    pg2["hereisoneagene"] <- 1e-14
    pg1["bereisisabgene"] <- 1e-7
    pg2["bereisisabgene"] <- 1e-7
    pg1["nnhsisthecgene"] <- 1e-14
    pg2["nnhsisthecgene"] <- 1e-3
    m1.pg1.a <- oncoSimulPop(pops,
                           fe,
                           mu = pg1,
                           muEF = m1,
                           finalTime = ft,
                           mutationPropGrowth = FALSE,
                           initSize = no,
                           initMutant = "hereisoneagene",
                           onlyCancer = FALSE)
    m1.pg2.a <- oncoSimulPop(pops,
                           fe,
                           mu = pg2,
                           muEF = m1,
                           finalTime = ft,
                           mutationPropGrowth = FALSE,
                           initSize = no,
                           initMutant = "hereisoneagene",
                           onlyCancer = FALSE)
    m1.pg1.c <- oncoSimulPop(pops,
                           fe,
                           mu = pg1,
                           muEF = m1,
                           finalTime = ft,
                           mutationPropGrowth = FALSE,
                           initSize = no,
                           initMutant = "nnhsisthecgene",
                           onlyCancer = FALSE)
    m1.pg2.c <- oncoSimulPop(pops,
                           fe,
                           mu = pg2,
                           muEF = m1,
                           finalTime = ft,
                           mutationPropGrowth = FALSE,
                           initSize = no,
                           initMutant = "nnhsisthecgene",
                           onlyCancer = FALSE)
    m1.pg1.b <- oncoSimulPop(pops,
                           fe,
                           mu = pg1,
                           muEF = m1,
                           finalTime = ft,
                           mutationPropGrowth = FALSE,
                           initSize = no,
                           initMutant = "bereisisabgene",
                           onlyCancer = FALSE)
    expect_true(medianNClones(m1.pg1.a) < medianNClones(m1.pg2.a))
    expect_true(medianNClones(m1.pg1.c) > medianNClones(m1.pg2.c))
    expect_true(medianNClones(m1.pg1.a) > medianNClones(m1.pg1.c))
    expect_true(medianNClones(m1.pg2.a) > medianNClones(m1.pg2.c))
    expect_true(medianNClones(m1.pg1.a) > medianNClones(m1.pg1.b))
    m2.pg1.a <- oncoSimulPop(pops,
                           fe,
                           mu = pg1,
                           muEF = m2,
                           finalTime = ft,
                           mutationPropGrowth = FALSE,
                           initSize = no,
                           initMutant = "hereisoneagene",
                           onlyCancer = FALSE)
    m2.pg2.a <- oncoSimulPop(pops,
                           fe,
                           mu = pg2,
                           muEF = m2,
                           finalTime = ft,
                           mutationPropGrowth = FALSE,
                           initSize = no,
                           initMutant = "hereisoneagene",
                           onlyCancer = FALSE)
    m2.pg1.c <- oncoSimulPop(pops,
                           fe,
                           mu = pg1,
                           muEF = m2,
                           finalTime = ft,
                           mutationPropGrowth = FALSE,
                           initSize = no,
                           initMutant = "nnhsisthecgene",
                           onlyCancer = FALSE)
    m2.pg2.c <- oncoSimulPop(pops,
                           fe,
                           mu = pg2,
                           muEF = m2,
                           finalTime = ft,
                           mutationPropGrowth = FALSE,
                           initSize = no,
                           initMutant = "nnhsisthecgene",
                           onlyCancer = FALSE)
    expect_true(medianNClones(m2.pg1.a) < medianNClones(m2.pg2.a))
    expect_true(medianNClones(m2.pg1.c) > medianNClones(m2.pg2.c))
    expect_true(medianNClones(m2.pg1.a) < medianNClones(m2.pg1.c))
    expect_true(medianNClones(m2.pg2.a) < medianNClones(m2.pg2.c))
    expect_true(medianNClones(m1.pg1.c) < medianNClones(m2.pg1.c))
    expect_true(medianNClones(m1.pg1.c) < medianNClones(m2.pg2.c))
    expect_true(medianNClones(m2.pg1.a) < medianNClones(m1.pg1.a))
    expect_true(medianNClones(m2.pg1.a) < medianNClones(m1.pg2.a))
    expect_true(medianNClones(m2.pg2.a) < medianNClones(m2.pg1.c))
    expect_true(medianNClones(m1.pg2.a) > medianNClones(m1.pg1.c))
})


test_that("McFL: Num clones: Mutator and var mut rates and init and really scrambled names", {
    pops <- 10
    ft <- 1
    lni <- 200
    no <- 2e5
    ni <- c(0, 0, 0, rep(0, lni))
    ## scramble around names
    names(ni) <- c("hereisoneagene",
                   "bereisisabgene",
                   "nnhsisthecgene",
                   replicate(lni,
                             paste(sample(letters, 12), collapse = "")))
    ni <- ni[order(names(ni))]
    fe <- allFitnessEffects(noIntGenes = ni)
    mutator1 <- mutator2 <- rep(1, lni + 3)
    pg1 <- pg2 <- rep(1e-7, lni + 3)
    ## scramble names of mutator and per-gene too
    names(mutator1) <- sample(names(ni))
    names(mutator2) <- sample(names(ni))
    names(pg1) <- sample(names(ni))
    names(pg2) <- sample(names(ni))
    mutator1["hereisoneagene"] <- 100
    mutator2["hereisoneagene"] <- 1
    mutator1["bereisisabgene"] <- 1
    mutator2["bereisisabgene"] <- 0.01
    mutator1["nnhsisthecgene"] <- 0.01
    mutator2["nnhsisthecgene"] <- 100
    m1 <- allMutatorEffects(noIntGenes = mutator1)
    m2 <- allMutatorEffects(noIntGenes = mutator2)
    pg1["hereisoneagene"] <- 1e-3
    pg2["hereisoneagene"] <- 1e-14
    pg1["bereisisabgene"] <- 1e-7
    pg2["bereisisabgene"] <- 1e-7
    pg1["nnhsisthecgene"] <- 1e-14
    pg2["nnhsisthecgene"] <- 1e-3
    m1.pg1.a <- oncoSimulPop(pops,
                           fe,
                           mu = pg1,
                           muEF = m1,
                           finalTime = ft, model = "McFL",
                           mutationPropGrowth = FALSE,
                           initSize = no,
                           initMutant = "hereisoneagene",
                           onlyCancer = FALSE)
    m1.pg2.a <- oncoSimulPop(pops,
                           fe,
                           mu = pg2,
                           muEF = m1,
                           finalTime = ft, model = "McFL",
                           mutationPropGrowth = FALSE,
                           initSize = no,
                           initMutant = "hereisoneagene",
                           onlyCancer = FALSE)
    m1.pg1.c <- oncoSimulPop(pops,
                           fe,
                           mu = pg1,
                           muEF = m1,
                           finalTime = ft, model = "McFL",
                           mutationPropGrowth = FALSE,
                           initSize = no,
                           initMutant = "nnhsisthecgene",
                           onlyCancer = FALSE)
    m1.pg2.c <- oncoSimulPop(pops,
                           fe,
                           mu = pg2,
                           muEF = m1,
                           finalTime = ft, model = "McFL",
                           mutationPropGrowth = FALSE,
                           initSize = no,
                           initMutant = "nnhsisthecgene",
                           onlyCancer = FALSE)
    m1.pg1.b <- oncoSimulPop(pops,
                           fe,
                           mu = pg1,
                           muEF = m1,
                           finalTime = ft, model = "McFL",
                           mutationPropGrowth = FALSE,
                           initSize = no,
                           initMutant = "bereisisabgene",
                           onlyCancer = FALSE)
    expect_true(medianNClones(m1.pg1.a) < medianNClones(m1.pg2.a))
    expect_true(medianNClones(m1.pg1.c) > medianNClones(m1.pg2.c))
    expect_true(medianNClones(m1.pg1.a) > medianNClones(m1.pg1.c))
    expect_true(medianNClones(m1.pg2.a) > medianNClones(m1.pg2.c))
    expect_true(medianNClones(m1.pg1.a) > medianNClones(m1.pg1.b))
    m2.pg1.a <- oncoSimulPop(pops,
                           fe,
                           mu = pg1,
                           muEF = m2,
                           finalTime = ft, model = "McFL",
                           mutationPropGrowth = FALSE,
                           initSize = no,
                           initMutant = "hereisoneagene",
                           onlyCancer = FALSE)
    m2.pg2.a <- oncoSimulPop(pops,
                           fe,
                           mu = pg2,
                           muEF = m2,
                           finalTime = ft, model = "McFL",
                           mutationPropGrowth = FALSE,
                           initSize = no,
                           initMutant = "hereisoneagene",
                           onlyCancer = FALSE)
    m2.pg1.c <- oncoSimulPop(pops,
                           fe,
                           mu = pg1,
                           muEF = m2,
                           finalTime = ft, model = "McFL",
                           mutationPropGrowth = FALSE,
                           initSize = no,
                           initMutant = "nnhsisthecgene",
                           onlyCancer = FALSE)
    m2.pg2.c <- oncoSimulPop(pops,
                           fe,
                           mu = pg2,
                           muEF = m2,
                           finalTime = ft, model = "McFL",
                           mutationPropGrowth = FALSE,
                           initSize = no,
                           initMutant = "nnhsisthecgene",
                           onlyCancer = FALSE)
    expect_true(medianNClones(m2.pg1.a) < medianNClones(m2.pg2.a))
    expect_true(medianNClones(m2.pg1.c) > medianNClones(m2.pg2.c))
    expect_true(medianNClones(m2.pg1.a) < medianNClones(m2.pg1.c))
    expect_true(medianNClones(m2.pg2.a) < medianNClones(m2.pg2.c))
    expect_true(medianNClones(m1.pg1.c) < medianNClones(m2.pg1.c))
    expect_true(medianNClones(m1.pg1.c) < medianNClones(m2.pg2.c))
    expect_true(medianNClones(m2.pg1.a) < medianNClones(m1.pg1.a))
    expect_true(medianNClones(m2.pg1.a) < medianNClones(m1.pg2.a))
    expect_true(medianNClones(m2.pg2.a) < medianNClones(m2.pg1.c))
    expect_true(medianNClones(m1.pg2.a) > medianNClones(m1.pg1.c))
})




test_that("Num clones: Mutator and var mut rates and init and really scrambled names, mutPropgrowth", {
    pops <- 10
    ft <- 1
    lni <- 200
    no <- 5e5
    ni <- c(0, 0, 0, rep(0, lni))
    ## scramble around names
    names(ni) <- c("hereisoneagene",
                   "bereisisabgene",
                   "nnhsisthecgene",
                   replicate(lni,
                             paste(sample(letters, 12), collapse = "")))
    ni <- ni[order(names(ni))]
    fe <- allFitnessEffects(noIntGenes = ni)
    mutator1 <- mutator2 <- rep(1, lni + 3)
    pg1 <- pg2 <- rep(1e-7, lni + 3)
    ## scramble names of mutator and per-gene too
    names(mutator1) <- sample(names(ni))
    names(mutator2) <- sample(names(ni))
    names(pg1) <- sample(names(ni))
    names(pg2) <- sample(names(ni))
    mutator1["hereisoneagene"] <- 100
    mutator2["hereisoneagene"] <- 1
    mutator1["bereisisabgene"] <- 1
    mutator2["bereisisabgene"] <- 0.01
    mutator1["nnhsisthecgene"] <- 0.01
    mutator2["nnhsisthecgene"] <- 100
    m1 <- allMutatorEffects(noIntGenes = mutator1)
    m2 <- allMutatorEffects(noIntGenes = mutator2)
    pg1["hereisoneagene"] <- 1e-3
    pg2["hereisoneagene"] <- 1e-14
    pg1["bereisisabgene"] <- 1e-7
    pg2["bereisisabgene"] <- 1e-7
    pg1["nnhsisthecgene"] <- 1e-14
    pg2["nnhsisthecgene"] <- 1e-3
    m1.pg1.a <- oncoSimulPop(pops,
                           fe,
                           mu = pg1,
                           muEF = m1,
                           finalTime = ft,
                           mutationPropGrowth = TRUE,
                           initSize = no,
                           initMutant = "hereisoneagene",
                           onlyCancer = FALSE)
    m1.pg2.a <- oncoSimulPop(pops,
                           fe,
                           mu = pg2,
                           muEF = m1,
                           finalTime = ft,
                           mutationPropGrowth = TRUE,
                           initSize = no,
                           initMutant = "hereisoneagene",
                           onlyCancer = FALSE)
    m1.pg1.c <- oncoSimulPop(pops,
                           fe,
                           mu = pg1,
                           muEF = m1,
                           finalTime = ft,
                           mutationPropGrowth = TRUE,
                           initSize = no,
                           initMutant = "nnhsisthecgene",
                           onlyCancer = FALSE)
    m1.pg2.c <- oncoSimulPop(pops,
                           fe,
                           mu = pg2,
                           muEF = m1,
                           finalTime = ft,
                           mutationPropGrowth = TRUE,
                           initSize = no,
                           initMutant = "nnhsisthecgene",
                           onlyCancer = FALSE)
    m1.pg1.b <- oncoSimulPop(pops,
                           fe,
                           mu = pg1,
                           muEF = m1,
                           finalTime = ft,
                           mutationPropGrowth = TRUE,
                           initSize = no,
                           initMutant = "bereisisabgene",
                           onlyCancer = FALSE)
    expect_true(medianNClones(m1.pg1.a) < medianNClones(m1.pg2.a))
    expect_true(medianNClones(m1.pg1.c) > medianNClones(m1.pg2.c))
    expect_true(medianNClones(m1.pg1.a) > medianNClones(m1.pg1.c))
    expect_true(medianNClones(m1.pg2.a) > medianNClones(m1.pg2.c))
    expect_true(medianNClones(m1.pg1.a) > medianNClones(m1.pg1.b))
    m2.pg1.a <- oncoSimulPop(pops,
                           fe,
                           mu = pg1,
                           muEF = m2,
                           finalTime = ft,
                           mutationPropGrowth = TRUE,
                           initSize = no,
                           initMutant = "hereisoneagene",
                           onlyCancer = FALSE)
    m2.pg2.a <- oncoSimulPop(pops,
                           fe,
                           mu = pg2,
                           muEF = m2,
                           finalTime = ft,
                           mutationPropGrowth = TRUE,
                           initSize = no,
                           initMutant = "hereisoneagene",
                           onlyCancer = FALSE)
    m2.pg1.c <- oncoSimulPop(pops,
                           fe,
                           mu = pg1,
                           muEF = m2,
                           finalTime = ft,
                           mutationPropGrowth = TRUE,
                           initSize = no,
                           initMutant = "nnhsisthecgene",
                           onlyCancer = FALSE)
    m2.pg2.c <- oncoSimulPop(pops,
                           fe,
                           mu = pg2,
                           muEF = m2,
                           finalTime = ft,
                           mutationPropGrowth = TRUE,
                           initSize = no,
                           initMutant = "nnhsisthecgene",
                           onlyCancer = FALSE)
    expect_true(medianNClones(m2.pg1.a) < medianNClones(m2.pg2.a))
    expect_true(medianNClones(m2.pg1.c) > medianNClones(m2.pg2.c))
    expect_true(medianNClones(m2.pg1.a) < medianNClones(m2.pg1.c))
    expect_true(medianNClones(m2.pg2.a) < medianNClones(m2.pg2.c))
    expect_true(medianNClones(m1.pg1.c) < medianNClones(m2.pg1.c))
    expect_true(medianNClones(m1.pg1.c) < medianNClones(m2.pg2.c))
    expect_true(medianNClones(m2.pg1.a) < medianNClones(m1.pg1.a))
    expect_true(medianNClones(m2.pg1.a) < medianNClones(m1.pg2.a))
    expect_true(medianNClones(m2.pg2.a) < medianNClones(m2.pg1.c))
    expect_true(medianNClones(m1.pg2.a) > medianNClones(m1.pg1.c))
})






test_that("McFL: Num clones: Mutator and var mut rates and init and really scrambled names, mutPropGrowt", {
    pops <- 10
    ft <- 1
    lni <- 200
    no <- 2e5
    ni <- c(0, 0, 0, rep(0, lni))
    ## scramble around names
    names(ni) <- c("hereisoneagene",
                   "bereisisabgene",
                   "nnhsisthecgene",
                   replicate(lni,
                             paste(sample(letters, 12), collapse = "")))
    ni <- ni[order(names(ni))]
    fe <- allFitnessEffects(noIntGenes = ni)
    mutator1 <- mutator2 <- rep(1, lni + 3)
    pg1 <- pg2 <- rep(1e-7, lni + 3)
    ## scramble names of mutator and per-gene too
    names(mutator1) <- sample(names(ni))
    names(mutator2) <- sample(names(ni))
    names(pg1) <- sample(names(ni))
    names(pg2) <- sample(names(ni))
    mutator1["hereisoneagene"] <- 100
    mutator2["hereisoneagene"] <- 1
    mutator1["bereisisabgene"] <- 1
    mutator2["bereisisabgene"] <- 0.01
    mutator1["nnhsisthecgene"] <- 0.01
    mutator2["nnhsisthecgene"] <- 100
    m1 <- allMutatorEffects(noIntGenes = mutator1)
    m2 <- allMutatorEffects(noIntGenes = mutator2)
    pg1["hereisoneagene"] <- 1e-3
    pg2["hereisoneagene"] <- 1e-14
    pg1["bereisisabgene"] <- 1e-7
    pg2["bereisisabgene"] <- 1e-7
    pg1["nnhsisthecgene"] <- 1e-14
    pg2["nnhsisthecgene"] <- 1e-3
    m1.pg1.a <- oncoSimulPop(pops,
                           fe,
                           mu = pg1,
                           muEF = m1,
                           finalTime = ft, model = "McFL",
                           mutationPropGrowth = TRUE,
                           initSize = no,
                           initMutant = "hereisoneagene",
                           onlyCancer = FALSE)
    m1.pg2.a <- oncoSimulPop(pops,
                           fe,
                           mu = pg2,
                           muEF = m1,
                           finalTime = ft, model = "McFL",
                           mutationPropGrowth = TRUE,
                           initSize = no,
                           initMutant = "hereisoneagene",
                           onlyCancer = FALSE)
    m1.pg1.c <- oncoSimulPop(pops,
                           fe,
                           mu = pg1,
                           muEF = m1,
                           finalTime = ft, model = "McFL",
                           mutationPropGrowth = TRUE,
                           initSize = no,
                           initMutant = "nnhsisthecgene",
                           onlyCancer = FALSE)
    m1.pg2.c <- oncoSimulPop(pops,
                           fe,
                           mu = pg2,
                           muEF = m1,
                           finalTime = ft, model = "McFL",
                           mutationPropGrowth = TRUE,
                           initSize = no,
                           initMutant = "nnhsisthecgene",
                           onlyCancer = FALSE)
    m1.pg1.b <- oncoSimulPop(pops,
                           fe,
                           mu = pg1,
                           muEF = m1,
                           finalTime = ft, model = "McFL",
                           mutationPropGrowth = TRUE,
                           initSize = no,
                           initMutant = "bereisisabgene",
                           onlyCancer = FALSE)
    expect_true(medianNClones(m1.pg1.a) < medianNClones(m1.pg2.a))
    expect_true(medianNClones(m1.pg1.c) > medianNClones(m1.pg2.c))
    expect_true(medianNClones(m1.pg1.a) > medianNClones(m1.pg1.c))
    expect_true(medianNClones(m1.pg2.a) > medianNClones(m1.pg2.c))
    expect_true(medianNClones(m1.pg1.a) > medianNClones(m1.pg1.b))
    m2.pg1.a <- oncoSimulPop(pops,
                           fe,
                           mu = pg1,
                           muEF = m2,
                           finalTime = ft, model = "McFL",
                           mutationPropGrowth = TRUE,
                           initSize = no,
                           initMutant = "hereisoneagene",
                           onlyCancer = FALSE)
    m2.pg2.a <- oncoSimulPop(pops,
                           fe,
                           mu = pg2,
                           muEF = m2,
                           finalTime = ft, model = "McFL",
                           mutationPropGrowth = TRUE,
                           initSize = no,
                           initMutant = "hereisoneagene",
                           onlyCancer = FALSE)
    m2.pg1.c <- oncoSimulPop(pops,
                           fe,
                           mu = pg1,
                           muEF = m2,
                           finalTime = ft, model = "McFL",
                           mutationPropGrowth = TRUE,
                           initSize = no,
                           initMutant = "nnhsisthecgene",
                           onlyCancer = FALSE)
    m2.pg2.c <- oncoSimulPop(pops,
                           fe,
                           mu = pg2,
                           muEF = m2,
                           finalTime = ft, model = "McFL",
                           mutationPropGrowth = TRUE,
                           initSize = no,
                           initMutant = "nnhsisthecgene",
                           onlyCancer = FALSE)
    expect_true(medianNClones(m2.pg1.a) < medianNClones(m2.pg2.a))
    expect_true(medianNClones(m2.pg1.c) > medianNClones(m2.pg2.c))
    expect_true(medianNClones(m2.pg1.a) < medianNClones(m2.pg1.c))
    expect_true(medianNClones(m2.pg2.a) < medianNClones(m2.pg2.c))
    expect_true(medianNClones(m1.pg1.c) < medianNClones(m2.pg1.c))
    expect_true(medianNClones(m1.pg1.c) < medianNClones(m2.pg2.c))
    expect_true(medianNClones(m2.pg1.a) < medianNClones(m1.pg1.a))
    expect_true(medianNClones(m2.pg1.a) < medianNClones(m1.pg2.a))
    expect_true(medianNClones(m2.pg2.a) < medianNClones(m2.pg1.c))
    expect_true(medianNClones(m1.pg2.a) > medianNClones(m1.pg1.c))
})





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
## I think very large per-gene mutations are not seen in expected props.
## humm... is this mutator effects?
## But remember we tested this in test.per-gene-mutation-rates?

## The following works fine
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
                       finalTime = 1,
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


test_that("Different freqs as they should be ordered and chisq, when s=0 and t = 1, and initMutant, many genotypes",{
    ## Occasionally, the c++ code can blow up as a large mu, etc

    lni <- 26  ## 16
    muvar2 <- rep(1e-5, lni)
    names(muvar2) <- c(replicate(lni,
                                 paste(sample(letters, 12), collapse = "")))
    muvar2[2] <- 2e-5 ## if this is, say, 5e-5, it fails. with
    names(muvar2)[3] <- "e"
    muvar2[2] <- 1e-3
    ni1 <- rep(0, lni)
    names(ni1) <- names(muvar2)
    fe1 <- allFitnessEffects(noIntGenes = ni1)
    no <- 1e5
    reps <- 500
    ## set.seed(6305)
    bb <- oncoSimulPop(reps,
                       fe1, mu = muvar2, onlyCancer = FALSE,
                       initSize = no,
                       finalTime = .5,
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

## I think this is what happens:

## if any mu is very large or any lni is very large, it can fail unless
## pops is very large. And having a large mutator effect is like having a
## very large mu.

## Play also with final time: make smaller than 1. It is birth and rate
## that compound processes and of course we have non-independent sampling
## (overdispersion)

test_that("Expect freq genotypes, mutator and var mut rates, small case", {

    ## Now play changing the mutation rate randomly
    
   ## Fails; is it the number of cells, of genotypes?
    pops <- 400
    ft <- .01
    lni <- 80 ## small, so not too many cells
    no <- 1e6
    ni <- c(0, 0, 0, rep(0, lni))
    ## scramble around names
    names(ni) <- c("hereisoneagene",
                   "bereisisabgene",
                   "nnhsisthecgene",
                   replicate(lni,
                             paste(sample(letters, 12), collapse = "")))
    ni <- ni[order(names(ni))]
    fe <- allFitnessEffects(noIntGenes = ni)
    mutator1 <- rep(1, lni + 3)
    ## pg1 <- rep(1e-5, lni + 3)
    pg1 <- runif(lni + 3, min = 1e-7, max = 1e-4)
    names(mutator1) <- sample(names(ni))
    names(pg1) <- sample(names(ni))
    mutator1["bereisisabgene"] <- 34    ## 53
    m1 <- allMutatorEffects(noIntGenes = mutator1)
    pg1["hereisoneagene"] <- 1e-3 ## 1e-3
    m1.pg1.b <- oncoSimulPop(pops,
                           fe,
                           mu = pg1,
                           muEF = m1,
                           finalTime = ft,
                           mutationPropGrowth = FALSE,
                           initSize = no,
                           initMutant ="bereisisabgene",
                           onlyCancer = FALSE)
    expect_true(sm("bereisisabgene", m1.pg1.b) == totalind(m1.pg1.b))
    enom("bereisisabgene", pg1)
    snom("bereisisabgene", m1.pg1.b)
    p.fail <- 1e-3
    expect_true(chisq.test(snom("bereisisabgene", m1.pg1.b),
                           p = pnom("bereisisabgene", pg1))$p.value > p.fail)











    
    ## OK
    lni <- 5  ## ok with 5
    pops <- 200
    ft <- 1
    no <- 5e5
    ni <- c(0, 0, 0, rep(0, lni))
    ## scramble around names
    names(ni) <- c("hereisoneagene",
                   "bereisisabgene",
                   "nnhsisthecgene",
                   replicate(lni,
                             paste(sample(letters, 12), collapse = "")))
    ni <- ni[order(names(ni))]
    fe <- allFitnessEffects(noIntGenes = ni)
    mutator1 <- rep(1, lni + 3)
    pg1 <- rep(1e-7, lni + 3)
    names(mutator1) <- sample(names(ni))
    names(pg1) <- sample(names(ni))
    mutator1["bereisisabgene"] <- 53
    m1 <- allMutatorEffects(noIntGenes = mutator1)
    pg1["hereisoneagene"] <- 1e-3 ## 1e-3
    m1.pg1.b <- oncoSimulPop(pops,
                           fe,
                           mu = pg1,
                           muEF = m1,
                           finalTime = ft,
                           mutationPropGrowth = FALSE,
                           initSize = no,
                           initMutant ="bereisisabgene",
                           onlyCancer = FALSE)
    expect_true(sm("bereisisabgene", m1.pg1.b) == totalind(m1.pg1.b))
    enom("bereisisabgene", pg1)
    snom("bereisisabgene", m1.pg1.b)
    p.fail <- 1e-3
    expect_true(chisq.test(snom("bereisisabgene", m1.pg1.b),
                           p = pnom("bereisisabgene", pg1))$p.value > p.fail)


    ## OK too
    pops <- 200
    ft <- 1
    lni <- 5 ## small, so not too many cells
    no <- 5e5
    ni <- c(0, 0, 0, rep(0, lni))
    ## scramble around names
    names(ni) <- c("hereisoneagene",
                   "bereisisabgene",
                   "nnhsisthecgene",
                   replicate(lni,
                             paste(sample(letters, 12), collapse = "")))
    ni <- ni[order(names(ni))]
    fe <- allFitnessEffects(noIntGenes = ni)
    mutator1 <- rep(1, lni + 3)
    pg1 <- rep(1e-7, lni + 3)
    names(mutator1) <- sample(names(ni))
    names(pg1) <- sample(names(ni))
    mutator1["bereisisabgene"] <- 53
    m1 <- allMutatorEffects(noIntGenes = mutator1)
    pg1["hereisoneagene"] <- 1e-3 ## 1e-3
    m1.pg1.b <- oncoSimulPop(pops,
                           fe,
                           mu = pg1,
                           muEF = m1,
                           finalTime = ft,
                           mutationPropGrowth = FALSE,
                           initSize = no,
                           initMutant ="bereisisabgene",
                           onlyCancer = FALSE)
    expect_true(sm("bereisisabgene", m1.pg1.b) == totalind(m1.pg1.b))
    enom("bereisisabgene", pg1)
    snom("bereisisabgene", m1.pg1.b)
    p.fail <- 1e-3
    expect_true(chisq.test(snom("bereisisabgene", m1.pg1.b),
                           p = pnom("bereisisabgene", pg1))$p.value > p.fail)



 
    

})






test_that("Expect freq genotypes: Mutator and var mut rates and init and really scrambled names", {

    pops <- 200
    ft <- 1
    lni <- 5 ## small, so not too many cells
    no <- 5e5
    ni <- c(0, 0, 0, rep(0, lni))
    ## scramble around names
    names(ni) <- c("hereisoneagene",
                   "bereisisabgene",
                   "nnhsisthecgene",
                   replicate(lni,
                             paste(sample(letters, 12), collapse = "")))
    ni <- ni[order(names(ni))]
    fe <- allFitnessEffects(noIntGenes = ni)
    mutator1 <- mutator2 <- rep(1, lni + 3)
    pg1 <- pg2 <- rep(1e-7, lni + 3)
    ## scramble names of mutator and per-gene too
    names(mutator1) <- sample(names(ni))
    names(mutator2) <- sample(names(ni))
    names(pg1) <- sample(names(ni))
    names(pg2) <- sample(names(ni))
    mutator1["hereisoneagene"] <- 100  ## 100
    mutator2["hereisoneagene"] <- 1
    mutator1["bereisisabgene"] <- 53
    mutator2["bereisisabgene"] <- 0.01
    mutator1["nnhsisthecgene"] <- 0.01
    mutator2["nnhsisthecgene"] <- 100
    m1 <- allMutatorEffects(noIntGenes = mutator1)
    m2 <- allMutatorEffects(noIntGenes = mutator2)
    pg1["hereisoneagene"] <- 5e-5 ## 1e-3
    pg2["hereisoneagene"] <- 1e-14
    pg1["bereisisabgene"] <- 1e-7
    pg2["bereisisabgene"] <- 1e-7
    pg1["nnhsisthecgene"] <- 1e-14
    pg2["nnhsisthecgene"] <- 1e-3
    ## something strange here
    ## I think very large per-gene mutations are not seen in expected props.
    ## humm...
    m1.pg1.b <- oncoSimulPop(pops,
                           fe,
                           mu = pg1,
                           muEF = m1,
                           finalTime = ft,
                           mutationPropGrowth = FALSE,
                           initSize = no,
                           initMutant ="bereisisabgene",
                           onlyCancer = FALSE)
    expect_true(sm("bereisisabgene", m1.pg1.b) == totalind(m1.pg1.b))
    enom("bereisisabgene", pg1)
    snom("bereisisabgene", m1.pg1.b)
    p.fail <- 1e-3
    expect_true(chisq.test(snom("bereisisabgene", m1.pg1.b),
                           p = pnom("bereisisabgene", pg1))$p.value > p.fail)

    
    
    m1.pg1.a <- oncoSimulPop(pops,
                           fe,
                           mu = pg1,
                           muEF = m1,
                           finalTime = ft,
                           mutationPropGrowth = FALSE,
                           initSize = no,
                           initMutant = "hereisoneagene",
                           onlyCancer = FALSE)
    expect_true(sm("hereisoneagene", m1.pg1.a) == totalind(m1.pg1.a))
    enom("hereisoneagene", pg1)
    snom("hereisoneagene", m1.pg1.a)
    p.fail <- 1e-3
    expect_true(chisq.test(snom("hereisoneagene", m1.pg1.a),
                           p = pnom("hereisoneagene", pg1))$p.value > p.fail)







    
    ## no change if different mutator rates, but of course need to increase
    ## pops and/or no as much smaller mut rates
    pops <- 500
    no <- 1e6
    m2.pg1.a <- oncoSimulPop(pops,
                           fe,
                           mu = pg1,
                           muEF = m2,
                           finalTime = ft,
                           mutationPropGrowth = FALSE,
                           initSize = no,
                           initMutant = "hereisoneagene",
                           onlyCancer = FALSE)

    expect_true(sm("hereisoneagene", m2.pg1.a) == totalind(m2.pg1.a))
    enom("hereisoneagene", pg1)
    snom("hereisoneagene", m2.pg1.a)
    p.fail <- 1e-3
    expect_true(chisq.test(snom("hereisoneagene", m2.pg1.a),
                           p = pnom("hereisoneagene", pg1))$p.value > p.fail)





    

    

    
    m1.pg2.a <- oncoSimulPop(pops,
                           fe,
                           mu = pg2,
                           muEF = m1,
                           finalTime = ft,
                           mutationPropGrowth = FALSE,
                           initSize = no,
                           initMutant = "hereisoneagene",
                           onlyCancer = FALSE)
    m1.pg1.c <- oncoSimulPop(pops,
                           fe,
                           mu = pg1,
                           muEF = m1,
                           finalTime = ft,
                           mutationPropGrowth = FALSE,
                           initSize = no,
                           initMutant = "nnhsisthecgene",
                           onlyCancer = FALSE)
    m1.pg2.c <- oncoSimulPop(pops,
                           fe,
                           mu = pg2,
                           muEF = m1,
                           finalTime = ft,
                           mutationPropGrowth = FALSE,
                           initSize = no,
                           initMutant = "nnhsisthecgene",
                           onlyCancer = FALSE)
    m1.pg1.b <- oncoSimulPop(pops,
                           fe,
                           mu = pg1,
                           muEF = m1,
                           finalTime = ft,
                           mutationPropGrowth = FALSE,
                           initSize = no,
                           initMutant = "bereisisabgene",
                           onlyCancer = FALSE)
    expect_true(medianNClones(m1.pg1.a) < medianNClones(m1.pg2.a))
    expect_true(medianNClones(m1.pg1.c) > medianNClones(m1.pg2.c))
    expect_true(medianNClones(m1.pg1.a) > medianNClones(m1.pg1.c))
    expect_true(medianNClones(m1.pg2.a) > medianNClones(m1.pg2.c))
    expect_true(medianNClones(m1.pg1.a) > medianNClones(m1.pg1.b))
    m2.pg1.a <- oncoSimulPop(pops,
                           fe,
                           mu = pg1,
                           muEF = m2,
                           finalTime = ft,
                           mutationPropGrowth = FALSE,
                           initSize = no,
                           initMutant = "hereisoneagene",
                           onlyCancer = FALSE)
    m2.pg2.a <- oncoSimulPop(pops,
                           fe,
                           mu = pg2,
                           muEF = m2,
                           finalTime = ft,
                           mutationPropGrowth = FALSE,
                           initSize = no,
                           initMutant = "hereisoneagene",
                           onlyCancer = FALSE)
    m2.pg1.c <- oncoSimulPop(pops,
                           fe,
                           mu = pg1,
                           muEF = m2,
                           finalTime = ft,
                           mutationPropGrowth = FALSE,
                           initSize = no,
                           initMutant = "nnhsisthecgene",
                           onlyCancer = FALSE)
    m2.pg2.c <- oncoSimulPop(pops,
                           fe,
                           mu = pg2,
                           muEF = m2,
                           finalTime = ft,
                           mutationPropGrowth = FALSE,
                           initSize = no,
                           initMutant = "nnhsisthecgene",
                           onlyCancer = FALSE)
    expect_true(medianNClones(m2.pg1.a) < medianNClones(m2.pg2.a))
    expect_true(medianNClones(m2.pg1.c) > medianNClones(m2.pg2.c))
    expect_true(medianNClones(m2.pg1.a) < medianNClones(m2.pg1.c))
    expect_true(medianNClones(m2.pg2.a) < medianNClones(m2.pg2.c))
    expect_true(medianNClones(m1.pg1.c) < medianNClones(m2.pg1.c))
    expect_true(medianNClones(m1.pg1.c) < medianNClones(m2.pg2.c))
    expect_true(medianNClones(m2.pg1.a) < medianNClones(m1.pg1.a))
    expect_true(medianNClones(m2.pg1.a) < medianNClones(m1.pg2.a))
    expect_true(medianNClones(m2.pg2.a) < medianNClones(m2.pg1.c))
    expect_true(medianNClones(m1.pg2.a) > medianNClones(m1.pg1.c))
})


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

## test with var mut rate,
## run all tests
## create new tests

## oncosimulPop
## docs:
##    - help
##  -fignete
##  - finish docs


## Modules same and different from fitness effects.


## check fail if mutator and fitness not both given in the FitAndMut functions.

## And use initMutant with different mutator effects






## ## Some checks of C++ code with mutator and per-gene mutation rates
## set.seed(2) 
## ft <- 600
## lni <- 10
## no <- 5e3
## ni <- c(0, 0, 0, rep(0, lni))
## ## scramble around names
## names(ni) <- c("hereisoneagene",
##                "bereisisabgene",
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
## mutator1["bereisisabgene"] <- 5
## mutator1["nnhsisthecgene"] <- 0.01
## m1 <- allMutatorEffects(noIntGenes = mutator1)
## pg1["hereisoneagene"] <- 1e-5
## pg1["bereisisabgene"] <- 1e-7
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
##                onlyCancer = FALSE)    
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




    









fe <- allFitnessEffects(noIntGenes = c("a" = 0.12,
                                       "c" = 0.14,
                                       "d" = 0.16,
                                       "e" = 0.11))
fm <- allMutatorEffects(noIntGenes = c("a" = 10,
                                       "c" = 5))
fm2 <- allMutatorEffects(noIntGenes = c("a" = .010,
                                           "c" = .05))
fm3 <- allMutatorEffects(noIntGenes = c("a" = 10,
                                        "b" = 3,
                                        "d" = 8,
                                        "c" = 5))
fm4 <- allMutatorEffects(noIntGenes = c("a" = .010,
                                        "b" = .03,
                                        "d" = .08,
                                        "c" = .05))




fm5 <- allMutatorEffects(noIntGenes = c("a" = 1e-6,
                                        "b" = 1e-6,
                                        "d" = 1e-6,
                                        "c" = 1e-6))
oncoSimulIndiv(fe, muEF = fm5, finalTime = 100, initSize = 1e5, onlyCancer = FALSE)


fe <- allFitnessEffects(noIntGenes = c("a" = 0.12,
                                       "c" = 0.14,
                                       "d" = 0.16,
                                       "e" = 0.11))
fm6 <- allMutatorEffects(noIntGenes = c("a" = 1e2,
                                        "b" = 1,
                                        "d" = 1,
                                        "c" = 1e2))
oncoSimulIndiv(fe, muEF = fm6, finalTime = 100, initSize = 1e5, onlyCancer = FALSE,
               verbosity = 6)



## Below: makes sense in terms of number of clones
fe <- allFitnessEffects(noIntGenes = c("a" = 0.12,
                                       "b" = 0.14,
                                       "c" = 0.16,
                                       "d" = 0.11))
fm6 <- allMutatorEffects(noIntGenes = c("a" = 10,
                                        "b" = 20,
                                        "c" = 30,
                                        "d" = 40))
oncoSimulIndiv(fe, muEF = fm6, finalTime =250,
               mutationPropGrowth = FALSE,
               initSize = 1e5,
               verbosity = 6,
               onlyCancer = FALSE)



fe <- allFitnessEffects(noIntGenes = c("a" = 0.12,
                                       "b" = 0.14,
                                       "c" = 0.16,
                                       "d" = 0.11))
fm6 <- allMutatorEffects(noIntGenes = c("a" = 1e-6,
                                        "b" = 1e-6,
                                        "c" = 1e-6,
                                        "d" = 1e-6))
oncoSimulIndiv(fe, muEF = fm6, finalTime =250,
               mutationPropGrowth = FALSE,
               initSize = 1e5,
               onlyCancer = FALSE)




oncoSimulIndiv(fe, muEF = fm, finalTime = 100, initSize = 1e5, onlyCancer = FALSE)
oncoSimulIndiv(fe, muEF = fm2, finalTime = 100, initSize = 1e5, onlyCancer = FALSE)
oncoSimulIndiv(fe, muEF = fm3, finalTime = 100, initSize = 1e5, onlyCancer = FALSE)
oncoSimulIndiv(fe, muEF = fm4, finalTime = 100, initSize = 1e5, onlyCancer = FALSE)




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
