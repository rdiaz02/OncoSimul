
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



test_that("Relative ordering of number of clones with init mutant of mutator effects", {
    pops <- 5
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
    pops <- 5
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


test_that("Mutator and var mut rates and init and really scrambled names", {
    
    pops <- 5
    ft <- 1
    lni <- 200
    no <- 5e3
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
    pg1["hereisoneagene"] <- 1e-5
    pg2["hereisoneagene"] <- 1e-10
    pg1["bereisisabgene"] <- 1e-7
    pg2["bereisisabgene"] <- 1e-7
    pg1["nnhsisthecgene"] <- 1e-10
    pg2["nnhsisthecgene"] <- 1e-5

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

    


    
    m2.pg1.a <- oncoSimulPop(pops,
                           fe,
                           mu = pg1,
                           muEF = m2,
                           finalTime = ft,
                           mutationPropGrowth = FALSE,
                           initSize = no,
                           initMutant = "hereisoneagene",
                           onlyCancer = FALSE)

    


    
})
 



## breaks an assertion
set.seed(2) ## not with set.seed(1)
pops <- 5
ft <- 1
lni <- 200
no <- 5e3
ni <- c(0, 0, 0, rep(0, lni))
## scramble around names
names(ni) <- c("hereisoneagene",
               "bereisisabgene",
               "nnhsisthecgene",
               replicate(lni,
                         paste(sample(letters, 12), collapse = "")))
ni <- ni[order(names(ni))]
ni <- sample(ni)
fe <- allFitnessEffects(noIntGenes = ni)
mutator1 <- rep(1, lni + 3)
pg1 <- rep(1e-7, lni + 3)  ## what breaks is using 1e-7 instead of 1e-9
## scramble names of mutator and per-gene too
names(mutator1) <- sample(names(ni))
names(pg1) <- sample(names(ni))
mutator1["hereisoneagene"] <- 100
mutator1["bereisisabgene"] <- 1
mutator1["nnhsisthecgene"] <- 0.01
m1 <- allMutatorEffects(noIntGenes = mutator1)
pg1["hereisoneagene"] <- 1e-5
pg1["bereisisabgene"] <- 1e-7
pg1["nnhsisthecgene"] <- 1e-10
m1.pg1.a <- oncoSimulIndiv(fe,
                           mu = pg1,
                           muEF = m1,
                           finalTime = ft,
                           mutationPropGrowth = FALSE,
                           initSize = no,
                           initMutant = "hereisoneagene",
                           onlyCancer = FALSE,
                           verbosity = 6)

mu.nop <- function(p, mu) {
    nni <- names(ni)[p]
    pmu <- which(names(mu) %in%  nni)
    mu[-pmu]
}

100 * sum(mu.nop(c(45, 64), pg1))

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

## Modify the printing of the genotype?





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
