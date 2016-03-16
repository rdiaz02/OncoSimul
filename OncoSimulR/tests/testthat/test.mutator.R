
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



## test with var mut rate,
## run all tests
## create new tests

## oncosimulPop
## docs:
##    - help
##  -fignete

## Fitness of 0, but mutator effects.

## Modules same and different from fitness effects.


## check fail if mutator and fitness not subseted in calls that use both.

## check fail if mutator and fitness not both given in the FitAndMut functions.


## test changes in genotype accumulation with mutator effects (1, >>1,
## <<1) and number of mutator genes


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


fm6 <- allMutatorEffects(noIntGenes = c("a" = 1e2,
                                        "b" = 1,
                                        "d" = 1,
                                        "c" = 1e2))
oncoSimulIndiv(fe, muEF = fm6, finalTime = 100, initSize = 1e5, onlyCancer = FALSE)


fe <- allFitnessEffects(noIntGenes = c("a" = 0.12,
                                       "b" = 0.14,
                                       "c" = 0.16,
                                       "d" = 0.11))
fm6 <- allMutatorEffects(noIntGenes = c("a" = 1e2,
                                        "b" = 1,
                                        "c" = 1,
                                        "d" = 1e2))
oncoSimulIndiv(fe, muEF = fm6, finalTime = 5,
               mutationPropGrowth = FALSE,
               initSize = 1e5,
               onlyCancer = FALSE, verbosity = 6)




oncoSimulIndiv(fe, muEF = fm, finalTime = 100, initSize = 1e5, onlyCancer = FALSE)
oncoSimulIndiv(fe, muEF = fm2, finalTime = 100, initSize = 1e5, onlyCancer = FALSE)
oncoSimulIndiv(fe, muEF = fm3, finalTime = 100, initSize = 1e5, onlyCancer = FALSE)
oncoSimulIndiv(fe, muEF = fm4, finalTime = 100, initSize = 1e5, onlyCancer = FALSE)

