## ##     A way to check is to see the output from the C++ code with the
## ##     verbosity option.

## ni <- rep(0.4, 20)
## names(ni) <- c("a", "b", "c", "d", paste0("n", 1:16))
## fe <- allFitnessEffects(noIntGenes = ni)
## set.seed(5) 
## oncoSimulIndiv(fe, finalTime =30,
##                mutationPropGrowth = TRUE,
##                initSize = 1e4,
##                mu = 1e-06,
##                verbosity = 6,
##                onlyCancer = FALSE)

## ###### Iteration 30.
## ## mutation
## ## child
## 1.4 * 1e-06 * 19

## ni <- rep(0.4, 20)
## names(ni) <- c("a", "b", "c", "d", paste0("n", 1:16))
## fe <- allFitnessEffects(noIntGenes = ni)
## set.seed(25) 
## oncoSimulIndiv(fe, finalTime =40,
##                mutationPropGrowth = TRUE,
##                initSize = 1e4,
##                mu = 1e-06,
##                verbosity = 6,
##                onlyCancer = FALSE)
## ## Iteration 48.
## ## Birth of child:
## 1.4 * 1.4
## ## Mutation of child
## 1.96 * 1e-06 * 18

RNGkind("L'Ecuyer-CMRG") ## for the mclapplies

## The tests below can occasionally fail (but that probability decreases
## as we increase number of pops), as they should.

## I fix the seed for now.
test_that("Ordering of number of clones with mutpropgrowth", {
    set.seed(1)
    pops <- 50
    lni <- 200
    no <- 5e3
    ni <- c(5, 2, rep(0, lni))
    names(ni) <- c("a", "b", paste0("n", seq.int(lni)))
    fe <- allFitnessEffects(noIntGenes = ni)
    nca <- oncoSimulPop(pops, fe, finalTime =1,
                        mutationPropGrowth = TRUE,
                        initSize = no,
                        initMutant = "a",
                        onlyCancer = FALSE, seed = NULL)
    ncb <- oncoSimulPop(pops, fe, finalTime =1,
                        mutationPropGrowth = TRUE,
                        initSize = no,
                        initMutant = "b",
                        onlyCancer = FALSE, seed = NULL)
    nca2 <- oncoSimulPop(pops, fe, finalTime =1,
                        mutationPropGrowth = FALSE,
                        initSize = no,
                        initMutant = "a",
                        onlyCancer = FALSE, seed = NULL)
    ncb2 <- oncoSimulPop(pops, fe, finalTime =1,
                        mutationPropGrowth = FALSE,
                        initSize = no,
                        initMutant = "b",
                        onlyCancer = FALSE, seed = NULL)
    ## I once saw a weird thing
    expect_true(var(summary(nca)$NumClones) > 1e-4)
    expect_true(var(summary(ncb)$NumClones) > 1e-4)
    expect_true(var(summary(nca2)$NumClones) > 1e-4)
    expect_true(var(summary(ncb2)$NumClones) > 1e-4)
    ## The real comparison
    expect_true( median(summary(nca)$NumClones) >
                 median(summary(nca2)$NumClones))
    expect_true( median(summary(nca)$NumClones) >
                 median(summary(ncb)$NumClones))
    expect_true( median(summary(ncb)$NumClones) >
                 median(summary(ncb2)$NumClones))
})

test_that("Ordering of number of clones with mutpropgrowth, McFL", {
    set.seed(2)
    pops <- 50
    lni <- 200
    no <- 5e3
    ni <- c(5, 2, rep(0, lni))
    names(ni) <- c("a", "b", paste0("n", seq.int(lni)))
    fe <- allFitnessEffects(noIntGenes = ni)
    nca <- oncoSimulPop(pops, fe, finalTime =1,
                        mutationPropGrowth = TRUE,
                        initSize = no, model = "McFL",
                        initMutant = "a",
                        onlyCancer = FALSE, seed = NULL)
    ncb <- oncoSimulPop(pops, fe, finalTime =1,
                        mutationPropGrowth = TRUE,
                        initSize = no, model = "McFL",
                        initMutant = "b",
                        onlyCancer = FALSE, seed = NULL)
    nca2 <- oncoSimulPop(pops, fe, finalTime =1,
                        mutationPropGrowth = FALSE,
                        initSize = no, model = "McFL",
                        initMutant = "a",
                        onlyCancer = FALSE, seed = NULL)
    ncb2 <- oncoSimulPop(pops, fe, finalTime =1,
                        mutationPropGrowth = FALSE,
                        initSize = no, model = "McFL",
                        initMutant = "b",
                        onlyCancer = FALSE, seed = NULL)
    ## I once saw a weird thing
    expect_true(var(summary(nca)$NumClones) > 1e-4)
    expect_true(var(summary(ncb)$NumClones) > 1e-4)
    expect_true(var(summary(nca2)$NumClones) > 1e-4)
    expect_true(var(summary(ncb2)$NumClones) > 1e-4)
    ## The real comparison
    expect_true( median(summary(nca)$NumClones) >
                 median(summary(nca2)$NumClones))
    expect_true( median(summary(nca)$NumClones) >
                 median(summary(ncb)$NumClones))
    expect_true( median(summary(ncb)$NumClones) >
                 median(summary(ncb2)$NumClones))
})
