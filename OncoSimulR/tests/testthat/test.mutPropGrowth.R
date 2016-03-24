RNGkind("L'Ecuyer-CMRG") ## for the mclapplies
## If crashes I want to see where: thus output seed.
## The tests below can occasionally fail (but that probability decreases
## as we increase number of pops), as they should.

cat("\n", date(), "\n") ## whole file takes about 6 seconds
test_that("Ordering of number of clones with mutpropgrowth", {
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n omp1: the seed is", pseed, "\n")
    pops <- 100
    lni <- 200
    no <- 5e3
    ni <- c(5, 2, rep(0, lni))
    names(ni) <- c("a", "b", paste0("n", seq.int(lni)))
    fe <- allFitnessEffects(noIntGenes = ni)
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n omp1a: the seed is", pseed, "\n")
    nca <- oncoSimulPop(pops, fe, finalTime =1,
                        mutationPropGrowth = TRUE,
                        initSize = no,
                        initMutant = "a",
                        onlyCancer = FALSE, seed = NULL, mc.cores = 2)
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n omp1b: the seed is", pseed, "\n")
    ncb <- oncoSimulPop(pops, fe, finalTime =1,
                        mutationPropGrowth = TRUE,
                        initSize = no,
                        initMutant = "b",
                        onlyCancer = FALSE, seed = NULL, mc.cores = 2)
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n omp1c: the seed is", pseed, "\n")
    nca2 <- oncoSimulPop(pops, fe, finalTime =1,
                        mutationPropGrowth = FALSE,
                        initSize = no,
                        initMutant = "a",
                        onlyCancer = FALSE, seed = NULL, mc.cores = 2)
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n omp1d: the seed is", pseed, "\n")
    ncb2 <- oncoSimulPop(pops, fe, finalTime =1,
                        mutationPropGrowth = FALSE,
                        initSize = no,
                        initMutant = "b",
                        onlyCancer = FALSE, seed = NULL, mc.cores = 2)
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
cat("\n", date(), "\n")

cat("\n", date(), "\n")
test_that("Ordering of number of clones with mutpropgrowth, McFL", {
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n omp2: the seed is", pseed, "\n")
    pops <- 100
    lni <- 200
    no <- 5e3
    ni <- c(5, 2, rep(0, lni))
    names(ni) <- c("a", "b", paste0("n", seq.int(lni)))
    fe <- allFitnessEffects(noIntGenes = ni)
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n omp2a: the seed is", pseed, "\n")
    nca <- oncoSimulPop(pops, fe, finalTime =1,
                        mutationPropGrowth = TRUE,
                        initSize = no, model = "McFL",
                        initMutant = "a",
                        onlyCancer = FALSE, seed = NULL, mc.cores = 2)
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n omp2b: the seed is", pseed, "\n")
    ncb <- oncoSimulPop(pops, fe, finalTime =1,
                        mutationPropGrowth = TRUE,
                        initSize = no, model = "McFL",
                        initMutant = "b",
                        onlyCancer = FALSE, seed = NULL, mc.cores = 2)
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n omp2c: the seed is", pseed, "\n")
    nca2 <- oncoSimulPop(pops, fe, finalTime =1,
                        mutationPropGrowth = FALSE,
                        initSize = no, model = "McFL",
                        initMutant = "a",
                        onlyCancer = FALSE, seed = NULL, mc.cores = 2)
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n omp2d: the seed is", pseed, "\n")
    ncb2 <- oncoSimulPop(pops, fe, finalTime =1,
                        mutationPropGrowth = FALSE,
                        initSize = no, model = "McFL",
                        initMutant = "b",
                        onlyCancer = FALSE, seed = NULL, mc.cores = 2)
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
cat("\n", date(), "\n")




