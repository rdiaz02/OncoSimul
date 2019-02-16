## The following tests is likely to crash in non-Linux and/or non-gcc
## machines (e.g., using clang), as the details depend on the random
## numbers generated in R and C++.

## In addition, changes in the C++ code can upset this. Say you change a
## detail about the RNG, or multiply by 2 mindummy, etc.
## So maybe I should remove these tests.

RNGkind("Mersenne-Twister")

date()
test_that("driverCounts: we get the right results", {
    ## Comparing against know results with something that used to give
    ## wrong ones

    ## I also check mue11$Genotypes
    ## to know what to expect
    
    set.seed(1)
    ni <- runif(10, min = -0.01, max = 0.1)
    names(ni) <- paste0("g", 2:11)
    fe31 <- allFitnessEffects(noIntGenes = ni,
                              drvNames = "g2") 
    set.seed(1)
    mue11 <- oncoSimulIndiv(fe31, 
                            mu = 1e-6,
                            initSize = 1e5,
                            model = "McFL",
                            detectionSize = 5e6,
                            finalTime = 5, 
                            onlyCancer = FALSE, detectionProb = NA)
    expect_true(mue11$MaxNumDrivers == 1)
    expect_true(mue11$MaxDriversLast == 0)

    set.seed(1)
    ni <- runif(10, min = -0.01, max = 0.1)
    names(ni) <- paste0("g", 2:11)
    fe31 <- allFitnessEffects(noIntGenes = ni,
                              drvNames = "g3") 
    set.seed(1)
    mue11 <- oncoSimulIndiv(fe31, 
                            mu = 1e-6,
                            initSize = 1e5,
                            model = "McFL",
                            detectionSize = 5e6,
                            finalTime = 5, 
                            onlyCancer = FALSE, detectionProb = NA)
    expect_true(mue11$MaxNumDrivers == 1)
    expect_true(mue11$MaxDriversLast == 1)

    set.seed(1)
    ni <- runif(10, min = -0.01, max = 0.1)
    names(ni) <- paste0("g", 2:11)
    fe31 <- allFitnessEffects(noIntGenes = ni,
                              drvNames = "g4") 
    set.seed(1)
    mue11 <- oncoSimulIndiv(fe31, 
                            mu = 1e-6,
                            initSize = 1e5,
                            model = "McFL",
                            detectionSize = 5e6,
                            finalTime = 5, 
                            onlyCancer = FALSE, detectionProb = NA)
    expect_true(mue11$MaxNumDrivers == 0)
    expect_true(mue11$MaxDriversLast == 0)

    set.seed(1)
    ni <- runif(10, min = -0.01, max = 0.1)
    names(ni) <- paste0("g", 2:11)
    fe31 <- allFitnessEffects(noIntGenes = ni,
                              drvNames = "g5") 
    set.seed(1)
    mue11 <- oncoSimulIndiv(fe31, 
                            mu = 1e-6,
                            initSize = 1e5,
                            model = "McFL",
                            detectionSize = 5e6,
                            finalTime = 5, 
                            onlyCancer = FALSE, detectionProb = NA)
    expect_true(mue11$MaxNumDrivers == 0)
    expect_true(mue11$MaxDriversLast == 0)

    set.seed(1)
    ni <- runif(10, min = -0.01, max = 0.1)
    names(ni) <- paste0("g", 2:11)
    fe31 <- allFitnessEffects(noIntGenes = ni,
                              drvNames = "g6") 
    set.seed(1)
    mue11 <- oncoSimulIndiv(fe31, 
                            mu = 1e-6,
                            initSize = 1e5,
                            model = "McFL",
                            detectionSize = 5e6,
                            finalTime = 5, 
                            onlyCancer = FALSE, detectionProb = NA)
    expect_true(mue11$MaxNumDrivers == 1)
    expect_true(mue11$MaxDriversLast == 0)

        set.seed(1)
    ni <- runif(10, min = -0.01, max = 0.1)
    names(ni) <- paste0("g", 2:11)
    fe31 <- allFitnessEffects(noIntGenes = ni,
                              drvNames = "g7") 
    set.seed(1)
    mue11 <- oncoSimulIndiv(fe31, 
                            mu = 1e-6,
                            initSize = 1e5,
                            model = "McFL",
                            detectionSize = 5e6,
                            finalTime = 5, 
                            onlyCancer = FALSE, detectionProb = NA)
    expect_true(mue11$MaxNumDrivers == 0)
    expect_true(mue11$MaxDriversLast == 0)

        set.seed(1)
    ni <- runif(10, min = -0.01, max = 0.1)
    names(ni) <- paste0("g", 2:11)
    fe31 <- allFitnessEffects(noIntGenes = ni,
                              drvNames = "g8") 
    set.seed(1)
    mue11 <- oncoSimulIndiv(fe31, 
                            mu = 1e-6,
                            initSize = 1e5,
                            model = "McFL",
                            detectionSize = 5e6,
                            finalTime = 5, 
                            onlyCancer = FALSE, detectionProb = NA)
    expect_true(mue11$MaxNumDrivers == 1)
    expect_true(mue11$MaxDriversLast == 1)

    set.seed(1)
    ni <- runif(10, min = -0.01, max = 0.1)
    names(ni) <- paste0("g", 2:11)
    fe31 <- allFitnessEffects(noIntGenes = ni,
                              drvNames = "g9") 
    set.seed(1)
    mue11 <- oncoSimulIndiv(fe31, 
                            mu = 1e-6,
                            initSize = 1e5,
                            model = "McFL",
                            detectionSize = 5e6,
                            finalTime = 5, 
                            onlyCancer = FALSE, detectionProb = NA)
    expect_true(mue11$MaxNumDrivers == 1)
    expect_true(mue11$MaxDriversLast == 0)

    set.seed(1)
    ni <- runif(10, min = -0.01, max = 0.1)
    names(ni) <- paste0("g", 2:11)
    fe31 <- allFitnessEffects(noIntGenes = ni,
                              drvNames = "g10") 
    set.seed(1)
    mue11 <- oncoSimulIndiv(fe31, 
                            mu = 1e-6,
                            initSize = 1e5,
                            model = "McFL",
                            detectionSize = 5e6,
                            finalTime = 5, 
                            onlyCancer = FALSE, detectionProb = NA)
    expect_true(mue11$MaxNumDrivers == 1)
    expect_true(mue11$MaxDriversLast == 1)

    set.seed(1)
    ni <- runif(10, min = -0.01, max = 0.1)
    names(ni) <- paste0("g", 2:11)
    fe31 <- allFitnessEffects(noIntGenes = ni,
                              drvNames = "g11") 
    set.seed(1)
    mue11 <- oncoSimulIndiv(fe31, 
                            mu = 1e-6,
                            initSize = 1e5,
                            model = "McFL",
                            detectionSize = 5e6,
                            finalTime = 5, 
                            onlyCancer = FALSE, detectionProb = NA)
    expect_true(mue11$MaxNumDrivers == 0)
    expect_true(mue11$MaxDriversLast == 0)

    set.seed(1)
    ni <- runif(10, min = -0.01, max = 0.1)
    names(ni) <- paste0("g", 2:11)
    fe31 <- allFitnessEffects(noIntGenes = ni,
                              drvNames = "g2") 
    set.seed(1)
    mue11 <- oncoSimulIndiv(fe31, 
                            mu = 1e-6,
                            initSize = 1e5,
                            model = "McFL",
                            detectionSize = 5e6,
                            finalTime = 5, 
                            onlyCancer = FALSE, detectionProb = NA)
    expect_true(mue11$MaxNumDrivers == 1)
    expect_true(mue11$MaxDriversLast == 0)

    set.seed(1)
    ni <- runif(10, min = -0.01, max = 0.1)
    names(ni) <- paste0("g", 2:11)
    fe31 <- allFitnessEffects(noIntGenes = ni,
                              drvNames = c("g2", "g4"))
    set.seed(1)
    mue11 <- oncoSimulIndiv(fe31, 
                            mu = 1e-6,
                            initSize = 1e5,
                            model = "McFL",
                            detectionSize = 5e6,
                            finalTime = 5, 
                            onlyCancer = FALSE, detectionProb = NA)
    expect_true(mue11$MaxNumDrivers == 1)
    expect_true(mue11$MaxDriversLast == 0)

    set.seed(1)
    ni <- runif(10, min = -0.01, max = 0.1)
    names(ni) <- paste0("g", 2:11)
    fe31 <- allFitnessEffects(noIntGenes = ni,
                              drvNames = c("g4", "g5"))
    set.seed(1)
    mue11 <- oncoSimulIndiv(fe31, 
                            mu = 1e-6,
                            initSize = 1e5,
                            model = "McFL",
                            detectionSize = 5e6,
                            finalTime = 5, 
                            onlyCancer = FALSE, detectionProb = NA)
    expect_true(mue11$MaxNumDrivers == 0)
    expect_true(mue11$MaxDriversLast == 0)

    set.seed(1)
    ni <- runif(10, min = -0.01, max = 0.1)
    names(ni) <- paste0("g", 2:11)
    fe31 <- allFitnessEffects(noIntGenes = ni,
                              drvNames = c("g3", "g10"))
    set.seed(1)
    mue11 <- oncoSimulIndiv(fe31, 
                            mu = 1e-6,
                            initSize = 1e5,
                            model = "McFL",
                            detectionSize = 5e6,
                            finalTime = 5, 
                            onlyCancer = FALSE, detectionProb = NA)
    expect_true(mue11$MaxNumDrivers == 1)
    expect_true(mue11$MaxDriversLast == 1)


    set.seed(1)
    ni <- runif(10, min = -0.01, max = 0.1)
    names(ni) <- paste0("g", 2:11)
    fe31 <- allFitnessEffects(noIntGenes = ni,
                              drvNames = c("g5", "g9"))
    set.seed(1124)
    mue11 <- oncoSimulIndiv(fe31, 
                            mu = 1e-6,
                            initSize = 1e5,
                            model = "McFL",
                            detectionSize = 1e8,
                            finalTime = 120, 
                            onlyCancer = FALSE, detectionProb = NA)
    expect_true(mue11$MaxNumDrivers == 2)
    expect_true(mue11$MaxDriversLast == 2)
    mue11

    set.seed(1)
    ni <- runif(10, min = -0.01, max = 0.1)
    names(ni) <- paste0("g", 2:11)
    fe31 <- allFitnessEffects(noIntGenes = ni,
                              drvNames = c("g9", "g5"))
    set.seed(1124)
    mue11 <- oncoSimulIndiv(fe31, 
                            mu = 1e-6,
                            initSize = 1e5,
                            model = "McFL",
                            detectionSize = 5e6,
                            finalTime = 130, 
                            onlyCancer = FALSE, detectionProb = NA)
    expect_true(mue11$MaxNumDrivers == 2)
    expect_true(mue11$MaxDriversLast == 1)
    mue11

    set.seed(1)
    ni <- runif(10, min = -0.01, max = 0.1)
    names(ni) <- paste0("g", 2:11)
    fe31 <- allFitnessEffects(noIntGenes = ni,
                              drvNames = c("g6", "g7", "g8", "g9"))
    set.seed(4)
    mue11 <- oncoSimulIndiv(fe31, 
                            mu = 1e-6,
                            initSize = 1e5,
                            model = "McFL",
                            detectionSize = 5e6,
                            finalTime = 90, 
                            onlyCancer = FALSE, detectionProb = NA)
    expect_true(mue11$MaxNumDrivers == 1)
    expect_true(mue11$MaxDriversLast == 1)
    mue11


    ## A different fe
    set.seed(1) ## for reproducibility
    ni <- c(rep(0, 7), runif(10, min = -0.01, max = 0.1))
    names(ni) <- c("a1", "a2", "b1", "b2", "b3", "c1", "c2",
                   paste0("g", 1:10))
    fe31 <- allFitnessEffects(noIntGenes = ni,
                              drvNames = c("g1"))
    set.seed(1) ## so that it is easy to reproduce
    mue11 <- oncoSimulIndiv(fe31, 
                            mu = 1e-6,
                            initSize = 1e6,
                            model = "McFL",
                            detectionSize = 5e6,
                            finalTime = 2,
                            onlyCancer = FALSE, detectionProb = NA)
    expect_true(mue11$MaxNumDrivers == 1)
    expect_true(mue11$MaxDriversLast == 1)
    
})
date()


set.seed(NULL)
