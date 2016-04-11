cat(paste("\n Starting driverCounts long at", date()))

date()
test_that("Runs without crashes", {
    
    iteration <- 1
    ni <- runif(10, min = -0.01, max = 0.1)
    names(ni) <- paste0("g", 2:11)
    ## each single one
    for(i in 2:11) {
        cat("\n doing iteration", iteration, "\n")
        ni[] <- runif(10, min = -0.01, max = 0.1)
        ni <- sample(ni)
        drvN <- paste0("g", i)
        fe31 <- allFitnessEffects(noIntGenes = ni,
                                  drvNames = drvN)
        expect_silent(uu <- oncoSimulPop(20,
                                         fe31, 
                                         mu = 1e-6,
                                         initSize = 1e5,
                                         model = "McFL",
                                         detectionSize = 5e6,
                                         finalTime = 5,
                                         keepEvery = 1,
                                         onlyCancer = FALSE,
                                         mc.cores = 2,
                                         sampleEvery = 0.03
                                         ))
        iteration <- iteration + 1
    }
    ## all pairs, trios, etc, shuffled order
    for(n in 2:10) {
        m <- combinations(10, n, 2:11)
        for(j in 1:nrow(m)) {
            cat("\n doing iteration", iteration, "\n")
            ni[] <- runif(10, min = -0.01, max = 0.1)
            ni <- sample(ni)
            drvN <- paste0("g", sample(m[j, ]))
            fe31 <- allFitnessEffects(noIntGenes = ni,
                                  drvNames = drvN)
            expect_silent(uu <- oncoSimulPop(10,
                                             fe31, 
                                             mu = 1e-6,
                                             initSize = 1e5,
                                             model = "McFL",
                                             detectionSize = 5e6,
                                             finalTime = 5,
                                             keepEvery = 1,
                                             onlyCancer = FALSE,
                                             mc.cores = 2,
                                             sampleEvery = 0.03
                                             ))
            iteration <- iteration + 1
        }
    }

})
date()





## The following tests is likely to crash in non-Linux and/or non-gcc
## machines (e.g., using clang), as the details depend on the random
## numbers generated in R and C++.


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
                            onlyCancer = FALSE)
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
                            onlyCancer = FALSE)
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
                            onlyCancer = FALSE)
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
                            onlyCancer = FALSE)
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
                            onlyCancer = FALSE)
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
                            onlyCancer = FALSE)
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
                            onlyCancer = FALSE)
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
                            onlyCancer = FALSE)
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
                            onlyCancer = FALSE)
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
                            onlyCancer = FALSE)
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
                            onlyCancer = FALSE)
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
                            onlyCancer = FALSE)
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
                            onlyCancer = FALSE)
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
                            onlyCancer = FALSE)
    expect_true(mue11$MaxNumDrivers == 1)
    expect_true(mue11$MaxDriversLast == 1)


    set.seed(1)
    ni <- runif(10, min = -0.01, max = 0.1)
    names(ni) <- paste0("g", 2:11)
    fe31 <- allFitnessEffects(noIntGenes = ni,
                              drvNames = c("g5", "g3"))
    set.seed(4)
    mue11 <- oncoSimulIndiv(fe31, 
                            mu = 1e-6,
                            initSize = 1e5,
                            model = "McFL",
                            detectionSize = 5e6,
                            finalTime = 90, 
                            onlyCancer = FALSE)
    expect_true(mue11$MaxNumDrivers == 2)
    expect_true(mue11$MaxDriversLast == 2)

    set.seed(1)
    ni <- runif(10, min = -0.01, max = 0.1)
    names(ni) <- paste0("g", 2:11)
    fe31 <- allFitnessEffects(noIntGenes = ni,
                              drvNames = c("g6", "g5"))
    set.seed(4)
    mue11 <- oncoSimulIndiv(fe31, 
                            mu = 1e-6,
                            initSize = 1e5,
                            model = "McFL",
                            detectionSize = 5e6,
                            finalTime = 90, 
                            onlyCancer = FALSE)
    expect_true(mue11$MaxNumDrivers == 2)
    expect_true(mue11$MaxDriversLast == 1)

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
                            onlyCancer = FALSE)
    expect_true(mue11$MaxNumDrivers == 1)
    expect_true(mue11$MaxDriversLast == 1)


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
                            onlyCancer = FALSE)
    expect_true(mue11$MaxNumDrivers == 1)
    expect_true(mue11$MaxDriversLast == 1)
    
})
date()



cat(paste("\n Ending driverCounts long at", date()))
