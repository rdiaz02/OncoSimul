## this works OK
runif(1)
the.seed <- .Random.seed
runif(1)

cat("\n set seed\n")
set.seed(1)
runif(3)

cat("\n reset the seed\n")
.Random.seed <- the.seed
runif(1)





runif(1)
the.seed <- .Random.seed
cat("after \n")
runif(1)


## test_that("this is a dummy test", {
##     cat("\n set seed\n")
##     set.seed(1)
##     runif(3)
##     expect_true(TRUE)

##     set.seed(1) ## for reproducibility
##     ni <- c(rep(0, 7), runif(10, min = -0.01, max = 0.1))
##     names(ni) <- c("a1", "a2", "b1", "b2", "b3", "c1", "c2",
##                    paste0("g", 1:10))
##     fe31 <- allFitnessEffects(noIntGenes = ni,
##                               drvNames = c("g1"))
##     set.seed(1) ## so that it is easy to reproduce
##     mue11 <- oncoSimulIndiv(fe31, 
##                             mu = 1e-6,
##                             initSize = 1e6,
##                             model = "McFL",
##                             detectionSize = 5e6,
##                             finalTime = 2,
##                             onlyCancer = FALSE)
##     expect_true(mue11$MaxNumDrivers == 1)
##     expect_true(mue11$MaxDriversLast == 1)

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

cat("\n reset the seed\n")
.Random.seed <- the.seed
runif(1)

save(file = "aSeed.RData", the.seed)

