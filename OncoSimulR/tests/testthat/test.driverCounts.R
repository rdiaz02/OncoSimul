cat(paste("\n Starting driverCounts at", date()))
RNGkind("Mersenne-Twister") ## we have some examples below with random

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
        expect_output(print(oncoSimulPop(20,
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
                                         )),
                      "Population of OncoSimul trajectories",
                      fixed = TRUE)
        iteration <- iteration + 1
    }
    ## some pairs, trios, etc, shuffled order; all are run in the long
    ## test; this is a paranoid testing set up
    ## The ones run here are the most dangerous ones.
    require(gtools)
    for(n in c(2, 3)) {
        m <- gtools::combinations(10, n, 2:11)
        for(j in sample(nrow(m), 45)) { ## all for 2, a sample of 45 for 3
            cat("\n doing iteration", iteration, "\n")
            ni[] <- runif(10, min = -0.01, max = 0.1)
            ni <- sample(ni)
            drvN <- paste0("g", sample(m[j, ]))
            fe31 <- allFitnessEffects(noIntGenes = ni,
                                  drvNames = drvN)
            expect_output(print(oncoSimulPop(10,
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
                                             )),
                          "Population of OncoSimul trajectories",
                          fixed = TRUE)
            iteration <- iteration + 1
        }
    }

})
date()

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


date()
test_that("Assertion is correct when nothing returned",{
    RNGkind("L'Ecuyer-CMRG") 
    set.seed(13)
    oi <- allFitnessEffects(orderEffects =
                                c("F > D" = -0.3, "D > F" = 0.4),
                            noIntGenes = runif(5, 0.01, 0.06),
                            geneToModule =
                                c("Root" = "Root",
                                  "F" = "f1, f2, f3",
                                  "D" = "d1, d2"),
                            drvNames = c("d1", "d2", "f1", "f2", "f3"))
    set.seed(13)
    expect_message(ou <- oncoSimulSample(1, 
                                  oi,
                                  sampleEvery = 0.03,
                                  onlyCancer = FALSE,
                                  model = "Bozic",
                                  mutationPropGrowth = TRUE,
                                  seed = NULL),
                  "Subjects by Genes", fixed = TRUE)



    RNGkind("Mersenne-Twister") 
    set.seed(13)
    oi <- allFitnessEffects(orderEffects =
                                c("F > D" = -0.3, "D > F" = 0.4),
                            noIntGenes = runif(5, 0.01, 0.06),
                            geneToModule =
                                c("Root" = "Root",
                                  "F" = "f1, f2, f3",
                                  "D" = "d1, d2"),
                            drvNames = c("d1", "d2", "f1", "f2", "f3"))
    set.seed(13)
    expect_message(ou2 <- oncoSimulSample(1, 
                    oi,
                    sampleEvery = 0.03,
                    onlyCancer = FALSE,
                    model = "Bozic",
                    mutationPropGrowth = TRUE,
                    seed = NULL),
                  "Subjects by Genes", fixed = TRUE)

} )
date()


cat(paste("\n Ending driverCounts at", date(), "\n"))
