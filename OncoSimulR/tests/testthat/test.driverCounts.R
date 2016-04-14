cat(paste("\n Starting driverCounts at", date()))
## RNGkind("Mersenne-Twister") ## we have some examples below with random

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


cat(paste("\n Ending driverCounts at", date(), "\n"))
