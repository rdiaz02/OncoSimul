cat(paste("\n Starting driverCounts long at", date()))
cat(paste("\n             a runif", runif(1), "\n"))
cat(paste("\n             a runif", runif(1), "\n"))
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
                                         onlyCancer = FALSE, detectionProb = NA,
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
                                             onlyCancer = FALSE, detectionProb = NA,
                                             mc.cores = 2,
                                             sampleEvery = 0.03
                                             ))
            iteration <- iteration + 1
        }
    }

})
date()



cat(paste("\n Ending driverCounts long at", date()))






