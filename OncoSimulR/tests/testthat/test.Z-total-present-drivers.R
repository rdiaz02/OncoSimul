inittime <- Sys.time()
cat(paste("\n Starting Z-total-present-drivers tests", date(), "\n"))


test_that("Count TotalPresentDrivers", {
    ## There once was a bug here (in the driverCounts C++ function)
    ## Testing we are OK
    ## Below, I use specific seeds for a few cases where I got
    ## errors. This actual cases might not be the same with other
    ## compilers. But we are checking anyway against the number obtained
    ## by other means.
    RNGkind("Mersenne-Twister")
    set.seed(1)
    lni <- 5 ## no fitness effects genes
    fni <- 50 ## fitness effects genes
    no <- 1e2
    ft <- 2
    s3 <- 2.0
    mu <- 1e-5 
    ni <- rep(0, lni)
    gn <- paste(paste0("a", 1:fni), collapse = ", ")
    set.seed(2)
    ## There are no drivers here
    f3 <- allFitnessEffects(epistasis = c("A" = s3),
                            geneToModule = c("A" = gn),
                            noIntGenes = ni)
    s3.ng <- oncoSimulIndiv(f3,
                            mu = mu,
                            mutationPropGrowth = FALSE,
                            finalTime =ft,
                            initSize = no,
                            sampleEvery = 0.03,
                            keepEvery = 1,
                            onlyCancer = FALSE, detectionProb = NA,
                            seed = NULL)
    s3.ng
    ## so this should be zero
    s3.ng$TotalPresentDrivers
    expect_true(s3.ng$TotalPresentDrivers ==
                length(strsplit(s3.ng$OccurringDrivers, ", ")[[1]]))
    set.seed(1)
    f4 <- allFitnessEffects(epistasis = c("A" = s3),
                            geneToModule = c("A" = gn),
                            noIntGenes = ni
                          , drvNames = c("a2", "a5")
                            )
    s4.ng <- oncoSimulIndiv(f4,
                            mu = mu,
                            mutationPropGrowth = FALSE,
                            finalTime =ft,
                            initSize = no,
                            sampleEvery = 0.03,
                            keepEvery = 1,
                            onlyCancer = FALSE, detectionProb = NA,
                            ## verbosity = 4,
                            seed = NULL)
    ## this should be zero
    s4.ng$TotalPresentDrivers
    s4.ng
    expect_true(s4.ng$TotalPresentDrivers ==
                length(strsplit(s4.ng$OccurringDrivers, ", ")[[1]]))
    ## and something with drivers
    set.seed(53)
    s3.l <- oncoSimulIndiv(f3,
                           mu = mu,
                           mutationPropGrowth = FALSE,
                           finalTime = 28,
                           initSize = no,
                           sampleEvery = 0.03,
                           keepEvery = 1,
                           onlyCancer = FALSE, detectionProb = NA,
                           ## verbosity = 4,
                           seed = NULL)
    s3.l
    expect_true(s3.l$TotalPresentDrivers ==
                length(strsplit(s3.l$OccurringDrivers, ", ")[[1]]))
    set.seed(53)
    s3.l <- oncoSimulIndiv(f3,
                           mu = mu,
                           mutationPropGrowth = FALSE,
                           finalTime = 35,
                           initSize = no,
                           sampleEvery = 0.03,
                           keepEvery = 1,
                           onlyCancer = FALSE, detectionProb = NA,
                           ## verbosity = 4,
                           seed = NULL)
    s3.l
    expect_true(s3.l$TotalPresentDrivers ==
                length(strsplit(s3.l$OccurringDrivers, ", ")[[1]]))
})
set.seed(NULL)
set.seed(NULL)
cat(paste("\n Ending Z-total-present-drivers tests", date(), "\n"))


cat(paste("  Took ", round(difftime(Sys.time(), inittime, units = "secs"), 2), "\n\n"))
rm(inittime)
