inittime <- Sys.time()
cat(paste("\n Starting test.Z-fixation at", date(), "\n"))

## Some tests below might only work on Linux because of compiler
## differences, because the rng is done in C++, etc.
## Note that the difference is in whether a certain code
## is exercised. The runs should work in all platforms, though.
## This is the same as in test.Z-oncoSimulIndiv.R


test_that("Three cases with fixation of genotypes" ,{
    load(system.file("testdata_fee.RData", package="OncoSimulR"))
    i <- 3
    ng <- 7
    RNGkind("Mersenne-Twister")

    set.seed(i)
    r3 <- oncoSimulIndiv( fp = fee,
                         model = "McFL",
                         initSize = 2000,
                         mu = 1e-4,
                         detectionSize = NA,
                         sampleEvery = .03,
                         keepEvery = 1,
                         finalTime = 50000,
                         fixation = unlist(feex[["labelled_peaks"]]), 
                         detectionDrivers = NA,
                         detectionProb = NA,
                         onlyCancer = TRUE,
                         max.num.tries = 500,
                         max.wall.time = 20, 
                         errorHitMaxTries = TRUE,
                         keepPhylog = FALSE)
    
    summary(r3)
    if(Sys.info()["sysname"] == "Linux") {
    ## stopping at ABCG, which is not a maximum, not a labelled peak
        expect_equal(
            c(A = 1, B = 1, C = 1, D = 0, E = 0, F = 0, G = 1),
            samplePop(r3, "last", "singleCell")[1, ])
    }
    feex$labelled_peaks
    ## And not a single genotype unique
    expect_true(r3$TotalPopSize > r3$LargestClone)
    
    ## stop on genotypes
    set.seed(i)
    r4 <- oncoSimulIndiv( fp = fee,
                     model = "McFL",
                     initSize = 2000,
                     mu = 1e-4,
                     detectionSize = NA,
                     sampleEvery = .03,
                     keepEvery = 1,
                     finalTime = 50000,
                     fixation = paste0("_,", unlist(feex[["labelled_peaks"]])),
                     detectionDrivers = NA,
                     detectionProb = NA,
                     onlyCancer = TRUE,
                     max.num.tries = 500,
                     max.wall.time = 20, 
                     errorHitMaxTries = TRUE,
                     keepPhylog = FALSE)
    summary(r4)
    ## now a local max is the single clone
    if(Sys.info()["sysname"] == "Linux") {
        expect_equal(
            c(A = 1, B = 1, C = 1, D = 0, E = 0, F = 1, G = 1),
            samplePop(r4, "last", "singleCell")[1, ])
    }
    expect_true(r4$TotalPopSize == r4$LargestClone)


    ## tolerance
    set.seed(i)
    r5 <- oncoSimulIndiv( fp = fee,
                         model = "McFL",
                         initSize = 2000,
                         mu = 1e-4,
                         detectionSize = NA,
                         sampleEvery = .03,
                         keepEvery = 1,
                         finalTime = 50000,
                         fixation = c(paste0("_,", unlist(feex[["labelled_peaks"]])),
                                      fixation_tolerance = 0.05),
                         detectionDrivers = NA,
                         detectionProb = NA,
                         onlyCancer = TRUE,
                         max.num.tries = 500,
                         max.wall.time = 20, 
                         errorHitMaxTries = TRUE,
                         keepPhylog = FALSE)
    summary(r5)
    if(Sys.info()["sysname"] == "Linux") {
    expect_equal(
        c(A = 1, B = 1, C = 1, D = 0, E = 0, F = 1, G = 1),
        samplePop(r5, "last", "singleCell")[1, ])
    }
    ## but other clones present too
    expect_true(r5$TotalPopSize > r5$LargestClone)
})


set.seed(NULL)
cat(paste("\n Ending test.Z-fixation at", date(), "\n"))
cat(paste("  Took ", round(difftime(Sys.time(), inittime, units = "secs"), 2), "\n\n"))
rm(inittime)
