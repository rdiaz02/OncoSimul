cat(paste("\n Starting sample-prob", date(), "\n"))


mcflEv <- function(p, s, initSize) {
    ## expects vectors for p and s
    K <- initSize/(exp(1) - 1)
    
    ## Expected number at equilibrium
    return( K * (exp(prod((1 + s)^p)) - 1))
}

mcflEvF <- function(F, initSize) {
    ## As mcflEv, but with fitness
    ## expects vectors for p and s
    K <- initSize/(exp(1) - 1)
    
    ## Expected number at equilibrium
    return( K * (exp(F) - 1))
}



test_that("some runs", {

    library(OncoSimulR)
    gi <- rep(0.1, 10)
    names(gi) <- letters[1:10]
    oi <- allFitnessEffects(noIntGenes = gi)

    s1 <- oncoSimulIndiv(oi,
                         model = "McFL",
                         initSize = 1000,
                         detectionSize = 1800,
                         n2 = NULL, p2 = NULL, cPDetect = NULL,
                         keepEvery = 5,
                         onlyCancer = TRUE,
                         detectionDrivers = 99)
    s1

    s2 <- oncoSimulIndiv(oi,
                         model = "McFL",
                         initSize = 1000,
                         detectionSize = 4800,
                         keepEvery = 5,
                         verbosity = 1,
                         p2 = .1,
                         checkSizePEvery = 10,
                         PDBaseline = 1100,
                         onlyCancer = TRUE,
                         detectionDrivers = 99)
    s2


    s3 <- oncoSimulIndiv(oi,
                         model = "McFL",
                         initSize = 1000,
                         detectionSize = 2200,
                         keepEvery = 5,
                         verbosity = 1,
                         cPDetect = 1e-7, n2 = NULL, p2 = NULL,
                         PDBaseline = 1100,
                         onlyCancer = TRUE,
                         detectionDrivers = 99)
    s3
    

    expect_output(print(oncoSimulIndiv(oi, initSize = 1,
                                 sampleEvery = 0.03,
                                 keepEvery = 5,
                                 onlyCancer = FALSE)),
                  "Individual OncoSimul trajectory", fixed = TRUE)
    
    expect_error(oncoSimulIndiv(oi, initSize = 1,
                                onlyCancer = FALSE, model = "McFL"),
                 "Using McFarland's model: K cannot be < 1",
                 fixed = TRUE)


})





cat(paste("\n Ending sample-prob tests", date(), "\n"))
