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

    set.seed(2)
    
    s2 <- oncoSimulIndiv(oi,
                         model = "McFL",
                         initSize = 1000,
                         detectionSize = 4800,
                         keepEvery = -9,
                         p2 = .1,
                         checkSizePEvery = 10,
                         finalTime = 100000,
                         verbosity = 1,
                         PDBaseline = 1100,
                         onlyCancer = TRUE,
                         detectionDrivers = 99)
    s2

    s4 <- oncoSimulPop(50,
                       oi,
                       model = "McFL",
                       initSize = 1000,
                       detectionSize = 4800,
                       keepEvery = -9,
                       p2 = .05,
                       finalTime = 100000,
                       checkSizePEvery = 50,
                       verbosity = 1,
                       PDBaseline = 1100,
                       onlyCancer = TRUE,
                       detectionDrivers = 99)
    s4

    s6 <- oncoSimulPop(50,
                       oi,
                       model = "McFL",
                       initSize = 20,
                       detectionSize = 4800,
                       keepEvery = -9,
                       p2 = .05,
                       finalTime = 100000,
                       checkSizePEvery = 10,
                       verbosity = 1,
                       onlyCancer = TRUE,
                       detectionDrivers = 99)
    s6


    gi2 <- rep(0, 10)
    names(gi2) <- letters[1:10]
    oi2 <- allFitnessEffects(noIntGenes = gi2)

    ## nicely exponential, as expected
    s5 <- oncoSimulPop(200,
                       oi2,
                       model = "McFL",
                       initSize = 1000,
                       detectionSize = 4800,
                       finalTime = 100000, ## crucial this is increased
                       keepEvery = -9,
                       p2 = .1,
                       checkSizePEvery = 2,
                       verbosity = 0,
                       PDBaseline = 1000,
                       onlyCancer = TRUE,
                       detectionDrivers = 99)
    s5


    

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
