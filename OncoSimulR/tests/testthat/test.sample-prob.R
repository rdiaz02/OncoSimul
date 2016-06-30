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
                         keepEvery = 5,
                         onlyCancer = TRUE,
                         detectionSize = 1.1e3,
                         detectionDrivers = 99)

    

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
