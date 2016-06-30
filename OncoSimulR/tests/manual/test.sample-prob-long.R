## This ain't particularly long, but adds time and is redundant.
## I leave it here for extra checks.


## Testing when no onlyCancer, but of course we have to exit sooner

cat(paste("\n Starting sample-prob-long", date(), "\n"))

p.value.threshold <- 0.05


test_that("Increasing cPDetect decreases time, Exp" , {
    gi <- rep(0.1, 10)
    names(gi) <- letters[1:10]
    oi <- allFitnessEffects(noIntGenes = gi)
    n <- 20
    max.tries <- 4  
    for(tries in 1:max.tries) {
        sa <- oncoSimulPop(n,
                           oi,
                           model = "Exp",
                           initSize = 1000,
                           keepEvery = -9,
                           p2 = NULL, n2 = NULL, cPDetect = 1e-5,
                           finalTime = 100000,
                           checkSizePEvery = 50,
                           PDBaseline = 500,
                           onlyCancer = FALSE,
                           detectionDrivers = 99)
        sb <- oncoSimulPop(n,
                           oi,
                           model = "Exp",
                           initSize = 1000,
                           keepEvery = -9,
                           p2 = NULL, n2 = NULL, cPDetect = .01,
                           finalTime = 100000,
                           checkSizePEvery = 50,
                           PDBaseline = 500,
                           onlyCancer = FALSE,
                           detectionDrivers = 99)
        ta <- unlist(lapply(sa, function(x) x$FinalTime))
        tb <- unlist(lapply(sb, function(x) x$FinalTime))         
        T1 <- (wilcox.test(ta, tb, alternative = "greater")$p.value < p.value.threshold)
        if(T1) break;
    }
    cat(paste("\n done tries", tries, "\n"))
    expect_true(T1)
})


test_that("Increasing p2 decreases time, Exp" , {
    gi <- rep(0.1, 10)
    names(gi) <- letters[1:10]
    oi <- allFitnessEffects(noIntGenes = gi)
    n <- 20
    max.tries <- 4  
    for(tries in 1:max.tries) {
        sa <- oncoSimulPop(n,
                           oi,
                           model = "Exp",
                           initSize = 2000,
                           keepEvery = -9,
                           p2 = .1, n2 = 8500, cPDetect = NULL,
                           finalTime = 100000,
                           checkSizePEvery = 10,
                           PDBaseline = 1100,
                           onlyCancer = FALSE,
                           detectionDrivers = 99)
        sb <- oncoSimulPop(n,
                           oi,
                           model = "Exp",
                           initSize = 2000,
                           keepEvery = -9,
                           p2 = .8, n2 = 8500, cPDetect = NULL,
                           finalTime = 100000,
                           checkSizePEvery = 10,
                           PDBaseline = 1100,
                           onlyCancer = FALSE,
                           detectionDrivers = 99)
        (ta <- unlist(lapply(sa, function(x) x$FinalTime)))
        (tb <- unlist(lapply(sb, function(x) x$FinalTime)))         
        T1 <- (wilcox.test(ta, tb, alternative = "greater")$p.value < p.value.threshold)
        if(T1) break;
    }
    cat(paste("\n done tries", tries, "\n"))
    expect_true(T1)
})



test_that("Increasing n2 increases time, Exp" , {
    gi <- rep(0.1, 10)
    names(gi) <- letters[1:10]
    oi <- allFitnessEffects(noIntGenes = gi)
    n <- 20
    max.tries <- 4  
    for(tries in 1:max.tries) {
        sa <- oncoSimulPop(n,
                           oi,
                           model = "Exp",
                           initSize = 5000,
                           keepEvery = -9,
                           p2 = .6, n2 = 7000, cPDetect = NULL,
                           finalTime = 100000,
                           checkSizePEvery = 50,
                           PDBaseline = 2100,
                           onlyCancer = FALSE,
                           detectionDrivers = 99)
        sb <- oncoSimulPop(n,
                           oi,
                           model = "Exp",
                           initSize = 5000,
                           keepEvery = -9,
                           p2 = .6, n2 = 2500, cPDetect = NULL,
                           finalTime = 100000,
                           checkSizePEvery = 50,
                           PDBaseline = 2100,
                           onlyCancer = FALSE,
                           detectionDrivers = 99)
        (ta <- unlist(lapply(sa, function(x) x$FinalTime)))
        (tb <- unlist(lapply(sb, function(x) x$FinalTime)))         
        T1 <- (wilcox.test(ta, tb, alternative = "greater")$p.value < p.value.threshold)
        if(T1) break;
    }
    cat(paste("\n done tries", tries, "\n"))
    expect_true(T1)
})




test_that("Increasing checkSizePEvery increases time, Exp" , {
    gi <- rep(0.1, 10)
    names(gi) <- letters[1:10]
    oi <- allFitnessEffects(noIntGenes = gi)
    n <- 20
    max.tries <- 4  
    for(tries in 1:max.tries) {
        sa <- oncoSimulPop(n,
                           oi,
                           model = "Exp",
                           initSize = 2000,
                           keepEvery = -9,
                           p2 = .1, n2 = 1500, cPDetect = NULL,
                           finalTime = 100000,
                           checkSizePEvery = 50,
                           PDBaseline = 1100,
                           onlyCancer = FALSE,
                           detectionDrivers = 99)
        sb <- oncoSimulPop(n,
                           oi,
                           model = "Exp",
                           initSize = 2000,
                           keepEvery = -9,
                           p2 = .1, n2 = 1500, cPDetect = NULL,
                           finalTime = 100000,
                           checkSizePEvery = 10,
                           PDBaseline = 1100,
                           onlyCancer = FALSE,
                           detectionDrivers = 99)
        (ta <- unlist(lapply(sa, function(x) x$FinalTime)))
        (tb <- unlist(lapply(sb, function(x) x$FinalTime)))         
        T1 <- (wilcox.test(ta, tb, alternative = "greater")$p.value < p.value.threshold)
        if(T1) break;
    }
    cat(paste("\n done tries", tries, "\n"))
    expect_true(T1)
})

cat(paste("\n Ending sample-prob-long tests", date(), "\n"))
