cat(paste("\n Starting sample-prob", date(), "\n"))

p.value.threshold <- 0.05

test_that("Increasing cPDetect decreases time" , {
    gi <- rep(0.1, 10)
    names(gi) <- letters[1:10]
    oi <- allFitnessEffects(noIntGenes = gi)
    n <- 15
    max.tries <- 4  
    for(tries in 1:max.tries) {
        sa <- oncoSimulPop(n,
                           oi,
                           model = "McFL",
                           initSize = 2000,
                           keepEvery = -9,
                           p2 = NULL, n2 = NULL, cPDetect = 1e-4,
                           finalTime = 100000,
                           checkSizePEvery = 20,
                           PDBaseline = 1100,
                           onlyCancer = TRUE,
                           detectionDrivers = 99)
        sb <- oncoSimulPop(n,
                           oi,
                           model = "McFL",
                           initSize = 2000,
                           keepEvery = -9,
                           p2 = NULL, n2 = NULL, cPDetect = 1e-2,
                           finalTime = 100000,
                           checkSizePEvery = 20,
                           PDBaseline = 1100,
                           onlyCancer = TRUE,
                           detectionDrivers = 99)
        ta <- unlist(lapply(sa, function(x) x$FinalTime))
        tb <- unlist(lapply(sb, function(x) x$FinalTime))         
        T1 <- (wilcox.test(ta, tb, alternative = "greater")$p.value < p.value.threshold)
        if(T1) break;
    }
    cat(paste("\n done tries", tries, "\n"))
    expect_true(T1)
})


test_that("Increasing p2 decreases time" , {
    gi <- rep(0.1, 10)
    names(gi) <- letters[1:10]
    oi <- allFitnessEffects(noIntGenes = gi)
    n <- 20
    max.tries <- 4  
    for(tries in 1:max.tries) {
        sa <- oncoSimulPop(n,
                           oi,
                           model = "McFL",
                           initSize = 2000,
                           keepEvery = -9,
                           p2 = .2, n2 = 3500, cPDetect = NULL,
                           finalTime = 100000,
                           checkSizePEvery = 13,
                           PDBaseline = 2000,
                           onlyCancer = TRUE,
                           detectionDrivers = 99)
        sb <- oncoSimulPop(n,
                           oi,
                           model = "McFL",
                           initSize = 2000,
                           keepEvery = -9,
                           p2 = .8, n2 = 3500, cPDetect = NULL,
                           finalTime = 100000,
                           checkSizePEvery = 13,
                           PDBaseline = 2000,
                           onlyCancer = TRUE,
                           detectionDrivers = 99)
        (ta <- unlist(lapply(sa, function(x) x$FinalTime)))
        (tb <- unlist(lapply(sb, function(x) x$FinalTime)))         
        T1 <- (wilcox.test(ta, tb, alternative = "greater")$p.value < p.value.threshold)
        if(T1) break;
    }
    cat(paste("\n done tries", tries, "\n"))
    expect_true(T1)
})



test_that("Increasing n2 increases time" , {
    gi <- rep(0.1, 10)
    names(gi) <- letters[1:10]
    oi <- allFitnessEffects(noIntGenes = gi)
    n <- 20
    max.tries <- 4  
    for(tries in 1:max.tries) {
        sa <- oncoSimulPop(n,
                           oi,
                           model = "McFL",
                           initSize = 2000,
                           keepEvery = -9,
                           p2 = .1, n2 = 5000, cPDetect = NULL,
                           finalTime = 100000,
                           checkSizePEvery = 20,
                           PDBaseline = 1500,
                           onlyCancer = TRUE,
                           detectionDrivers = 99)
        sb <- oncoSimulPop(n,
                           oi,
                           model = "McFL",
                           initSize = 2000,
                           keepEvery = -9,
                           p2 = .1, n2 = 2001, cPDetect = NULL,
                           finalTime = 100000,
                           checkSizePEvery = 20,
                           PDBaseline = 1500,
                           onlyCancer = TRUE,
                           detectionDrivers = 99)
        (ta <- unlist(lapply(sa, function(x) x$FinalTime)))
        (tb <- unlist(lapply(sb, function(x) x$FinalTime)))         
        T1 <- (wilcox.test(ta, tb, alternative = "greater")$p.value < p.value.threshold)
        if(T1) break;
    }
    cat(paste("\n done tries", tries, "\n"))
    expect_true(T1)
})




test_that("Increasing checkSizePEvery increases time" , {
    gi <- rep(0.1, 10)
    names(gi) <- letters[1:10]
    oi <- allFitnessEffects(noIntGenes = gi)
    n <- 20
    max.tries <- 4  
    for(tries in 1:max.tries) {
        sa <- oncoSimulPop(n,
                           oi,
                           model = "McFL",
                           initSize = 2000,
                           keepEvery = -9,
                           p2 = .1, n2 = 1500, cPDetect = NULL,
                           finalTime = 100000,
                           checkSizePEvery = 50,
                           PDBaseline = 1100,
                           onlyCancer = TRUE,
                           detectionDrivers = 99)
        sb <- oncoSimulPop(n,
                           oi,
                           model = "McFL",
                           initSize = 2000,
                           keepEvery = -9,
                           p2 = .1, n2 = 1500, cPDetect = NULL,
                           finalTime = 100000,
                           checkSizePEvery = 10,
                           PDBaseline = 1100,
                           onlyCancer = TRUE,
                           detectionDrivers = 99)
        (ta <- unlist(lapply(sa, function(x) x$FinalTime)))
        (tb <- unlist(lapply(sb, function(x) x$FinalTime)))         
        T1 <- (wilcox.test(ta, tb, alternative = "greater")$p.value < p.value.threshold)
        if(T1) break;
    }
    cat(paste("\n done tries", tries, "\n"))
    expect_true(T1)
})



## Same, with Exp


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
                           onlyCancer = TRUE,
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
                           onlyCancer = TRUE,
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
                           onlyCancer = TRUE,
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
                           onlyCancer = TRUE,
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
    n <- 30
    max.tries <- 4  
    for(tries in 1:max.tries) {
        sa <- oncoSimulPop(n,
                           oi,
                           model = "Exp",
                           initSize = 2000,
                           keepEvery = -9,
                           p2 = .1, n2 = 9000, cPDetect = NULL,
                           finalTime = 100000,
                           checkSizePEvery = 20,
                           PDBaseline = 2100,
                           onlyCancer = TRUE,
                           detectionDrivers = 99)
        sb <- oncoSimulPop(n,
                           oi,
                           model = "Exp",
                           initSize = 2000,
                           keepEvery = -9,
                           p2 = .1, n2 = 2200, cPDetect = NULL,
                           finalTime = 100000,
                           checkSizePEvery = 20,
                           PDBaseline = 2100,
                           onlyCancer = TRUE,
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
                           onlyCancer = TRUE,
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
                           onlyCancer = TRUE,
                           detectionDrivers = 99)
        (ta <- unlist(lapply(sa, function(x) x$FinalTime)))
        (tb <- unlist(lapply(sb, function(x) x$FinalTime)))         
        T1 <- (wilcox.test(ta, tb, alternative = "greater")$p.value < p.value.threshold)
        if(T1) break;
    }
    cat(paste("\n done tries", tries, "\n"))
    expect_true(T1)
})








## And there is no need for fitness effects

test_that("Increasing cPDetect decreases time" , {
    gi <- rep(0.0,  10)
    names(gi) <- letters[1:10]
    oi <- allFitnessEffects(noIntGenes = gi)
    n <- 15
    max.tries <- 4  
    for(tries in 1:max.tries) {
        sa <- oncoSimulPop(n,
                           oi,
                           model = "McFL",
                           initSize = 2000,
                           keepEvery = -9,
                           p2 = NULL, n2 = NULL, cPDetect = 1e-4,
                           finalTime = 100000,
                           checkSizePEvery = 10,
                           PDBaseline = 1100,
                           onlyCancer = FALSE,
                           detectionDrivers = 99)
        sb <- oncoSimulPop(n,
                           oi,
                           model = "McFL",
                           initSize = 2000,
                           keepEvery = -9,
                           p2 = NULL, n2 = NULL, cPDetect = 1e-2,
                           finalTime = 100000,
                           checkSizePEvery = 1,
                           PDBaseline = 1100,
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


test_that("Increasing p2 decreases time" , {
    gi <- rep(0.0,  10)
    names(gi) <- letters[1:10]
    oi <- allFitnessEffects(noIntGenes = gi)
    n <- 20
    max.tries <- 4  
    for(tries in 1:max.tries) {
        sa <- oncoSimulPop(n,
                           oi,
                           model = "McFL",
                           initSize = 2000,
                           keepEvery = -9,
                           p2 = .1, n2 = 3500, cPDetect = NULL,
                           finalTime = 100000,
                           checkSizePEvery = 10,
                           PDBaseline = 1100,
                           onlyCancer = FALSE,
                           detectionDrivers = 99)
        sb <- oncoSimulPop(n,
                           oi,
                           model = "McFL",
                           initSize = 2000,
                           keepEvery = -9,
                           p2 = .6, n2 = 3500, cPDetect = NULL,
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



test_that("Increasing n2 increases time" , {
    gi <- rep(0.0,  10)
    names(gi) <- letters[1:10]
    oi <- allFitnessEffects(noIntGenes = gi)
    n <- 20
    max.tries <- 4  
    for(tries in 1:max.tries) {
        sa <- oncoSimulPop(n,
                           oi,
                           model = "McFL",
                           initSize = 2000,
                           keepEvery = -9,
                           p2 = .1, n2 = 5000, cPDetect = NULL,
                           finalTime = 100000,
                           checkSizePEvery = 5,
                           PDBaseline = 1500,
                           onlyCancer = FALSE,
                           detectionDrivers = 99)
        sb <- oncoSimulPop(n,
                           oi,
                           model = "McFL",
                           initSize = 2000,
                           keepEvery = -9,
                           p2 = .1, n2 = 2001, cPDetect = NULL,
                           finalTime = 100000,
                           checkSizePEvery = 5,
                           PDBaseline = 1500,
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




test_that("Increasing checkSizePEvery increases time" , {
    gi <- rep(0.0,  10)
    names(gi) <- letters[1:10]
    oi <- allFitnessEffects(noIntGenes = gi)
    n <- 20
    max.tries <- 4  
    for(tries in 1:max.tries) {
        sa <- oncoSimulPop(n,
                           oi,
                           model = "McFL",
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
                           model = "McFL",
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










    ## gi2 <- rep(0, 5)
    ## names(gi2) <- letters[1:5]
    ## oi2 <- allFitnessEffects(noIntGenes = gi2)
    ## ## nicely exponential, as expected
    ## s5 <- oncoSimulPop(100,
    ##                    oi2,
    ##                    model = "McFL",
    ##                    initSize = 1000,
    ##                    detectionSize = 4800,
    ##                    finalTime = 100000, ## crucial this is increased
    ##                    keepEvery = -9,
    ##                    p2 = .1,
    ##                    checkSizePEvery = 2,
    ##                    verbosity = 0,
    ##                    PDBaseline = 1000,
    ##                    onlyCancer = TRUE,
    ##                    detectionDrivers = 99)
    ## s5
    ## hist(unlist(lapply(s5, function(x) x$FinalTime)))

cat(paste("\n Ending sample-prob tests", date(), "\n"))
