inittime <- Sys.time()
cat(paste("\n Starting sample-prob", date(), "\n"))

p.value.threshold <- 1e-4

## a McFL version in long tests and also below
date()
test_that("Increasing cPDetect decreases time, Exp" , {
    gi <- rep(0.1, 10)
    names(gi) <- letters[1:10]
    oi <- allFitnessEffects(noIntGenes = gi)
    n <- 400 ## I once saw a failure in BioC, Windows 
    max.tries <- 4 
    for(tries in 1:max.tries) {
        sa <- oncoSimulPop(n,
                           oi,
                           model = "Exp",
                           initSize = 1000,
                           keepEvery = -9,
                           detectionProb = c(p2 = NULL, n2 = NULL,
                                             cPDetect = 0.05),
                           finalTime = NA,
                           onlyCancer = FALSE,
                           detectionDrivers = 99, mc.cores = 2)
        sb <- oncoSimulPop(n,
                           oi,
                           model = "Exp",
                           initSize = 1000,
                           keepEvery = -9,
                           detectionProb = c(p2 = NULL, n2 = NULL,
                                             cPDetect = 2),
                           finalTime = NA,
                           onlyCancer = FALSE,
                           detectionDrivers = 99, mc.cores = 2)
        ta <- unlist(lapply(sa, function(x) x$FinalTime))
        tb <- unlist(lapply(sb, function(x) x$FinalTime))
        print(suppressWarnings(wilcox.test(ta, tb, alternative = "greater")$p.value))
        T1 <- suppressWarnings(wilcox.test(ta, tb, alternative = "greater")$p.value < p.value.threshold)
        if(T1) break;
    }
    cat(paste("\n done tries", tries, "\n"))
    expect_true(T1)
})
date()

## McFL in long
date()
test_that("Increasing p2 decreases time, Exp" , {
    gi <- rep(0.1, 10)
    names(gi) <- letters[1:10]
    oi <- allFitnessEffects(noIntGenes = gi)
    n <- 200
    max.tries <- 4  
    for(tries in 1:max.tries) {
        sa <- oncoSimulPop(n,
                           oi,
                           model = "Exp",
                           initSize = 3000,
                           keepEvery = NA,
                           detectionProb = c(p2 = .005, n2 = 5000, checkSizePEvery = 1,  cPDetect = NULL),
                           finalTime = NA,
                           onlyCancer = FALSE,
                           detectionDrivers = 99, mc.cores = 2)
        sb <- oncoSimulPop(n,
                           oi,
                           model = "Exp",
                           initSize = 3000,
                           keepEvery = NA,
                           detectionProb = c(p2 = .8, n2 = 5000, checkSizePEvery = 1,  cPDetect = NULL),
                           finalTime = NA,
                           onlyCancer = FALSE,
                           detectionDrivers = 99, mc.cores = 2)
        (ta <- unlist(lapply(sa, function(x) x$FinalTime)))
        (tb <- unlist(lapply(sb, function(x) x$FinalTime)))
        print(suppressWarnings(wilcox.test(ta, tb, alternative = "greater")$p.value))
        T1 <- suppressWarnings(wilcox.test(ta, tb, alternative = "greater")$p.value < p.value.threshold)
        if(T1) break;
    }
    cat(paste("\n done tries", tries, "\n"))
    expect_true(T1)
})
date()




date()
test_that("Increasing n2 increases time" , {
    gi <- rep(0.1, 10)
    names(gi) <- letters[1:10]
    oi <- allFitnessEffects(noIntGenes = gi)
    n <- 30 ## 70
    max.tries <- 4  
    for(tries in 1:max.tries) {
        sa <- oncoSimulPop(n,
                           oi,
                           model = "McFL",
                           initSize = 2000,
                           keepEvery = NA,
                           detectionProb = c(p2 = .1, n2 = 4000,  checkSizePEvery = 1,
                                             PDBaseline = 1900, cPDetect = NA),
                           finalTime = NA, detectionSize = NA,
                           onlyCancer = TRUE,
                           detectionDrivers = NA, mc.cores = 2)
        sb <- oncoSimulPop(n,
                           oi,
                           model = "McFL",
                           initSize = 2000,
                           keepEvery = NA,
                           detectionProb = c(p2 = .1, n2 = 2001, checkSizePEvery = 1,
                                             PDBaseline = 1900, cPDetect = NA),
                           finalTime = NA, detectionSize = NA,
                           onlyCancer = TRUE,
                           detectionDrivers = NA, mc.cores = 2)
        (ta <- unlist(lapply(sa, function(x) x$FinalTime)))
        (tb <- unlist(lapply(sb, function(x) x$FinalTime)))
        print(suppressWarnings(wilcox.test(ta, tb, alternative = "greater")$p.value))
        T1 <- suppressWarnings(wilcox.test(ta, tb, alternative = "greater")$p.value < p.value.threshold)
        if(T1) break;
    }
    cat(paste("\n done tries", tries, "\n"))
    expect_true(T1)
})
date()



date()
test_that("Increasing checkSizePEvery increases time" , {
    gi <- rep(0.1, 10)
    names(gi) <- letters[1:10]
    oi <- allFitnessEffects(noIntGenes = gi)
    n <- 30
    max.tries <- 4  
    for(tries in 1:max.tries) {
        sa <- oncoSimulPop(n,
                           oi,
                           model = "McFL",
                           initSize = 2000,
                           keepEvery = NA,
                           detectionProb = c(p2 = .1, n2 = 1500, checkSizePEvery = 20,
                                             PDBaseline = 1100, cPDetect = NA),
                           finalTime = NA, detectionSize = NA,
                           onlyCancer = TRUE,
                           detectionDrivers = NA, mc.cores = 2)
        sb <- oncoSimulPop(n,
                           oi,
                           model = "McFL",
                           initSize = 2000,
                           keepEvery = NA,
                           detectionProb = c(p2 = .1, n2 = 1500, checkSizePEvery = 1,
                                             PDBaseline = 1100, cPDetect = NA),
                           finalTime = NA, detectionSize = NA,
                           onlyCancer = TRUE,
                           detectionDrivers = NA, mc.cores = 2)
        (ta <- unlist(lapply(sa, function(x) x$FinalTime)))
        (tb <- unlist(lapply(sb, function(x) x$FinalTime)))
        print(suppressWarnings(wilcox.test(ta, tb, alternative = "greater")$p.value))
        T1 <- suppressWarnings(wilcox.test(ta, tb, alternative = "greater")$p.value < p.value.threshold)
        if(T1) break;
    }
    cat(paste("\n done tries", tries, "\n"))
    expect_true(T1)
})
date()



date()
test_that("Increasing cPDetect decreases time, Exp" , {
    gi <- rep(0.1, 10)
    names(gi) <- letters[1:10]
    oi <- allFitnessEffects(noIntGenes = gi)
    n <- 50
    max.tries <- 4  
    for(tries in 1:max.tries) {
        sa <- oncoSimulPop(n,
                           oi,
                           model = "Exp",
                           initSize = 1000,
                           keepEvery = NA,
                           detectionProb = c(p2 = NA, n2 = NA, checkSizePEvery = 5,
                                             PDBaseline = 500, cPDetect = 0.01), ## 1e-5),
                           finalTime = NA, detectionSize = NA,
                           onlyCancer = TRUE,
                           detectionDrivers = NA, mc.cores = 2)
        sb <- oncoSimulPop(n,
                           oi,
                           model = "Exp",
                           initSize = 1000,
                           keepEvery = NA,
                           detectionProb = c(p2 = NA, n2 = NA, checkSizePEvery = 5,
                                             PDBaseline = 500, cPDetect = 0.2), ## .01),
                           finalTime = NA, detectionSize = NA,
                           onlyCancer = TRUE,
                           detectionDrivers = NA, mc.cores = 2)
        ta <- unlist(lapply(sa, function(x) x$FinalTime))
        tb <- unlist(lapply(sb, function(x) x$FinalTime))
        print(suppressWarnings(wilcox.test(ta, tb, alternative = "greater")$p.value))
        T1 <- suppressWarnings(wilcox.test(ta, tb, alternative = "greater")$p.value < p.value.threshold)
        if(T1) break;
    }
    cat(paste("\n done tries", tries, "\n"))
    expect_true(T1)
})
date()

date()
test_that("Increasing p2 decreases time, Exp" , {
    gi <- rep(0.2, 10)
    names(gi) <- letters[1:10]
    oi <- allFitnessEffects(noIntGenes = gi)
    n <- 50
    max.tries <- 4  
    for(tries in 1:max.tries) {
        sa <- oncoSimulPop(n,
                           oi,
                           model = "Exp",
                           initSize = 2000,
                           keepEvery = NA,
                           detectionProb = c(p2 = .1, n2 = 8500, checkSizePEvery = 2,
                                             PDBaseline = 1100, cPDetect = NA),
                           finalTime = NA, detectionSize = NA,
                           onlyCancer = TRUE,
                           detectionDrivers = NA, mc.cores = 2)
        sb <- oncoSimulPop(n,
                           oi,
                           model = "Exp",
                           initSize = 2000,
                           keepEvery = NA,
                           detectionProb = c(p2 = .8, n2 = 8500, checkSizePEvery = 2,
                                             PDBaseline = 1100, cPDetect = NA),
                           finalTime = NA, detectionSize = NA,
                           onlyCancer = TRUE,
                           detectionDrivers = NA, mc.cores = 2)
        (ta <- unlist(lapply(sa, function(x) x$FinalTime)))
        (tb <- unlist(lapply(sb, function(x) x$FinalTime)))         
        print(suppressWarnings(wilcox.test(ta, tb, alternative = "greater")$p.value))
        T1 <- suppressWarnings(wilcox.test(ta, tb, alternative = "greater")$p.value < p.value.threshold)
        if(T1) break;
    }
    cat(paste("\n done tries", tries, "\n"))
    expect_true(T1)
})
date()


date()
test_that("Increasing n2 increases time, Exp" , {
    gi <- rep(0.1, 10)
    names(gi) <- letters[1:10]
    oi <- allFitnessEffects(noIntGenes = gi)
    n <- 50
    max.tries <- 4  
    for(tries in 1:max.tries) {
        sa <- oncoSimulPop(n,
                           oi,
                           model = "Exp",
                           initSize = 2000,
                           keepEvery = NA,
                           detectionProb = c(p2 = .21, n2 = 9000, PDBaseline = 2100, cPDetect = NA),
                           finalTime = NA, detectionSize = NA,
                           onlyCancer = TRUE,
                           detectionDrivers = NA, mc.cores = 2)
        sb <- oncoSimulPop(n,
                           oi,
                           model = "Exp",
                           initSize = 2000,
                           keepEvery = NA,
                           detectionProb = c(p2 = .21, n2 = 2101, PDBaseline = 2100, cPDetect = NA),
                           finalTime = NA, detectionSize = NA,
                           onlyCancer = TRUE,
                           detectionDrivers = NA, mc.cores = 2)
        (ta <- unlist(lapply(sa, function(x) x$FinalTime)))
        (tb <- unlist(lapply(sb, function(x) x$FinalTime)))
        print(suppressWarnings(wilcox.test(ta, tb, alternative = "greater")$p.value))
        T1 <- suppressWarnings(wilcox.test(ta, tb, alternative = "greater")$p.value < p.value.threshold)
        if(T1) break;
    }
    cat(paste("\n done tries", tries, "\n"))
    expect_true(T1)
})
date()



date()
test_that("Increasing checkSizePEvery increases time, Exp" , {
    gi <- rep(0.1, 10)
    names(gi) <- letters[1:10]
    oi <- allFitnessEffects(noIntGenes = gi)
    n <- 70
    max.tries <- 4  
    for(tries in 1:max.tries) {
        sa <- oncoSimulPop(n,
                           oi,
                           model = "Exp",
                           initSize = 2000,
                           keepEvery = NA,
                           detectionProb = c(p2 = .1, n2 = 1500, checkSizePEvery = 50,
                                             PDBaseline = 1100, cPDetect = NA),
                           finalTime = NA, detectionSize = NA,
                           onlyCancer = TRUE,
                           detectionDrivers = NA, mc.cores = 2)
        sb <- oncoSimulPop(n,
                           oi,
                           model = "Exp",
                           initSize = 2000,
                           keepEvery = NA,
                           detectionProb = c(p2 = .1, n2 = 1500, checkSizePEvery = 10,
                                             PDBaseline = 1100, cPDetect = NA),
                           finalTime = NA, detectionSize = NA,
                           onlyCancer = TRUE,
                           detectionDrivers = NA, mc.cores = 2)
        (ta <- unlist(lapply(sa, function(x) x$FinalTime)))
        (tb <- unlist(lapply(sb, function(x) x$FinalTime)))
        print(suppressWarnings(wilcox.test(ta, tb, alternative = "greater")$p.value))
        T1 <- suppressWarnings(wilcox.test(ta, tb, alternative = "greater")$p.value < p.value.threshold)
        if(T1) break;
    }
    cat(paste("\n done tries", tries, "\n"))
    expect_true(T1)
})
date()


## And there is no need for fitness effects
date()
test_that("Increasing cPDetect decreases time" , {
    gi <- rep(0.0,  10)
    names(gi) <- letters[1:10]
    oi <- allFitnessEffects(noIntGenes = gi)
    n <- 35 ## 75
    max.tries <- 4  
    for(tries in 1:max.tries) {
        sa <- oncoSimulPop(n,
                           oi,
                           model = "McFL",
                           initSize = 2000,
                           keepEvery = NA,
                           detectionProb = c(p2 = NA, n2 = NA, checkSizePEvery = 2,
                                             PDBaseline = 1100,  cPDetect = 0.1), ## 1e-4),
                           finalTime = NA, detectionSize = NA,
                           onlyCancer = FALSE,
                           detectionDrivers = NA, mc.cores = 2)
        sb <- oncoSimulPop(n,
                           oi,
                           model = "McFL",
                           initSize = 2000,
                           keepEvery = NA,
                           detectionProb = c(p2 = NA, n2 = NA,checkSizePEvery = 2,
                                             PDBaseline = 1100, cPDetect = 0.9), ## 1e-2),
                           finalTime = NA, detectionSize = NA,
                           onlyCancer = FALSE,
                           detectionDrivers = NA, mc.cores = 2)
        ta <- unlist(lapply(sa, function(x) x$FinalTime))
        tb <- unlist(lapply(sb, function(x) x$FinalTime))
        print(suppressWarnings(wilcox.test(ta, tb, alternative = "greater")$p.value))
        T1 <- suppressWarnings(wilcox.test(ta, tb, alternative = "greater")$p.value < p.value.threshold)
        if(T1) break;
    }
    cat(paste("\n done tries", tries, "\n"))
    expect_true(T1)
})
date()

date()
test_that("Increasing p2 decreases time" , {
    gi <- rep(0.0,  10)
    names(gi) <- letters[1:10]
    oi <- allFitnessEffects(noIntGenes = gi)
    n <- 40
    max.tries <- 4  
    for(tries in 1:max.tries) {
        sa <- oncoSimulPop(n,
                           oi,
                           model = "McFL",
                           initSize = 2000,
                           keepEvery = NA,
                           detectionProb = c(p2 = .1, n2 = 3500,checkSizePEvery = 2,
                                             PDBaseline = 1100, cPDetect = NA),
                           finalTime = NA, detectionSize = NA,
                           onlyCancer = FALSE,
                           detectionDrivers = NA, mc.cores = 2)
        sb <- oncoSimulPop(n,
                           oi,
                           model = "McFL",
                           initSize = 2000,
                           keepEvery = NA,
                           detectionProb = c(p2 = .6, n2 = 3500,checkSizePEvery = 2,
                                             PDBaseline = 1100, cPDetect = NA),
                           finalTime = NA, detectionSize = NA,
                           onlyCancer = FALSE,
                           detectionDrivers = NA, mc.cores = 2)
        (ta <- unlist(lapply(sa, function(x) x$FinalTime)))
        (tb <- unlist(lapply(sb, function(x) x$FinalTime)))
        print(suppressWarnings(wilcox.test(ta, tb, alternative = "greater")$p.value))
        T1 <- suppressWarnings(wilcox.test(ta, tb, alternative = "greater")$p.value < p.value.threshold)
        if(T1) break;
    }
    cat(paste("\n done tries", tries, "\n"))
    expect_true(T1)
})
date()


date()
test_that("Increasing n2 increases time" , {
    gi <- rep(0.0,  10)
    names(gi) <- letters[1:10]
    oi <- allFitnessEffects(noIntGenes = gi)
    n <- 40
    max.tries <- 4  
    for(tries in 1:max.tries) {
        sa <- oncoSimulPop(n,
                           oi,
                           model = "McFL",
                           initSize = 2000,
                           keepEvery = NA,
                           detectionProb = c(p2 = .15, n2 = 7000, checkSizePEvery = 1,
                                             PDBaseline = 1100, cPDetect = NA),
                           finalTime = NA, detectionSize = NA,
                           onlyCancer = FALSE,
                           detectionDrivers = NA, mc.cores = 2)
        sb <- oncoSimulPop(n,
                           oi,
                           model = "McFL",
                           initSize = 2000,
                           keepEvery = NA,
                           detectionProb = c(p2 = .15, n2 = 2001, checkSizePEvery = 1,
                                             PDBaseline = 1100, cPDetect = NA),
                           finalTime = NA, detectionSize = NA,
                           onlyCancer = FALSE,
                           detectionDrivers = NA, mc.cores = 2)
        (ta <- unlist(lapply(sa, function(x) x$FinalTime)))
        (tb <- unlist(lapply(sb, function(x) x$FinalTime)))
        print(suppressWarnings(wilcox.test(ta, tb, alternative = "greater")$p.value))
        T1 <- suppressWarnings(wilcox.test(ta, tb, alternative = "greater")$p.value < p.value.threshold)
        if(T1) break;
    }
    cat(paste("\n done tries", tries, "\n"))
    expect_true(T1)
})
date()



date()
test_that("Increasing checkSizePEvery increases time" , {
    gi <- rep(0.0,  10)
    names(gi) <- letters[1:10]
    oi <- allFitnessEffects(noIntGenes = gi)
    n <- 25
    max.tries <- 4  
    for(tries in 1:max.tries) {
        sa <- oncoSimulPop(n,
                           oi,
                           model = "McFL",
                           initSize = 2000,
                           keepEvery = NA,
                           detectionProb = c(p2 = .1, n2 = 1500, checkSizePEvery = 10,
                                             PDBaseline = 1100, cPDetect = NA),
                           finalTime = NA, detectionSize = NA,
                           onlyCancer = FALSE,
                           detectionDrivers = NA, mc.cores = 2)
        sb <- oncoSimulPop(n,
                           oi,
                           model = "McFL",
                           initSize = 2000,
                           keepEvery = NA,
                           detectionProb = c(p2 = .1, n2 = 1500, checkSizePEvery = 1,
                                             PDBaseline = 1100, cPDetect = NA),
                           finalTime = NA, detectionSize = NA,
                           onlyCancer = FALSE,
                           detectionDrivers = NA, mc.cores = 2)
        (ta <- unlist(lapply(sa, function(x) x$FinalTime)))
        (tb <- unlist(lapply(sb, function(x) x$FinalTime)))
        print(suppressWarnings(wilcox.test(ta, tb, alternative = "greater")$p.value))
        T1 <- suppressWarnings(wilcox.test(ta, tb, alternative = "greater")$p.value < p.value.threshold)
        if(T1) break;
    }
    cat(paste("\n done tries", tries, "\n"))
    expect_true(T1)
})
date()





date()
test_that("Exercise the default option and other substitutions/defaults" , {
    gi <- rep(0.1, 10)
    names(gi) <- letters[1:10]
    oi <- allFitnessEffects(noIntGenes = gi)
    expect_output(print(oncoSimulIndiv(
        oi,
        model = "Exp",
        initSize = 2000,
        keepEvery = NA,
        detectionProb = "default",
        finalTime = NA, detectionSize = NA,
        onlyCancer = TRUE,
        detectionDrivers = NA)),
        "Individual OncoSimul trajectory",
        fixed = TRUE)
    expect_output(print(oncoSimulIndiv(
        oi,
        model = "Exp",
        initSize = 2000,
        keepEvery = NA,
        detectionProb = c(cPDetect = 0.001),
        finalTime = NA, detectionSize = NA,
        onlyCancer = TRUE,
        detectionDrivers = NA)),
        "Individual OncoSimul trajectory",
        fixed = TRUE)
    expect_output(print(oncoSimulIndiv(
        oi,
        model = "Exp",
        initSize = 2000,
        keepEvery = NA,
        detectionProb = c(p2 = .9, n2 = 3000),
        finalTime = NA, detectionSize = NA,
        onlyCancer = TRUE,
        detectionDrivers = NA)),
        "Individual OncoSimul trajectory",
        fixed = TRUE)
    expect_output(print(oncoSimulIndiv(
        oi,
        model = "Exp",
        initSize = 2000,
        keepEvery = NA,
        extraTime = 2,
        detectionProb = c(p2 = .9, n2 = 3000),
        finalTime = NA, detectionSize = NA,
        onlyCancer = TRUE,
        detectionDrivers = NA)),
        "Individual OncoSimul trajectory",
        fixed = TRUE)
    expect_output(print(oncoSimulIndiv(
        oi,
        model = "Exp",
        initSize = 2000,
        keepEvery = NA,
        detectionProb = c(p2 = .9, n2 = 3000),
        finalTime = NA, detectionSize = NA,
        onlyCancer = TRUE,
        verbosity = 1,
        detectionDrivers = NA)),
        "Individual OncoSimul trajectory",
        fixed = TRUE)
    expect_output(print(oncoSimulIndiv(
        oi,
        model = "Exp",
        initSize = 2000,
        keepEvery = NA,
        detectionProb = c(PDBaseline = 2002),
        finalTime = NA, detectionSize = NA,
        onlyCancer = TRUE,
        verbosity = 1,
        detectionDrivers = NA)),
        "Individual OncoSimul trajectory",
        fixed = TRUE)
    expect_output(print(oncoSimulIndiv(
        oi,
        model = "Exp",
        initSize = 2000,
        keepEvery = NA,
        detectionProb = c(checkSizePEvery = 31),
        finalTime = NA, detectionSize = NA,
        onlyCancer = TRUE,
        verbosity = 1,
        detectionDrivers = NA)),
        "Individual OncoSimul trajectory",
        fixed = TRUE)
    expect_output(print(oncoSimulIndiv(
        oi,
        model = "Exp",
        initSize = 2000,
        keepEvery = NA,
        detectionProb = c(PDBaseline = 2002, checkSizePEvery = 17),
        finalTime = NA, detectionSize = NA,
        onlyCancer = TRUE,
        verbosity = 1,
        detectionDrivers = NA)),
        "Individual OncoSimul trajectory",
        fixed = TRUE)
    expect_output(print(oncoSimulIndiv(
        oi,
        model = "Exp",
        initSize = 2000,
        keepEvery = NA,
        detectionProb = c(n2 = 4000, p2 = 0.85, checkSizePEvery = 17),
        finalTime = NA, detectionSize = NA,
        onlyCancer = TRUE,
        verbosity = 1,
        detectionDrivers = NA)),
        "Individual OncoSimul trajectory",
        fixed = TRUE)
    expect_output(print(oncoSimulIndiv(
        oi,
        model = "Exp",
        initSize = 2000,
        keepEvery = NA,
        detectionProb = c(cPDetect = 0.001, checkSizePEvery = 17),
        finalTime = NA, detectionSize = NA,
        onlyCancer = TRUE,
        verbosity = 1,
        detectionDrivers = NA)),
        "Individual OncoSimul trajectory",
        fixed = TRUE)
    expect_output(print(oncoSimulIndiv(
        oi,
        model = "Exp",
        initSize = 2000,
        keepEvery = NA,
        detectionProb = c(cPDetect = 0.001, PDBaseline = 2030),
        finalTime = NA, detectionSize = NA,
        onlyCancer = TRUE,
        verbosity = 1,
        detectionDrivers = NA)),
        "Individual OncoSimul trajectory",
        fixed = TRUE)
    expect_output(print(oncoSimulIndiv(
        oi,
        model = "Exp",
        initSize = 2000, verbosity = -3,
        keepEvery = NA,
        detectionProb = c(cPDetect = 0.001),
        finalTime = NA, detectionSize = NA,
        onlyCancer = TRUE,
        detectionDrivers = NA)),
        "Individual OncoSimul trajectory",
        fixed = TRUE)
    expect_output(print(oncoSimulIndiv(
        oi,
        model = "Exp",
        initSize = 2000, verbosity = -3,
        keepEvery = NA,
        detectionProb = c(p2 = .9, n2 = 3000),
        finalTime = NA, detectionSize = NA,
        onlyCancer = TRUE,
        detectionDrivers = NA)),
        "Individual OncoSimul trajectory",
        fixed = TRUE)
    expect_output(print(oncoSimulIndiv(
        oi,
        model = "Exp",
        initSize = 2000, verbosity = -3,
        keepEvery = NA,
        detectionProb = c(PDBaseline = 2002),
        finalTime = NA, detectionSize = NA,
        onlyCancer = TRUE,
        detectionDrivers = NA)),
        "Individual OncoSimul trajectory",
        fixed = TRUE)
    expect_output(print(oncoSimulIndiv(
        oi,
        model = "Exp",
        initSize = 2000, verbosity = -3,
        keepEvery = NA,
        detectionProb = c(checkSizePEvery = 31),
        finalTime = NA, detectionSize = NA,
        onlyCancer = TRUE,
        detectionDrivers = NA)),
        "Individual OncoSimul trajectory",
        fixed = TRUE)
})
date()

date()
test_that("Fails as expected" , {
    data(examplePosets)
    p701 <- examplePosets[["p701"]]
    expect_error(oncoSimulIndiv(p701,
                                detectionProb = "default"),
                 "detectionProb cannot be used in v.1 objects",
                 fixed = TRUE)
    gi <- rep(0.1, 10)
    names(gi) <- letters[1:10]
    oi <- allFitnessEffects(noIntGenes = gi)
    expect_error(oncoSimulIndiv(oi,
                                model = "Exp",
                                initSize = 2000,
                                keepEvery = NA,
                                detectionProb = c(cPDete = 0.1, n2 = 3000, p2 = 0.9),
                                finalTime = NA, detectionSize = NA,
                                onlyCancer = TRUE,
                                detectionDrivers = NA),
                 "Names of some components of detectionProb are not recognized",
                 fixed = TRUE)
    gi <- rep(0.1, 10)
    names(gi) <- letters[1:10]
    oi <- allFitnessEffects(noIntGenes = gi)
    expect_error(oncoSimulIndiv(oi,
                                model = "Exp",
                                initSize = 2000,
                                keepEvery = NA,
                                detectionProb = c(cPDetect = 0.1, n2 = 3000, p2 = 0.9),
                                finalTime = NA, detectionSize = NA,
                                onlyCancer = TRUE,
                                detectionDrivers = NA),
                 "Specify only cPDetect",
                 fixed = TRUE)
    gi <- rep(0.1, 10)
    names(gi) <- letters[1:10]
    oi <- allFitnessEffects(noIntGenes = gi)
    expect_error(oncoSimulIndiv(oi,
                                model = "Exp",
                                initSize = 2000,
                                keepEvery = NA,
                                detectionProb = c(cPDetect = 0.1, n2 = 3000),
                                finalTime = NA, detectionSize = NA,
                                onlyCancer = TRUE,
                                detectionDrivers = NA),
                 "Specify only cPDetect",
                 fixed = TRUE)
    gi <- rep(0.1, 10)
    names(gi) <- letters[1:10]
    oi <- allFitnessEffects(noIntGenes = gi)
    expect_error(oncoSimulIndiv(oi,
                                model = "Exp",
                                initSize = 2000,
                                keepEvery = NA,
                                detectionProb = c(n2 = 3000),
                                finalTime = NA, detectionSize = NA,
                                onlyCancer = TRUE,
                                detectionDrivers = NA),
                 "If you pass one of n2 or p2, you must also",
                 fixed = TRUE)
    gi <- rep(0.1, 10)
    names(gi) <- letters[1:10]
    oi <- allFitnessEffects(noIntGenes = gi)
    expect_error(oncoSimulIndiv(oi,
                                model = "Exp",
                                initSize = 2000,
                                keepEvery = NA,
                                detectionProb = c(n2 = 3000, p2 = 1.1),
                                finalTime = NA, detectionSize = NA,
                                onlyCancer = TRUE,
                                detectionDrivers = NA),
                 "p2 >= 1",
                 fixed = TRUE)
    gi <- rep(0.1, 10)
    names(gi) <- letters[1:10]
    oi <- allFitnessEffects(noIntGenes = gi)
    expect_error(oncoSimulIndiv(oi,
                                model = "Exp",
                                initSize = 2000,
                                keepEvery = NA,
                                detectionProb = c(n2 = 3000, p2 = -.3),
                                finalTime = NA, detectionSize = NA,
                                onlyCancer = TRUE,
                                detectionDrivers = NA),
                 "p2 <= 0",
                 fixed = TRUE)
    gi <- rep(0.1, 10)
    names(gi) <- letters[1:10]
    oi <- allFitnessEffects(noIntGenes = gi)
    expect_error(oncoSimulIndiv(oi,
                                model = "Exp",
                                initSize = 2000,
                                keepEvery = NA,
                                detectionProb = c(n2 = 3000, p2 = .3, PDBaseline = 5000),
                                finalTime = NA, detectionSize = NA,
                                onlyCancer = TRUE,
                                detectionDrivers = NA),
                 "n2 <= PDBaseline",
                 fixed = TRUE)
    gi <- rep(0.1, 10)
    names(gi) <- letters[1:10]
    oi <- allFitnessEffects(noIntGenes = gi)
    expect_error(oncoSimulIndiv(oi,
                                model = "Exp",
                                initSize = 2000,
                                keepEvery = NA,
                                detectionProb = c(n2 = 3000, p2 = .3, PDBaseline = -3),
                                finalTime = NA, detectionSize = NA,
                                onlyCancer = TRUE,
                                detectionDrivers = NA),
                 "PDBaseline <= 0",
                 fixed = TRUE)
    expect_error(oncoSimulIndiv(oi,
                                model = "Exp",
                                initSize = 2000,
                                keepEvery = NA,
                                detectionProb = c(n2 = 3000, p2 = .3, PDBaseline = 0),
                                finalTime = NA, detectionSize = NA,
                                onlyCancer = TRUE,
                                detectionDrivers = NA),
                 "PDBaseline <= 0",
                 fixed = TRUE)
    expect_error(oncoSimulIndiv(oi,
                                model = "Exp",
                                initSize = 2000,
                                keepEvery = NA,
                                detectionProb = NA,
                                finalTime = NA,
                                detectionSize = NA,
                                onlyCancer = TRUE,
                                detectionDrivers = NA),
                 "At least one stopping condition should be given",
                 fixed = TRUE)
    expect_error(oncoSimulIndiv(oi,
                                model = "Exp",
                                initSize = 2000,
                                keepEvery = NA,
                                detectionProb = NA,
                                finalTime = NA,
                                detectionSize = NA,
                                onlyCancer = FALSE,
                                detectionDrivers = NA),
                 "At least one stopping condition should be given",
                 fixed = TRUE)
})
date()






    ## gi2 <- rep(0, 5)
    ## names(gi2) <- letters[1:5]
    ## oi2 <- allFitnessEffects(noIntGenes = gi2)
    ## ## nicely exponential, as expected
    ## s5 <- oncoSimulPop(100,
    ##                    oi2,
    ##                    model = "McFL",
    ##                    initSize = 1000,
    ##                    detectionSize = 4800,
    ##                    finalTime = NA, detectionSize = NA,
    ##                    keepEvery = -9,
    ##                    p2 = .1,
    ##                    checkSizePEvery = 2,
    ##                    verbosity = 0,
    ##                    PDBaseline = 1000,
    ##                    onlyCancer = TRUE,
    ##                    detectionDrivers = NA)
    ## s5
    ## hist(unlist(lapply(s5, function(x) x$FinalTime)))

cat(paste("\n Ending sample-prob tests", date(), "\n"))
cat(paste(" Took ", round(difftime(Sys.time(), inittime, units = "secs"), 2), "\n\n\n"))
rm(inittime)
