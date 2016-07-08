## Some tests below might only work on Linux because of compiler
## differences, because the rng is done in C++, etc.
## Note that the difference is in whether a certain code
## is exercised. The runs should work in all platforms, though.

test_that("exercise no positions left for mutation, updating in null mut, old format", {
    ## Do not do the capture output from oncoSimulPop,
    ## as that comes from mclapply and it is a mess.
    RNGkind("Mersenne-Twister")
    set.seed(1)
    p1 <- cbind(1L, 2L)
    st <- capture.output({
        pp1 <- oncoSimulIndiv(p1,
                              sh = 0,
                              initSize = 1e5,
                              sampleEvery = 0.02,
                              detectionSize = 1e9,
                              model = "Exp",
                              finalTime = 2000,
                              extraTime = 3.17,
                              onlyCancer = FALSE,
                              seed = NULL)
    })
    if(Sys.info()["sysname"] == "Linux") {
        expect_true(any(grepl("mutation = 0", st)))
        expect_true(any(grepl("updating in null mutation", st)))
    }
    expect_output(print(pp1),
                   "Individual OncoSimul trajectory", fixed = TRUE)
})


test_that("exercise mu > 1, old format", {
    RNGkind("Mersenne-Twister")
    set.seed(2)
    data(examplePosets)
    p701 <- examplePosets[["p701"]]
    st <- capture.output(pp1 <- oncoSimulIndiv(p701,
                          mu = 0.7,
                          sh = 0,
                          initSize = 1e5,
                          sampleEvery = 0.02,
                          detectionSize = 1e6,
                          model = "Exp",
                          finalTime = 2000,
                          onlyCancer = FALSE,
                          seed = NULL))
    
    expect_true(any(grepl("mutation > 1", st)))
    
    expect_output(print(pp1),
                  "Individual OncoSimul trajectory", fixed = TRUE)
})

test_that("using old poset format, hitting wall time", {
    RNGkind("Mersenne-Twister")
    set.seed(1)
    data(examplePosets)
    p701 <- examplePosets[["p701"]]
    pet <- oncoSimulIndiv(p701, sh = 0,
                          initSize = 1e2,
                          detectionSize = 5e8,
                          model = "McFL",
                          finalTime = 1e6,
                          extraTime = 3.17,
                          max.wall.time = 0.0001,
                          onlyCancer = FALSE,
                          seed = NULL)
    expect_true(pet$HittedWallTime)
})


test_that("using old poset format, verbose exercise iteration", {
    RNGkind("Mersenne-Twister")
    set.seed(1)
    data(examplePosets)
    p701 <- examplePosets[["p701"]]
    ## st <- capture.output(
        p1 <- oncoSimulIndiv(p701, sh = 0,
                             initSize = 1e4,
                             detectionSize = 1e7,
                             model = "Exp",
                             finalTime = 1994,
                             extraTime = 3.17,
                             onlyCancer = FALSE,
                             verbosity = 2,
                             seed = NULL)
    ## )
    expect_output(print(p1), "Individual OncoSimul", fixed = TRUE)
})



test_that("using old poset format, exercising minDetectDrv", {
    RNGkind("Mersenne-Twister")
    set.seed(1)
    data(examplePosets)
    p701 <- examplePosets[["p701"]]
    p1 <- oncoSimulIndiv(p701, sh = 0,
                         initSize = 1e2,
                         detectionSize = 5e8,
                         model = "McFL",
                         finalTime = 1e6,
                         extraTime = 233.17,
                         minDetectDrvCloneSz = 148,
                         detectionDrivers = 0,
                         onlyCancer = TRUE,
                         seed = NULL)
    expect_true(max(p1$pops.by.time[, -1]) >= 148)
    expect_silent(p1)
})

test_that("exercise no positions left for mutation, updating in null mut, new format", {
    RNGkind("Mersenne-Twister")
    set.seed(2)
    gg <- rep(0.01, 2)
    names(gg) <- letters[1:2]
    ii <- allFitnessEffects(noIntGenes = gg)
    st <- capture.output(pp1 <- oncoSimulIndiv(ii,
                          initSize = 5e6,
                          sampleEvery = 0.02,
                          detectionSize = 5e8,
                          model = "Exp",
                          keepEvery = NA,
                          finalTime = 2000,
                          extraTime = 3.17,
                          onlyCancer = FALSE,
                          seed = NULL,
                          verbosity = 1))
    expect_output(print(pp1),
                  "Individual OncoSimul", fixed = TRUE)
    if(Sys.info()["sysname"] == "Linux") {
        expect_true(any(grepl("mutation = 0", st)))
        expect_true(any(grepl("updating in null mutation", st)))
    }
})


test_that("exercise mu > 1, new format", {
    RNGkind("Mersenne-Twister")
    set.seed(2)
    gg <- rep(0.01, 3)
    names(gg) <- letters[1:3]
    ii <- allFitnessEffects(noIntGenes = gg)
    st <- capture.output(
        pp1 <- oncoSimulIndiv(ii,
                              mu = 0.7,
                              sh = 0,
                              initSize = 1e3,
                              sampleEvery = 0.02,
                              detectionSize = 1.1e3,
                              model = "Exp",
                              finalTime = 2000,
                              onlyCancer = FALSE,
                              seed = NULL))
    expect_true(any(grepl("mutation > 1", st)))
    expect_output(print(pp1),
                  "Individual OncoSimul trajectory", fixed = TRUE)
})



set.seed(NULL)
