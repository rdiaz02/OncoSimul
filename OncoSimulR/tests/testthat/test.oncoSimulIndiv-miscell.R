inittime <- Sys.time()
cat(paste("\n Starting oncoSimulIndiv-miscell tests", date(), "\n"))

test_that("sampleEvery must have a value", {
     oi <- allFitnessEffects(orderEffects =
                                c("F > D" = -0.3, "D > F" = 0.4),
                            noIntGenes = rexp(5, 10),
                            geneToModule =
                                c("Root" = "Root",
                                  "F" = "f1, f2, f3",
                                  "D" = "d1, d2") )
     expect_error(oncoSimulIndiv(pi, sampleEvery = NA),
                  "sampleEvery cannot be NULL or NA",
                  fixed = TRUE)
})

## RNGkind("Mersenne-Twister")
test_that("can start from 1 individual but error if McFL", {
    oi <- allFitnessEffects(orderEffects =
                                c("F > D" = -0.3, "D > F" = 0.4),
                            noIntGenes = rexp(5, 10),
                            geneToModule =
                                c("Root" = "Root",
                                  "F" = "f1, f2, f3",
                                  "D" = "d1, d2") )

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

test_that("oncoSimulSample and oncoSimulPop require >= 1 indiv", {
    pancr <- allFitnessEffects(data.frame(parent = c("Root", rep("KRAS", 4), "SMAD4", "CDNK2A", 
                                                        "TP53", "TP53", "MLL3"),
                                                    child = c("KRAS","SMAD4", "CDNK2A", 
                                                        "TP53", "MLL3",
                                                        rep("PXDN", 3), rep("TGFBR2", 2)),
                                                    s = 0.05,
                                                    sh = -0.3,
                                                    typeDep = "MN"))
    expect_error(pS <- oncoSimulSample(0, pancr),
                 "Nindiv must be >= 1",
                 fixed = TRUE)
    expect_error(pS <- oncoSimulPop(0, pancr),
                 "Nindiv must be >= 1",
                 fixed = TRUE)
    expect_error(pS <- oncoSimulSample(-3, pancr),
                 "Nindiv must be >= 1",
                 fixed = TRUE)
    expect_error(pS <- oncoSimulPop(-4, pancr),
                 "Nindiv must be >= 1",
                 fixed = TRUE)
})

test_that("numPassengers no effect with fitnessEffects objects", {
    oi <- allFitnessEffects(orderEffects =
                                c("F > D" = -0.3, "D > F" = 0.4),
                            noIntGenes = rexp(5, 10),
                            geneToModule =
                                c("Root" = "Root",
                                  "F" = "f1, f2, f3",
                                  "D" = "d1, d2") )
    expect_warning(oi1 <- oncoSimulIndiv(oi,
                                         sampleEvery = 0.03,
                                         keepEvery = 5,
                                         numPassengers = 10),
                   "Specifying numPassengers", fixed = TRUE)
})

test_that("initSize cannot be less than 1", {
    np <- 20
    s <- 0.015
    spp <- 0.01
    nd <- 5
    mcf1 <- allFitnessEffects(noIntGenes = rep(spp, np),
                              drvNames = integer(0))
    expect_error(oncoSimulIndiv(mcf1,
                                initSize = 0),
                 "initSize < 1", fixed = TRUE)
    expect_error(oncoSimulIndiv(mcf1,
                                initSize = -3),
                 "initSize < 1", fixed = TRUE)
})

test_that("samplePop with oncoSimulIndiv object", {
              oi <- allFitnessEffects(orderEffects =
               c("F > D" = -0.3, "D > F" = 0.4),
               noIntGenes = rexp(5, 10),
                          geneToModule =
                              c("Root" = "Root",
                                "F" = "f1, f2, f3",
                                "D" = "d1, d2") )
              oiI1 <- oncoSimulIndiv(oi,
                                     sampleEvery = 0.03,
                                     keepEvery = 5,
                                     model = "Exp")
              expect_message(out <- samplePop(oiI1),
                             "Subjects by Genes matrix of 1 subjects and 10 genes.",
                             fixed = TRUE)
              expect_true(ncol(out) == 10)
              expect_true(nrow(out) == 1)
          })


test_that("McFl warning with small sampleEvery", {
              oi <- allFitnessEffects(orderEffects =
               c("F > D" = -0.3, "D > F" = 0.4),
               noIntGenes = rexp(5, 10),
                          geneToModule =
                              c("Root" = "Root",
                                "F" = "f1, f2, f3",
                                "D" = "d1, d2") )
              expect_warning(oiI1 <- oncoSimulIndiv(oi, model = "McFL",
                                                    sampleEvery = 0.06,
                                                    keepEvery = 5,
                                                    onlyCancer = FALSE),
                             "With the McFarland model you often want smaller sampleEvery")
          })

test_that("mu < 0 error", {
              oi <- allFitnessEffects(orderEffects =
               c("F > D" = -0.3, "D > F" = 0.4),
               noIntGenes = rexp(5, 10),
                          geneToModule =
                              c("Root" = "Root",
                                "F" = "f1, f2, f3",
                                "D" = "d1, d2") )
              expect_error(oiI1 <- oncoSimulIndiv(oi, mu = -1),
                           "mutation rate (mu) is negative",
                           fixed = TRUE)                                                     
          })


test_that("keepEevery and sampleEvery consistency", {
              oi <- allFitnessEffects(orderEffects =
               c("F > D" = -0.3, "D > F" = 0.4),
               noIntGenes = rexp(5, 10),
                          geneToModule =
                              c("Root" = "Root",
                                "F" = "f1, f2, f3",
                                "D" = "d1, d2") )
              expect_warning(oiI1 <- oncoSimulIndiv(oi,
                                                    finalTime = 0.5,
                                                    keepEvery = 0.25,
                                                    onlyCancer = FALSE,
                                                    sampleEvery = 0.3),
                             "setting keepEvery <- sampleEvery",
                             fixed = TRUE)                                                     
          })


test_that("using old poset format: error if need param null", {
              data(examplePosets)
              p701 <- examplePosets[["p701"]]
              expect_error(oncoSimulIndiv(p701, s = NULL),
                           "You must specify all of")                                                     
          })

test_that("using old poset format: error if need param null", {
              data(examplePosets)
              p701 <- examplePosets[["p701"]]
              expect_error(oncoSimulIndiv(p701, sh = NULL),
                           "You must specify all of")                                                     
          })

test_that("using old poset format: error if need param null", {
              data(examplePosets)
              p701 <- examplePosets[["p701"]]
              expect_error(oncoSimulIndiv(p701, numPassengers = NULL),
                           "You must specify all of")                                                     
          })

test_that("verbosity options", {
              oi <- allFitnessEffects(orderEffects =
                                          c("F > D" = -0.3, "D > F" = 0.4),
                                      noIntGenes = rexp(5, 10),
                                      geneToModule =
                                          c("Root" = "Root",
                                            "F" = "f1, f2, f3",
                                            "D" = "d1, d2") )
              expect_output(oncoSimulIndiv(oi, verbosity = 2,
                                           detectionSize = 1e3,
                                           sampleEvery = 0.03,
                                           keepEvery = 2,
                                           finalTime = 3,
                                           onlyCancer = FALSE),
                            "Total Pop Size = ", fixed = TRUE)
          })


test_that("printing oncosimul object", {
              oi <- allFitnessEffects(orderEffects =
                                          c("F > D" = -0.3, "D > F" = 0.4),
                                      noIntGenes = rexp(5, 10),
                                      geneToModule =
                                          c("Root" = "Root",
                                            "F" = "f1, f2, f3",
                                            "D" = "d1, d2") )
              out <- oncoSimulIndiv(oi,
                                    sampleEvery = 0.03,
                                    keepEvery = 2,
                                    detectionSize = 1e3,
                                    finalTime = 3,
                                    onlyCancer = FALSE)
              expect_output(print(out),
                            "Individual OncoSimul trajectory with call")
          })


test_that("printing oncosimul pop object", {
              oi <- allFitnessEffects(orderEffects =
                                          c("F > D" = -0.3, "D > F" = 0.4),
                                      noIntGenes = rexp(5, 10),
                                      geneToModule =
                                          c("Root" = "Root",
                                            "F" = "f1, f2, f3",
                                            "D" = "d1, d2") )
              out <- oncoSimulPop(4,
                                  oi,
                                  sampleEvery = 0.03,
                                  keepEvery = 2,
                                  finalTime = 3,
                                  detectionSize = 1e3,
                                  onlyCancer = FALSE, mc.cores = 2)
              expect_output(print(out),
                            "Population of OncoSimul trajectories of 4 individuals")
          })



test_that("exercising oncoSimulSample, old format", {
              data(examplePosets)
              p701 <- examplePosets[["p701"]]
              expect_message(ofw <- oncoSimulSample(2, p701,
                                                    sampleEvery = 0.03,
                                                    detectionSize = 1e3,
                                                    finalTime = 3,
                                                    onlyCancer = FALSE,
                                                    showProgress = TRUE),
                             "Successfully sampled 2 individuals")
              expect_message(ofs <- oncoSimulSample(2, p701,
                                                    sampleEvery = 0.03,
                                                    detectionSize = 1e3,
                                                    finalTime = 3,
                                                    onlyCancer = FALSE,
                                                    typeSample = "single"),
                             "Successfully sampled 2 individuals")
              expect_equal(dim(ofw$popSample), c(2, 7))
              expect_equal(dim(ofs$popSample), c(2, 7))
          })


test_that("exercising oncoSimulSample, new format", {
              pancr <- allFitnessEffects(data.frame(parent = c("Root", rep("KRAS", 4), "SMAD4", "CDNK2A", 
                                                        "TP53", "TP53", "MLL3"),
                                                    child = c("KRAS","SMAD4", "CDNK2A", 
                                                        "TP53", "MLL3",
                                                        rep("PXDN", 3), rep("TGFBR2", 2)),
                                                    s = 0.05,
                                                    sh = -0.3,
                                                    typeDep = "MN"))
              expect_message(pS <- oncoSimulSample(2, pancr, sampleEvery = 0.03,
                                                    detectionSize = 1e3,
                                                    finalTime = 3,
                                                    onlyCancer = FALSE),
                             "Successfully sampled 2 individuals")
              expect_message(
                  pSs <- oncoSimulSample(2,
                                         pancr,
                                         sampleEvery = 0.03,
                                         detectionSize = 1e3,
                                         finalTime = 3,
                                         onlyCancer = FALSE,
                                         typeSample = "single",
                                         showProgress = TRUE),
                  "Successfully sampled 2 individuals")
              expect_equal(dim(pS$popSample), c(2, 7))
              expect_equal(dim(pSs$popSample), c(2, 7))
          })



test_that("check error unknown timeSample", {
    data(examplePosets)
    p701 <- examplePosets[["p701"]]
    r1 <- oncoSimulIndiv(p701, onlyCancer = TRUE, max.num.tries = 5000)
    expect_error(samplePop(r1, timeSample = "uniformo"), 
                 "Unknown timeSample option")
    expect_error(samplePop(r1, timeSample = "uni"), 
                 "Unknown timeSample option")
    expect_error(samplePop(r1, timeSample = "lasto"), 
                 "Unknown timeSample option")
    expect_error(samplePop(r1, timeSample = "whole"), 
                 "Unknown timeSample option")
    expect_error(samplePop(r1, timeSample = "single"), 
                 "Unknown timeSample option")
    expect_error(samplePop(r1, timeSample = "cucu"), 
                 "Unknown timeSample option")
})

test_that("check error unknown typeSample", {
    data(examplePosets)
    p701 <- examplePosets[["p701"]]
    r1 <- oncoSimulIndiv(p701, onlyCancer = TRUE, max.num.tries = 5000)
    expect_error(samplePop(r1, typeSample = "uniformo"), 
                 "Unknown typeSample option")
    expect_error(samplePop(r1, typeSample = "uni"), 
                 "Unknown typeSample option")
    expect_error(samplePop(r1, typeSample = "lasto"), 
                 "Unknown typeSample option")
    expect_error(samplePop(r1, typeSample = "uniform"), 
                 "Unknown typeSample option")
    expect_error(samplePop(r1, typeSample = "last"), 
                 "Unknown typeSample option")
    expect_error(samplePop(r1, typeSample = "wholo"), 
                 "Unknown typeSample option")
    expect_error(samplePop(r1, typeSample = "whola"), 
                 "Unknown typeSample option")
    expect_error(samplePop(r1, typeSample = "singlo"), 
                 "Unknown typeSample option")
    expect_error(samplePop(r1, typeSample = "cell"), 
                 "Unknown typeSample option")
    expect_error(samplePop(r1, typeSample = "cucu"), 
                 "Unknown typeSample option")
})


test_that("oncosimul sample without drivers", {
    np <- 20
    s <- 0.015
    spp <- 0.01
    nd <- 5
    mcf1 <- allFitnessEffects(noIntGenes = rep(spp, np),
                              drvNames = integer(0))
    mcf2 <- allFitnessEffects(noIntGenes = c(rep(s, nd), rep(spp, np)),
                              drvNames = character(0))
    expect_message(o1 <- oncoSimulSample(2, mcf1, sampleEvery = 0.03,
                                                    detectionSize = 1e3,
                                                    finalTime = 3,
                                                    onlyCancer = FALSE),
                   "Successfully sampled 2 individuals")
    expect_message(o2 <- oncoSimulSample(5, mcf2, sampleEvery = 0.03,
                                                    detectionSize = 1e3,
                                                    finalTime = 3,
                                                    onlyCancer = FALSE),
                   "Successfully sampled 5 individuals")
})


test_that("samplePop, a few examples with time last and whole pop", {
    ## Why am I testing this? No particular reason. But once I got
    ## confused and thought this function, which is the basis of
    ## oncoSimulSample, had a bug. The bug was in my understanding.
    initSize <- 10 ## as it is passed below, in the structure; the rest
                   ## for the same reason
    K <- 1
    mu <- 1e-6
    sampleEvery <- 0.1
    ni <- c(3, 1, rep(0, 5))
    names(ni) <- c("a", "b", paste0("n", seq.int(5)))
    fe <- allFitnessEffects(noIntGenes = ni)
    fp <- fe
    numPassengers <- 0
    seed <- NULL
    verbosity <- 0
    initMutant <- NULL
    finalTime <- 3.0
    max.memory <- 2000
    mutationPropGrowth <- TRUE
    max.wall.time.total <- 600
    possibleAttempts <- max.num.tries.total <- 500 * 10
    onlyCancer <- FALSE
    keepPhylog <- FALSE
    s1 <- structure(list(structure(list(pops.by.time = structure(c(3, 62925, 
2, 4, 1, 57), .Dim = c(1L, 6L)), NumClones = 5, TotalPopSize = 62989, 
    Genotypes = structure(c(1L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 1L, 
    0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 1L, 0L, 0L, 1L, 0L, 0L, 
    0L, 0L, 1L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 1L), .Dim = c(7L, 
    5L)), GenotypesWDistinctOrderEff = list(1L, 1:2, c(1L, 5L
    ), c(1L, 6L), c(1L, 7L)), GenotypesLabels = c("a", "a, b", 
    "a, n3", "a, n4", "a, n5"), MaxNumDrivers = 0L, MaxDriversLast = 0L, 
    NumDriversLargestPop = 0L, LargestClone = 62925, PropLargestPopLast = 0.998983949578498, 
    FinalTime = 3, NumIter = 34L, HittedWallTime = FALSE, HittedMaxTries = FALSE, 
    TotalPresentDrivers = 0L, CountByDriver = integer(0), OccurringDrivers = "", 
    PerSampleStats = structure(c(62989, 62925, 0.998983949578498, 
    0, 0), .Dim = c(1L, 5L)), other = structure(list(attemptsUsed = 1L, 
        errorMF = -99, errorMF_size = 0,  minDMratio = 25000, 
        minBMratio = 166666.666666667,  PhylogDF = structure(list(
            parent = structure(integer(0), .Label = character(0), class = "factor"), 
            child = structure(integer(0), .Label = character(0), class = "factor"), 
            time = numeric(0)), .Names = c("parent", "child", 
        "time"), row.names = integer(0), class = "data.frame"), 
        UnrecoverExcept = FALSE), .Names = c("attemptsUsed", 
    "errorMF", "errorMF_size",  "minDMratio", "minBMratio", 
     "PhylogDF", "UnrecoverExcept")), Drivers = integer(0), 
    geneNames = c("a", "b", "n1", "n2", "n3", "n4", "n5")), .Names = c("pops.by.time", 
"NumClones", "TotalPopSize", "Genotypes", "GenotypesWDistinctOrderEff", 
"GenotypesLabels", "MaxNumDrivers", "MaxDriversLast", "NumDriversLargestPop", 
"LargestClone", "PropLargestPopLast", "FinalTime", "NumIter", 
"HittedWallTime", "HittedMaxTries", "TotalPresentDrivers", "CountByDriver", 
"OccurringDrivers", "PerSampleStats", "other", "Drivers", "geneNames"
), class = c("oncosimul", "oncosimul2"), call = oncoSimulIndiv(fp = fp, 
    model = "Exp", numPassengers = numPassengers, mu = mu, detectionSize = 1e9, detectionDrivers = 99, 
    sampleEvery = sampleEvery, initSize = initSize, s = s, sh = sh, 
    K = K, keepEvery = -9, minDetectDrvCloneSz = "auto", extraTime = 0, 
    finalTime = finalTime, onlyCancer = onlyCancer, keepPhylog = keepPhylog, 
    mutationPropGrowth = mutationPropGrowth, max.memory = max.memory, 
    max.wall.time = max.wall.time.total, max.num.tries = possibleAttempts, 
    errorHitWallTime = TRUE, errorHitMaxTries = TRUE, verbosity = verbosity, 
    initMutant = initMutant, seed = seed)), structure(list(pops.by.time = structure(c(3, 
16000), .Dim = 1:2), NumClones = 1, TotalPopSize = 16000, Genotypes = structure(c(1L, 
0L, 0L, 0L, 0L, 0L, 0L), .Dim = c(7L, 1L)), GenotypesWDistinctOrderEff = list(
    1L), GenotypesLabels = "a", MaxNumDrivers = 0L, MaxDriversLast = 0L, 
    NumDriversLargestPop = 0L, LargestClone = 16000, PropLargestPopLast = 1, 
    FinalTime = 3, NumIter = 30L, HittedWallTime = FALSE, HittedMaxTries = FALSE, 
    TotalPresentDrivers = 0L, CountByDriver = integer(0), OccurringDrivers = "", 
    PerSampleStats = structure(c(16000, 16000, 1, 0, 0), .Dim = c(1L, 
    5L)), other = structure(list(attemptsUsed = 1L, errorMF = -99, 
        errorMF_size = 0,  minDMratio = 41666.6666666667, 
        minBMratio = 166666.666666667,  PhylogDF = structure(list(
            parent = structure(integer(0), .Label = character(0), class = "factor"), 
            child = structure(integer(0), .Label = character(0), class = "factor"), 
            time = numeric(0)), .Names = c("parent", "child", 
        "time"), row.names = integer(0), class = "data.frame"), 
        UnrecoverExcept = FALSE), .Names = c("attemptsUsed", 
    "errorMF", "errorMF_size",  "minDMratio", "minBMratio", 
     "PhylogDF", "UnrecoverExcept")), Drivers = integer(0), 
    geneNames = c("a", "b", "n1", "n2", "n3", "n4", "n5")), .Names = c("pops.by.time", 
"NumClones", "TotalPopSize", "Genotypes", "GenotypesWDistinctOrderEff", 
"GenotypesLabels", "MaxNumDrivers", "MaxDriversLast", "NumDriversLargestPop", 
"LargestClone", "PropLargestPopLast", "FinalTime", "NumIter", 
"HittedWallTime", "HittedMaxTries", "TotalPresentDrivers", "CountByDriver", 
"OccurringDrivers", "PerSampleStats", "other", "Drivers", "geneNames"
), class = c("oncosimul", "oncosimul2"), call = oncoSimulIndiv(fp = fp, 
    model = "Exp", numPassengers = numPassengers, mu = mu, detectionSize = 1e9, detectionDrivers = 99, 
    sampleEvery = sampleEvery, initSize = initSize, s = s, sh = sh, 
    K = K, keepEvery = -9, minDetectDrvCloneSz = "auto", extraTime = 0, 
    finalTime = finalTime, onlyCancer = onlyCancer, keepPhylog = keepPhylog, 
    mutationPropGrowth = mutationPropGrowth, max.memory = max.memory, 
    max.wall.time = max.wall.time.total, max.num.tries = possibleAttempts, 
    errorHitWallTime = TRUE, errorHitMaxTries = TRUE, verbosity = verbosity, 
    initMutant = initMutant, seed = seed)), structure(list(pops.by.time = structure(c(3, 
65810, 4), .Dim = c(1L, 3L)), NumClones = 2, TotalPopSize = 65814, 
    Genotypes = structure(c(1L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 
    0L, 0L, 1L, 0L, 0L), .Dim = c(7L, 2L)), GenotypesWDistinctOrderEff = list(
        1L, c(1L, 5L)), GenotypesLabels = c("a", "a, n3"), MaxNumDrivers = 0L, 
    MaxDriversLast = 0L, NumDriversLargestPop = 0L, LargestClone = 65810, 
    PropLargestPopLast = 0.999939222657793, FinalTime = 3, NumIter = 31L, 
    HittedWallTime = FALSE, HittedMaxTries = FALSE, TotalPresentDrivers = 0L, 
    CountByDriver = integer(0), OccurringDrivers = "", PerSampleStats = structure(c(65814, 
    65810, 0.999939222657793, 0, 0), .Dim = c(1L, 5L)), other = structure(list(
        attemptsUsed = 1L, errorMF = -99, errorMF_size = 0,  
        minDMratio = 41666.6666666667, minBMratio = 166666.666666667, 
         PhylogDF = structure(list(parent = structure(integer(0), .Label = character(0), class = "factor"), 
            child = structure(integer(0), .Label = character(0), class = "factor"), 
            time = numeric(0)), .Names = c("parent", "child", 
        "time"), row.names = integer(0), class = "data.frame"), 
        UnrecoverExcept = FALSE), .Names = c("attemptsUsed", 
    "errorMF", "errorMF_size",  "minDMratio", "minBMratio", 
     "PhylogDF", "UnrecoverExcept")), Drivers = integer(0), 
    geneNames = c("a", "b", "n1", "n2", "n3", "n4", "n5")), .Names = c("pops.by.time", 
"NumClones", "TotalPopSize", "Genotypes", "GenotypesWDistinctOrderEff", 
"GenotypesLabels", "MaxNumDrivers", "MaxDriversLast", "NumDriversLargestPop", 
"LargestClone", "PropLargestPopLast", "FinalTime", "NumIter", 
"HittedWallTime", "HittedMaxTries", "TotalPresentDrivers", "CountByDriver", 
"OccurringDrivers", "PerSampleStats", "other", "Drivers", "geneNames"
), class = c("oncosimul", "oncosimul2"), call = oncoSimulIndiv(fp = fp, 
    model = "Exp", numPassengers = numPassengers, mu = mu, detectionSize = 1e9, detectionDrivers = 99, 
    sampleEvery = sampleEvery, initSize = initSize, s = s, sh = sh, 
    K = K, keepEvery = -9, minDetectDrvCloneSz = "auto", extraTime = 0, 
    finalTime = finalTime, onlyCancer = onlyCancer, keepPhylog = keepPhylog, 
    mutationPropGrowth = mutationPropGrowth, max.memory = max.memory, 
    max.wall.time = max.wall.time.total, max.num.tries = possibleAttempts, 
    errorHitWallTime = TRUE, errorHitMaxTries = TRUE, verbosity = verbosity, 
    initMutant = initMutant, seed = seed)), structure(list(pops.by.time = structure(c(3, 
51024), .Dim = 1:2), NumClones = 1, TotalPopSize = 51024, Genotypes = structure(c(1L, 
0L, 0L, 0L, 0L, 0L, 0L), .Dim = c(7L, 1L)), GenotypesWDistinctOrderEff = list(
    1L), GenotypesLabels = "a", MaxNumDrivers = 0L, MaxDriversLast = 0L, 
    NumDriversLargestPop = 0L, LargestClone = 51024, PropLargestPopLast = 1, 
    FinalTime = 3, NumIter = 30L, HittedWallTime = FALSE, HittedMaxTries = FALSE, 
    TotalPresentDrivers = 0L, CountByDriver = integer(0), OccurringDrivers = "", 
    PerSampleStats = structure(c(51024, 51024, 1, 0, 0), .Dim = c(1L, 
    5L)), other = structure(list(attemptsUsed = 1L, errorMF = -99, 
        errorMF_size = 0,  minDMratio = 41666.6666666667, 
        minBMratio = 166666.666666667,  PhylogDF = structure(list(
            parent = structure(integer(0), .Label = character(0), class = "factor"), 
            child = structure(integer(0), .Label = character(0), class = "factor"), 
            time = numeric(0)), .Names = c("parent", "child", 
        "time"), row.names = integer(0), class = "data.frame"), 
        UnrecoverExcept = FALSE), .Names = c("attemptsUsed", 
    "errorMF", "errorMF_size",  "minDMratio", "minBMratio", 
     "PhylogDF", "UnrecoverExcept")), Drivers = integer(0), 
    geneNames = c("a", "b", "n1", "n2", "n3", "n4", "n5")), .Names = c("pops.by.time", 
"NumClones", "TotalPopSize", "Genotypes", "GenotypesWDistinctOrderEff", 
"GenotypesLabels", "MaxNumDrivers", "MaxDriversLast", "NumDriversLargestPop", 
"LargestClone", "PropLargestPopLast", "FinalTime", "NumIter", 
"HittedWallTime", "HittedMaxTries", "TotalPresentDrivers", "CountByDriver", 
"OccurringDrivers", "PerSampleStats", "other", "Drivers", "geneNames"
), class = c("oncosimul", "oncosimul2"), call = oncoSimulIndiv(fp = fp, 
    model = "Exp", numPassengers = numPassengers, mu = mu, detectionSize = 1e9, detectionDrivers = 99, 
    sampleEvery = sampleEvery, initSize = initSize, s = s, sh = sh, 
    K = K, keepEvery = -9, minDetectDrvCloneSz = "auto", extraTime = 0, 
    finalTime = finalTime, onlyCancer = onlyCancer, keepPhylog = keepPhylog, 
    mutationPropGrowth = mutationPropGrowth, max.memory = max.memory, 
    max.wall.time = max.wall.time.total, max.num.tries = possibleAttempts, 
    errorHitWallTime = TRUE, errorHitMaxTries = TRUE, verbosity = verbosity, 
    initMutant = initMutant, seed = seed)), structure(list(pops.by.time = structure(c(3, 
64375, 1), .Dim = c(1L, 3L)), NumClones = 2, TotalPopSize = 64376, 
    Genotypes = structure(c(1L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 
    1L, 0L, 0L, 0L, 0L), .Dim = c(7L, 2L)), GenotypesWDistinctOrderEff = list(
        1L, c(1L, 3L)), GenotypesLabels = c("a", "a, n1"), MaxNumDrivers = 0L, 
    MaxDriversLast = 0L, NumDriversLargestPop = 0L, LargestClone = 64375, 
    PropLargestPopLast = 0.999984466260718, FinalTime = 3, NumIter = 31L, 
    HittedWallTime = FALSE, HittedMaxTries = FALSE, TotalPresentDrivers = 0L, 
    CountByDriver = integer(0), OccurringDrivers = "", PerSampleStats = structure(c(64376, 
    64375, 0.999984466260718, 0, 0), .Dim = c(1L, 5L)), other = structure(list(
        attemptsUsed = 1L, errorMF = -99, errorMF_size = 0,  
        minDMratio = 41666.6666666667, minBMratio = 166666.666666667, 
         PhylogDF = structure(list(parent = structure(integer(0), .Label = character(0), class = "factor"), 
            child = structure(integer(0), .Label = character(0), class = "factor"), 
            time = numeric(0)), .Names = c("parent", "child", 
        "time"), row.names = integer(0), class = "data.frame"), 
        UnrecoverExcept = FALSE), .Names = c("attemptsUsed", 
    "errorMF", "errorMF_size",  "minDMratio", "minBMratio", 
     "PhylogDF", "UnrecoverExcept")), Drivers = integer(0), 
    geneNames = c("a", "b", "n1", "n2", "n3", "n4", "n5")), .Names = c("pops.by.time", 
"NumClones", "TotalPopSize", "Genotypes", "GenotypesWDistinctOrderEff", 
"GenotypesLabels", "MaxNumDrivers", "MaxDriversLast", "NumDriversLargestPop", 
"LargestClone", "PropLargestPopLast", "FinalTime", "NumIter", 
"HittedWallTime", "HittedMaxTries", "TotalPresentDrivers", "CountByDriver", 
"OccurringDrivers", "PerSampleStats", "other", "Drivers", "geneNames"
), class = c("oncosimul", "oncosimul2"), call = oncoSimulIndiv(fp = fp, 
    model = "Exp", numPassengers = numPassengers, mu = mu, detectionSize = 1e9, detectionDrivers = 99, 
    sampleEvery = sampleEvery, initSize = initSize, s = s, sh = sh, 
    K = K, keepEvery = -9, minDetectDrvCloneSz = "auto", extraTime = 0, 
    finalTime = finalTime, onlyCancer = onlyCancer, keepPhylog = keepPhylog, 
    mutationPropGrowth = mutationPropGrowth, max.memory = max.memory, 
    max.wall.time = max.wall.time.total, max.num.tries = possibleAttempts, 
    errorHitWallTime = TRUE, errorHitMaxTries = TRUE, verbosity = verbosity, 
    initMutant = initMutant, seed = seed))), class = "oncosimulpop")
cn <- c("a", "b", paste0("n", 1:5))
m1 <- matrix(c(rep(1, 5), rep(0, 30)),
             ncol = 7,
             dimnames = list(c(NULL), cn))
expect_identical( samplePop(s1, timeSample = "last", typeSample = "whole",
                            thresholdWhole = 0.5),
                 m1
                 )
m1 <- matrix(1, nrow = 5, 
             ncol = 7,
             dimnames = list(c(NULL), cn))
expect_identical( samplePop(s1, timeSample = "last", typeSample = "whole",
                            thresholdWhole = 0),
                 m1
                 )
## from s1[[5]], the two possible differentiating thresholds, right on the
## equality
m1 <- matrix(c(rep(1, 5), rep(0, 30)),
             ncol = 7,
             dimnames = list(c(NULL), cn))
m1[rbind(c(1, 2), c(1, 5), c(1, 6), c(1, 7), c(3, 5), c(5, 3))] <- 1
expect_identical(samplePop(s1, timeSample = "last", typeSample = "whole",
                           thresholdWhole = 1/64376),
                 m1)
m1 <- matrix(c(rep(1, 5), rep(0, 30)),
             ncol = 7,
             dimnames = list(c(NULL), cn))
m1[rbind(c(1, 2), c(1, 5), c(1, 6), c(1, 7), c(3, 5))] <- 1
expect_identical(samplePop(s1, timeSample = "last", typeSample = "whole",
                           thresholdWhole = 1/64375),
                 m1)
## from s1[[3]], the two possible differentiating thresholds, right on the
## equality; also affects first population
m1 <- matrix(c(rep(1, 5), rep(0, 30)),
             ncol = 7,
             dimnames = list(c(NULL), cn))
m1[rbind( c(1, 5),  c(1, 7), c(3, 5))] <- 1
expect_identical(samplePop(s1, timeSample = "last", typeSample = "whole",
                           thresholdWhole = 4/65814),
                 m1)
m1 <- matrix(c(rep(1, 5), rep(0, 30)),
             ncol = 7,
             dimnames = list(c(NULL), cn))
m1[rbind(c(1, 7))] <- 1
expect_identical(samplePop(s1, timeSample = "last", typeSample = "whole",
                           thresholdWhole = 5/65814),
                 m1)
## from s1[[1]], the two possible differentiating thresholds, right on the
## equality; also affects first population
m1 <- matrix(c(rep(1, 5), rep(0, 30)),
             ncol = 7,
             dimnames = list(c(NULL), cn))
m1[rbind(c(1, 2), c(1, 5), c(1, 6), c(1, 7), c(3, 5))] <- 1
expect_identical(samplePop(s1, timeSample = "last", typeSample = "whole",
                           thresholdWhole = 1/62989),
                 m1)
m1 <- matrix(c(rep(1, 5), rep(0, 30)),
             ncol = 7,
             dimnames = list(c(NULL), cn))
m1[rbind(c(1, 2), c(1, 5), c(1, 7), c(3, 5))] <- 1
expect_identical(samplePop(s1, timeSample = "last", typeSample = "whole",
                           thresholdWhole = 2/62989),
                 m1)
m1 <- matrix(c(rep(1, 5), rep(0, 30)),
             ncol = 7,
             dimnames = list(c(NULL), cn))
m1[rbind(c(1, 5),  c(1, 7))] <- 1
expect_identical(samplePop(s1, timeSample = "last", typeSample = "whole",
                           thresholdWhole = 4/62989),
                 m1)
m1 <- matrix(c(rep(1, 5), rep(0, 30)),
             ncol = 7,
             dimnames = list(c(NULL), cn))
m1[rbind(c(1, 7))] <- 1
expect_identical(samplePop(s1, timeSample = "last", typeSample = "whole",
                           thresholdWhole = 57/62989),
                 m1)
m1 <- matrix(c(rep(1, 5), rep(0, 30)),
             ncol = 7,
             dimnames = list(c(NULL), cn))
expect_identical(samplePop(s1, timeSample = "last", typeSample = "whole",
                           thresholdWhole = 58/62989),
                 m1)
})









test_that("exercising oncoSimulIndiv, new format, extra time", {
    pancr <- allFitnessEffects(data.frame(parent = c("Root", rep("KRAS", 4), "SMAD4", "CDNK2A", 
                                                     "TP53", "TP53", "MLL3"),
                                          child = c("KRAS","SMAD4", "CDNK2A", 
                                                    "TP53", "MLL3",
                                                    rep("PXDN", 3), rep("TGFBR2", 2)),
                                          s = 0.05,
                                          sh = -0.3,
                                          typeDep = "MN"))
    expect_silent(
        pSs <- oncoSimulIndiv(pancr,
                              sampleEvery = 0.03,
                              detectionSize = 1e3,
                              finalTime = 6,
                              extraTime = 5,
                              onlyCancer = FALSE))
    expect_output(print(pSs), "Individual OncoSimul trajectory",
                  fixed = TRUE)
})



test_that("exercising oncoSimulIndiv, hit max ram", {
    p1 <- allFitnessEffects(noIntGenes = rep(.1, 10))
    expect_output(pSs <- oncoSimulIndiv(p1,
                                        initSize = 1e4,
                                        sampleEvery = 3,
                                        detectionSize = 1e5,
                                        finalTime = 1000,
                                        extraTime = 5,
                                        onlyCancer = FALSE,
                                        max.memory = .01),
                  "Return outNS object > maxram",
                  fixed = TRUE)
})



test_that("exercising oncoSimulIndiv, verbosity", {
    ii <- rep(.1, 20)
    names(ii) <- letters[1:20]
    p1 <- allFitnessEffects(noIntGenes = ii)
    expect_output(pSs <- oncoSimulIndiv(p1,
                                        initMutant = "a",
                                        model = "McFL",
                                        initSize = 1e2,
                                        sampleEvery = 1,
                                        detectionSize = 1e10,
                                        finalTime = 2000,
                                        extraTime = 5,
                                        verbosity = 6,
                                        onlyCancer = FALSE),
                  "Looping through", fixed = TRUE)
    ## This is too much: can take a minute.
    ## ii <- rep(.01, 1000)
    ## names(ii) <- paste0("n", 1:1000)
    ## p1 <- allFitnessEffects(noIntGenes = ii)
    ## expect_output(pSs <- oncoSimulIndiv(p1,
    ##                                     initMutant = "n1",
    ##                                     model = "Exp",
    ##                                     initSize = 1e6,
    ##                                     sampleEvery = 2,
    ##                                     detectionSize = 1e7,
    ##                                     finalTime = 50,
    ##                                     extraTime = 5,
    ##                                     verbosity = 2,
    ##                                     onlyCancer = FALSE),
    ##               "Looping through", fixed = TRUE)
})



test_that("oncoSimulIndiv miscell C++ warnings", {
    pancri <- allFitnessEffects(
        data.frame(parent = c("Root", rep("KRAS", 4), "SMAD4", "CDNK2A", 
                              "TP53", "TP53", "MLL3"),
                   child = c("KRAS","SMAD4", "CDNK2A", 
                             "TP53", "MLL3",
                             rep("PXDN", 3), rep("TGFBR2", 2)),
                   s = Inf,
                   sh = Inf,
                   typeDep = "MN"))
    ## This are messages in R, not R warnings, as from Rcpp::Rcout
    expect_output(oncoSimulIndiv(pancri,
                                  initSize = 1,
                                  keepEvery = 5,
                                  onlyCancer = FALSE),
                   "at least one sh is positive infinite",
                   fixed = TRUE)
    expect_output(oncoSimulIndiv(pancri,
                                  initSize = 1,
                                  keepEvery = 5,
                                  onlyCancer = FALSE),
                   "at least one s is infinite",
                   fixed = TRUE)
})


test_that("using old poset format: Bozic and sh < 0", {
              data(examplePosets)
              p701 <- examplePosets[["p701"]]
              expect_silent(oncoSimulIndiv(p701, sh = -0.3,
                                           model = "Bozic",
                                           onlyCancer = FALSE))
})


test_that("using old poset format, extra time", {
              data(examplePosets)
              p701 <- examplePosets[["p701"]]
              expect_silent(oncoSimulIndiv(p701, sh = 0.3,
                                           initSize = 1000,
                                           detectionSize = 1100,
                                           model = "McFL",
                                           finalTime = 1000,
                                           extraTime = 3.17,
                                           onlyCancer = FALSE))
              expect_silent(oncoSimulIndiv(p701, sh = 0,
                                           model = "Exp",
                                           initSize = 1e4,
                                           detectionSize = 1e6,
                                           extraTime = 10,
                                           onlyCancer = FALSE))
              expect_silent(oncoSimulIndiv(p701, sh = -0.01,
                                           model = "Bozic",
                                           initSize = 1e4,
                                           detectionSize = 1e6,
                                           extraTime = 10,
                                           onlyCancer = FALSE))
})

test_that("using old poset format, no initMutant", {
    data(examplePosets)
    p701 <- examplePosets[["p701"]]
    expect_warning(p1 <- oncoSimulIndiv(p701, sh = 0.3,
                                        initSize = 1000,
                                        detectionSize = 1100,
                                        model = "Exp",
                                        finalTime = 1000,
                                        extraTime = 3.17,
                                        initMutant = 2,
                                        onlyCancer = FALSE),
                   "With the old poset format you can no longer use initMutant",
                   fixed = TRUE)
    expect_warning(p1 <- oncoSimulIndiv(p701, sh = 0.3,
                                        initSize = 1000,
                                        detectionSize = 1100,
                                        model = "Exp",
                                        finalTime = 1000,
                                        extraTime = 3.17,
                                        initMutant = c(2, 5),
                                        onlyCancer = FALSE),
                   "With the old poset format you can no longer use initMutant",
                   fixed = TRUE)
})

test_that("using old poset format, exercise verbosity", {
    data(examplePosets)
    p701 <- examplePosets[["p701"]]
    expect_output(oncoSimulIndiv(p701, sh = 0,
                                 initSize = 10000,
                                 detectionSize = 30000,
                                 model = "Exp",
                                 finalTime = 2000,
                                 extraTime = 3.17,
                                 verbosity = 10,
                                 onlyCancer = FALSE),
                  "Looping", fixed = TRUE)
})



test_that("exercising verbosity, new format", {
    gg <- rep(0.1, 20)
    names(gg) <- letters[1:20]
    ii <- allFitnessEffects(noIntGenes = gg)
    expect_output(
        pSs <- oncoSimulIndiv(ii,
                              sampleEvery = .5,
                              initSize = 1e4,
                              detectionSize = 5e5,
                              finalTime = 1000,
                              extraTime = 5,
                              keepEvery = 50,
                              verbosity = 10,
                              onlyCancer = FALSE),
        "Looping",
        fixed = TRUE)
})



test_that("old format: at most 64 genes", {
    p1 <- cbind(1L, 2L)
    expect_error(p1 <- oncoSimulIndiv(p1,
                        numPassengers = 66,
                        sh = 0,
                        initSize = 1e5,
                        sampleEvery = 0.02,
                        detectionSize = 1e9,
                        model = "Exp",
                        finalTime = 2000,
                        extraTime = 3.17,
                        onlyCancer = FALSE,
                        seed = NULL),
                 "Largest possible number of genes",
                 fixed = TRUE)
})



test_that("samplePop deals with failures in simuls", {
    
    fe <- allFitnessEffects(noIntGenes = c(-0.1, -0.2, -0.3))
    uu <- oncoSimulIndiv(fe, max.wall.time = 0.2, max.num.tries = 5)
    uup <- oncoSimulPop(4, fe, max.wall.time = 0.2, max.num.tries = 5,
                        mc.cores = 2)
    expect_warning(uus <- samplePop(uu),
                   "It looks like this simulation never completed",
                   fixed = TRUE)
    expect_warning(uups <- samplePop(uup),
                   "It looks like this simulation never completed",
                   fixed = TRUE)
    ## And it works when only some fail
    
    fe2 <- allFitnessEffects(noIntGenes = c(0.1, 0.2, 0.3))
    uu2 <- oncoSimulIndiv(fe2)
    u3 <- list(uu, uu2)
    class(u3) <- "oncosimulpop"
    expect_warning(uu3ps <- samplePop(u3),
                   "It looks like this simulation never completed",
                   fixed = TRUE)
})



test_that("summary.oncosimulepop deals with failures in simuls", {
    
    fe <- allFitnessEffects(noIntGenes = c(-0.1, -0.2, -0.3))
    uup <- oncoSimulPop(4, fe, max.wall.time = 0.2, max.num.tries = 5, mc.cores = 2)
    expect_warning(uus <- summary(uup),
                   "All simulations failed",
                   fixed = TRUE)
    ## And it works when only some fail
    fe2 <- allFitnessEffects(noIntGenes = c(0.1, 0.2, 0.3))
    uu2 <- oncoSimulPop(2, fe2, mc.cores = 2)
    uu <- oncoSimulPop(2, fe, max.wall.time = 0.2, max.num.tries = 5, mc.cores = 2)
    u3 <- c(uu, uu2)
    class(u3) <- "oncosimulpop"
    expect_warning(uu3ps <- summary(u3),
                   "Some simulations seem to have failed",
                   fixed = TRUE)
})


test_that("AND_DrvProbExit warnings and errors work" , {
    fe <- allFitnessEffects(noIntGenes = c("A" = 0.1, "B" = 0.2, "C" = 0.3),
                        drvNames = c("A", "B", "C"))

    expect_warning(u <- oncoSimulIndiv(fe, detectionDrivers = 1,
                        detectionProb = "default",
                        AND_DrvProbExit = TRUE),
                   "With AND_DrvProbExit = TRUE, detectionSize is ignored",
                   fixed = TRUE)

    expect_error(u <- oncoSimulIndiv(fe, detectionDrivers = NA,
                                     detectionProb = "default",
                                     AND_DrvProbExit = TRUE),
                 "AND_DrvProbExit is TRUE: both of detectionProb",
                 fixed = TRUE)


    data(examplePosets)
    p701 <- examplePosets[["p701"]]
    expect_error(u <- oncoSimulIndiv(p701, detectionDrivers = 1,
                                     detectionProb = "default",
                                     detectionSize = NA,
                                     AND_DrvProbExit = TRUE),
                 "The AND_DrvProbExit = TRUE setting is invalid with the old poset",
                 fixed = TRUE)
}
)


test_that("AND_DrvProbExit exercising and test it works" , {
    ## RNGkind("Mersenne-Twister")
    ## set.seed(1)

    fe <- allFitnessEffects(noIntGenes = c("A" = 0.1, "B" = 0.2, "C" = 0.3),
                        drvNames = c("A", "B"))
    ## set.seed(1)
    u0 <- oncoSimulIndiv(fe, detectionDrivers = 1,
                        detectionProb = "default",
                        AND_DrvProbExit = FALSE,
                        detectionSize = NA)
    ## set.seed(1)
    u <- oncoSimulIndiv(fe, detectionDrivers = 1,
                        detectionProb = "default",
                        AND_DrvProbExit = TRUE,
                        detectionSize = NA)
    expect_true(u$TotalPresentDrivers >= 1)
    ## set.seed(1)
    u2 <- oncoSimulIndiv(fe, detectionDrivers = 1,
                        detectionProb = "default",
                        AND_DrvProbExit = TRUE,
                        minDetectDrvCloneSz = 500,
                        detectionSize = NA)
    expect_true(u2$TotalPresentDrivers >= 1)
    ## set.seed(2)

    ## set.seed(2)
    ## m00 <- oncoSimulIndiv(fe, detectionDrivers = NA, model = "McFL",
    ##                     detectionProb = "default",
    ##                     AND_DrvProbExit = FALSE,
    ##                     detectionSize = NA,
    ##                     initMutant = "C")
    ## set.seed(2)
    m <- oncoSimulIndiv(fe, detectionDrivers = 1, model = "McFL",
                        detectionProb = "default",
                        AND_DrvProbExit = TRUE,
                        detectionSize = NA)
    expect_true(m$TotalPresentDrivers >= 1)
    ## set.seed(2)
    m1 <- oncoSimulIndiv(fe, detectionDrivers = 1, model = "McFL",
                        detectionProb = "default",
                        AND_DrvProbExit = TRUE,
                        minDetectDrvCloneSz = 10,
                        detectionSize = NA)
    expect_true(m1$TotalPresentDrivers >= 1)
    m1p <- m1$pops.by.time[, -c(1, 2), drop = FALSE]
    expect_true(sum(m1p[nrow(m1p), , drop = FALSE]) >= 10)
    ## set.seed(2)
    m2 <- oncoSimulIndiv(fe, detectionDrivers = 1, model = "McFL",
                        detectionProb = "default",
                        AND_DrvProbExit = TRUE,
                        minDetectDrvCloneSz = 100,
                        detectionSize = NA)
    expect_true(m2$TotalPresentDrivers >= 1)
    m2p <- m2$pops.by.time[, -c(1, 2), drop = FALSE]
    expect_true(sum(m2p[nrow(m2p), , drop = FALSE]) >= 100)
    ## set.seed(2)
    ## slow
    ## m3 <- oncoSimulIndiv(fe, detectionDrivers = 1, model = "McFL",
    ##                     detectionProb = "default",
    ##                     AND_DrvProbExit = TRUE,
    ##                     minDetectDrvCloneSz = 1200,
    ##                     detectionSize = NA)
    m3 <- oncoSimulIndiv(fe, detectionDrivers = 1, model = "McFL",
                         detectionProb = "default",
                         extraTime = 3,
                         AND_DrvProbExit = TRUE,
                         minDetectDrvCloneSz = 100,
                         detectionSize = NA)
    expect_true(m3$TotalPresentDrivers >= 1)
    m3p <- m3$pops.by.time[, -c(1, 2), drop = FALSE]
    expect_true(sum(m3p[nrow(m3p), , drop = FALSE]) >= 100)
    m3 <- oncoSimulIndiv(fe, detectionDrivers = 1, model = "McFL",
                         detectionProb = "default",
                         extraTime = 10,
                         AND_DrvProbExit = TRUE,
                         minDetectDrvCloneSz = 100,
                         detectionSize = NA)
    expect_true(m3$TotalPresentDrivers >= 1)
    m3p <- m3$pops.by.time[, -c(1, 2), drop = FALSE]
    expect_true(sum(m3p[nrow(m3p), , drop = FALSE]) >= 100)
}
)




test_that("exercising oncoSimulIndiv, max size warning", {
    p1 <- allFitnessEffects(noIntGenes = rep(.1, 10))
    expect_output(oncoSimulIndiv(p1, initSize = 1.5e15, verbosity = 1,
                                 onlyCancer = FALSE, mu= 1e-7))
    data(examplePosets)
    p701 <- examplePosets[["p701"]]
    expect_output(oncoSimulIndiv(p701, initSize = 4.1e15, verbosity = 1,
                                 onlyCancer = FALSE, mu= 1e-7))
})




cat(paste("\n Ending oncoSimulIndiv-miscell tests", date(), "\n"))







cat(paste("  Took ", round(difftime(Sys.time(), inittime, units = "secs"), 2), "\n\n"))
rm(inittime)
