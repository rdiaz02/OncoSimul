inittime <- Sys.time()
cat(paste("\n Starting samplePop tests", date(), "\n"))

## RNGkind("Mersenne-Twister")

test_that("exercising the sampling code, v1 objects", {
              data(examplePosets)
              p705 <- examplePosets[["p705"]]
              r1 <- oncoSimulIndiv(p705,
                                   sampleEvery = 0.03, keepEvery = 1
                                   )
              p1 <- oncoSimulPop(4, p705, numPassengers = 30, mc.cores = 2,
                                 sampleEvery = 0.03, keepEvery = 1)
              expect_message(samplePop(p1  ),
                            "Subjects by Genes matrix of 4 subjects and 37 genes")
              expect_message(samplePop(p1, typeSample = "single"),
                             "Subjects by Genes matrix of 4 subjects and 37 genes")
              expect_message(samplePop(p1, typeSample = "whole", timeSample = "unif"),
                             "Subjects by Genes matrix of 4 subjects and 37 genes")
              expect_message(samplePop(p1, typeSample = "single", timeSample = "unif"),
                            "Subjects by Genes matrix of 4 subjects and 37 genes")
              expect_message(samplePop(r1, typeSample = "whole", timeSample = "last"),
                             "Subjects by Genes matrix of 1 subjects and 7 genes")
              expect_message(samplePop(r1, typeSample = "single", timeSample = "last"),
                            "Subjects by Genes matrix of 1 subjects and 7 genes")
              expect_message(samplePop(r1, typeSample = "whole", timeSample = "unif"),
                            "Subjects by Genes matrix of 1 subjects and 7 genes")
              expect_message(samplePop(r1, typeSample = "single", timeSample = "uniform"),
                            "Subjects by Genes matrix of 1 subjects and 7 genes")
          })

test_that("exercising the sampling code, v1 objects, possibly NAs", {
              data(examplePosets)
              p705 <- examplePosets[["p705"]]
              r1 <- oncoSimulIndiv(p705, numPassengers = 30,
                                   sampleEvery = 0.03, keepEvery = 1,
                                   onlyCancer = FALSE)
              p1 <- oncoSimulPop(4, p705, mc.cores = 2, onlyCancer = FALSE,
                                 sampleEvery = 0.03, keepEvery = 1)
              expect_message(samplePop(p1),
                            "Subjects by Genes matrix of 4 subjects and 7 genes")
              expect_message(samplePop(p1, typeSample = "whole",
                                       timeSample = "last"),
                             "Subjects by Genes matrix of 4 subjects and 7 genes")
              expect_message(samplePop(p1, typeSample = "whole",
                                       timeSample = "unif"),
                             "Subjects by Genes matrix of 4 subjects and 7 genes")
              expect_message(samplePop(p1, typeSample = "single",
                                       timeSample = "unif"),
                             "Subjects by Genes matrix of 4 subjects and 7 genes")
              expect_message(samplePop(p1, typeSample = "single",
                                       timeSample = "last"),
                             "Subjects by Genes matrix of 4 subjects and 7 genes")
              expect_message(samplePop(p1, typeSample = "single"),
                             "Subjects by Genes matrix of 4 subjects and 7 genes")
              expect_message(samplePop(r1),
                             "Subjects by Genes matrix of 1 subjects and 37 genes")
              expect_message(samplePop(r1, typeSample = "whole",
                                       timeSample = "last"),
                             "Subjects by Genes matrix of 1 subjects and 37 genes")
              expect_message(samplePop(r1, typeSample = "whole",
                                       timeSample = "unif"),
                            "Subjects by Genes matrix of 1 subjects and 37 genes")
              expect_message(samplePop(r1, typeSample = "single",
                                       timeSample = "last"),
                            "Subjects by Genes matrix of 1 subjects and 37 genes")
              expect_message(samplePop(r1, typeSample = "single",
                                       timeSample = "uniform"),
                            "Subjects by Genes matrix of 1 subjects and 37 genes")
          })

test_that("exercising the sampling code, v1 objects, possibly NAs, more", {
              data(examplePosets)
              p1101 <- examplePosets[["p1101"]]
              r1 <- oncoSimulIndiv(p1101, onlyCancer = FALSE,
                                   sampleEvery = 0.03, keepEvery = 1)
              p1 <- oncoSimulPop(4, p1101, mc.cores = 2, onlyCancer = FALSE,
                                 sampleEvery = 0.03, keepEvery = 1)
              expect_message(samplePop(p1),
                            "Subjects by Genes matrix of 4 subjects and 11 genes")
              expect_message(samplePop(p1, typeSample = "whole",
                                       timeSample = "last"),
                             "Subjects by Genes matrix of 4 subjects and 11 genes")
              expect_message(samplePop(p1, typeSample = "whole",
                                       timeSample = "unif"),
                             "Subjects by Genes matrix of 4 subjects and 11 genes")
              expect_message(samplePop(p1, typeSample = "single",
                                       timeSample = "unif"),
                             "Subjects by Genes matrix of 4 subjects and 11 genes")
              expect_message(samplePop(p1, typeSample = "single",
                                       timeSample = "last"),
                             "Subjects by Genes matrix of 4 subjects and 11 genes")
              expect_message(samplePop(p1, typeSample = "single"),
                             "Subjects by Genes matrix of 4 subjects and 11 genes")
              expect_message(samplePop(r1),
                             "Subjects by Genes matrix of 1 subjects and 11 genes")
              expect_message(samplePop(r1, typeSample = "whole",
                                       timeSample = "last"),
                             "Subjects by Genes matrix of 1 subjects and 11 genes")
              expect_message(samplePop(r1, typeSample = "whole",
                                       timeSample = "unif"),
                            "Subjects by Genes matrix of 1 subjects and 11 genes")
              expect_message(samplePop(r1, typeSample = "single",
                                       timeSample = "last"),
                            "Subjects by Genes matrix of 1 subjects and 11 genes")
              expect_message(samplePop(r1, typeSample = "single",
                                       timeSample = "uniform"),
                            "Subjects by Genes matrix of 1 subjects and 11 genes")
          })

test_that("exercising the sampling code, v2 objects, possibly NAs in output", {
              oi <- allFitnessEffects(orderEffects =
                                          c("F > D" = -0.3, "D > F" = 0.4),
                                      noIntGenes = rexp(5, 10),
                                      geneToModule =
                                          c("Root" = "Root",
                                            "F" = "f1, f2, f3",
                                            "D" = "d1, d2") )
              o1 <- oncoSimulIndiv(oi, detectionSize = 1e4,
                                   onlyCancer = FALSE,
                                   sampleEvery = 0.03, keepEvery = 1)
              o4 <- oncoSimulPop(2,
                                 oi, 
                                 detectionSize = 1e4,
                                 mc.cores = 2,
                                 onlyCancer = FALSE,
                                 sampleEvery = 0.03, keepEvery = 1)
              ## many of them are generating warnings, because sampling
              ## with pop size of 0. That is OK.
              expect_message(samplePop(o1),
                             "Subjects by Genes matrix of 1 subjects and 10 genes")
              expect_message(samplePop(o1, typeSample = "single",
                                       timeSample = "last"),
                            "Subjects by Genes matrix of 1 subjects and 10 genes")
              expect_message(samplePop(o1, typeSample = "single",
                                       timeSample = "unif"),
                            "Subjects by Genes matrix of 1 subjects and 10 genes")
              expect_message(samplePop(o1, typeSample = "whole",
                                       timeSample = "unif"),
                             "Subjects by Genes matrix of 1 subjects and 10 genes")
              expect_message(samplePop(o1, typeSample = "whole",
                                       timeSample = "last"),
                             "Subjects by Genes matrix of 1 subjects and 10 genes")
              expect_message(samplePop(o4),
                                           "Subjects by Genes matrix of 2 subjects and 10 genes")
              expect_message(samplePop(o4, typeSample = "single",
                                       timeSample = "last"),
                            "Subjects by Genes matrix of 2 subjects and 10 genes")
              expect_message(samplePop(o4, typeSample = "single",
                                       timeSample = "uniform"),
                             "Subjects by Genes matrix of 2 subjects and 10 genes")
              expect_message(samplePop(o4, timeSample = "unif",
                                       typeSample = "whole"),
                             "Subjects by Genes matrix of 2 subjects and 10 genes")
              expect_message(samplePop(o4, timeSample = "last",
                                       typeSample = "whole"),
                             "Subjects by Genes matrix of 2 subjects and 10 genes")
          })

test_that("exercising the sampling code, v2 objects, possibly NAs in output, more", {
              cs <-  data.frame(parent = c(rep("Root", 4), "a", "b", "d", "e", "c"),
                                child = c("a", "b", "d", "e", "c", "c", rep("g", 3)),
                                s = 0.1,
                                sh = -0.9,
                                typeDep = "MN")
              cbn1 <- allFitnessEffects(cs)
              o1 <- oncoSimulIndiv(cbn1, detectionSize = 1e4,
                                   onlyCancer = FALSE,
                                   sampleEvery = 0.03, keepEvery = 1)
              o4 <- oncoSimulPop(4,
                                 cbn1, 
                                 detectionSize = 1e4,
                                 onlyCancer = FALSE,
                                 mc.cores = 2,
                                 sampleEvery = 0.03, keepEvery = 1)
              ## many of them are generating warnings, because sampling
              ## with pop size of 0. That is OK.
              expect_message(samplePop(o1),
                             "Subjects by Genes matrix of 1 subjects and 6 genes")
              expect_message(samplePop(o1, typeSample = "single",
                                       timeSample = "last"),
                            "Subjects by Genes matrix of 1 subjects and 6 genes")
              expect_message(samplePop(o1, typeSample = "single",
                                       timeSample = "unif"),
                            "Subjects by Genes matrix of 1 subjects and 6 genes")
              expect_message(samplePop(o1, typeSample = "whole",
                                       timeSample = "unif"),
                             "Subjects by Genes matrix of 1 subjects and 6 genes")
              expect_message(samplePop(o1, typeSample = "whole",
                                       timeSample = "last"),
                             "Subjects by Genes matrix of 1 subjects and 6 genes")
              expect_message(samplePop(o4),
                                           "Subjects by Genes matrix of 4 subjects and 6 genes")
              expect_message(samplePop(o4, typeSample = "single",
                                       timeSample = "last"),
                            "Subjects by Genes matrix of 4 subjects and 6 genes")
              expect_message(samplePop(o4, typeSample = "single",
                                       timeSample = "uniform"),
                             "Subjects by Genes matrix of 4 subjects and 6 genes")
              expect_message(samplePop(o4, timeSample = "unif",
                                       typeSample = "whole"),
                             "Subjects by Genes matrix of 4 subjects and 6 genes")
              expect_message(samplePop(o4, timeSample = "last",
                                       typeSample = "whole"),
                             "Subjects by Genes matrix of 4 subjects and 6 genes")
          })


test_that("exercising the sampling code, v2 objects", {
              oi <- allFitnessEffects(orderEffects =
                                          c("F > D" = -0.3, "D > F" = 0.4),
                                      noIntGenes = rexp(5, 10),
                                      geneToModule =
                                          c("Root" = "Root",
                                            "F" = "f1, f2, f3",
                                            "D" = "d1, d2") )
              o1 <- oncoSimulIndiv(oi, detectionSize = 1e4,
                                   onlyCancer = TRUE,
                                   max.num.tries = 5000,
                                   sampleEvery = 0.03, keepEvery = 1)
              o4 <- oncoSimulPop(2,
                                 oi, 
                                 detectionSize = 1e4,
                                 onlyCancer = TRUE,
                                 mc.cores = 2,
                                 max.num.tries = 5000,
                                 sampleEvery = 0.03, keepEvery = 1)
              expect_message(samplePop(o1),
                             "Subjects by Genes matrix of 1 subjects and 10 genes")
              expect_message(samplePop(o1, typeSample = "single",
                                       timeSample = "last"),
                            "Subjects by Genes matrix of 1 subjects and 10 genes")
              expect_message(samplePop(o1, typeSample = "single",
                                       timeSample = "unif"),
                            "Subjects by Genes matrix of 1 subjects and 10 genes")
              expect_message(samplePop(o1, typeSample = "whole",
                                       timeSample = "unif"),
                             "Subjects by Genes matrix of 1 subjects and 10 genes")
              expect_message(samplePop(o1, typeSample = "whole",
                                       timeSample = "last"),
                             "Subjects by Genes matrix of 1 subjects and 10 genes")
              expect_message(samplePop(o4),
                                           "Subjects by Genes matrix of 2 subjects and 10 genes")
              expect_message(samplePop(o4, typeSample = "single",
                                       timeSample = "last"),
                            "Subjects by Genes matrix of 2 subjects and 10 genes")
              expect_message(samplePop(o4, typeSample = "single",
                                       timeSample = "uniform"),
                             "Subjects by Genes matrix of 2 subjects and 10 genes")
              expect_message(samplePop(o4, timeSample = "unif",
                                       typeSample = "whole"),
                             "Subjects by Genes matrix of 2 subjects and 10 genes")
              expect_message(samplePop(o4, timeSample = "last",
                                       typeSample = "whole"),
                             "Subjects by Genes matrix of 2 subjects and 10 genes")
              expect_message(samplePop(o4, typeSample = "single-nowt",
                                       timeSample = "last"),
                            "Subjects by Genes matrix of 2 subjects and 10 genes")
              expect_message(samplePop(o4, typeSample = "single-nowt",
                                       timeSample = "uniform"),
                             "Subjects by Genes matrix of 2 subjects and 10 genes")
})

test_that("exercising the sampling code, v2 objects, more", {
              cs <-  data.frame(parent = c(rep("Root", 4), "a", "b", "d", "e", "c"),
                                child = c("a", "b", "d", "e", "c", "c", rep("g", 3)),
                                s = 0.1,
                                sh = -0.9,
                                typeDep = "MN")
              cbn1 <- allFitnessEffects(cs)
              o1 <- oncoSimulIndiv(cbn1, detectionSize = 1e4,
                                   onlyCancer = TRUE,
                                   max.num.tries = 5000,
                                   sampleEvery = 0.03, keepEvery = 1)
              o4 <- oncoSimulPop(4,
                                 cbn1, 
                                 detectionSize = 1e4,
                                 onlyCancer = TRUE,
                                 mc.cores = 2,
                                 max.num.tries = 5000,
                                 sampleEvery = 0.03, keepEvery = 1)
              expect_message(samplePop(o1),
                             "Subjects by Genes matrix of 1 subjects and 6 genes")
              expect_message(samplePop(o1, typeSample = "single",
                                       timeSample = "last"),
                            "Subjects by Genes matrix of 1 subjects and 6 genes")
              expect_message(samplePop(o1, typeSample = "single",
                                       timeSample = "unif"),
                            "Subjects by Genes matrix of 1 subjects and 6 genes")
              expect_message(samplePop(o1, typeSample = "whole",
                                       timeSample = "unif"),
                             "Subjects by Genes matrix of 1 subjects and 6 genes")
              expect_message(samplePop(o1, typeSample = "whole",
                                       timeSample = "last"),
                             "Subjects by Genes matrix of 1 subjects and 6 genes")
              expect_message(samplePop(o4),
                                           "Subjects by Genes matrix of 4 subjects and 6 genes")
              expect_message(samplePop(o4, typeSample = "single",
                                       timeSample = "last"),
                            "Subjects by Genes matrix of 4 subjects and 6 genes")
              expect_message(samplePop(o4, typeSample = "single",
                                       timeSample = "uniform"),
                             "Subjects by Genes matrix of 4 subjects and 6 genes")
              expect_message(samplePop(o4, timeSample = "unif",
                                       typeSample = "whole"),
                             "Subjects by Genes matrix of 4 subjects and 6 genes")
              expect_message(samplePop(o4, timeSample = "last",
                                       typeSample = "whole"),
                             "Subjects by Genes matrix of 4 subjects and 6 genes")
              expect_message(samplePop(o4, typeSample = "single-nowt",
                                       timeSample = "last"),
                            "Subjects by Genes matrix of 4 subjects and 6 genes")
              expect_message(samplePop(o4, typeSample = "single-nowt",
                                       timeSample = "uniform"),
                             "Subjects by Genes matrix of 4 subjects and 6 genes")
})



test_that("exercising sampling code, single sampled period", {
    o3init <- allFitnessEffects(orderEffects = c(
                            "M > D > F" = 0.99,
                            "D > M > F" = 0.2,
                            "D > M"     = 0.1,
                            "M > D"     = 0.9),
                        noIntGenes = c("u" = 0.01, 
                                       "v" = 0.01,
                                       "w" = 0.001,
                                       "x" = 0.0001,
                                       "y" = -0.0001,
                                       "z" = 0.001),
                        geneToModule =
                            c("Root" = "Root",
                              "M" = "m",
                              "F" = "f",
                              "D" = "d") )
    for(i in 1:10) {
        i1 <- oncoSimulIndiv(o3init, model = "Exp", initSize = 1e8,
                             mu = 1e-3, finalTime = .01,
                             detectionDrivers = 2,
                             onlyCancer = FALSE,
                             sampleEvery = 0.03, keepEvery = 1
                             )
        expect_message(spi <- samplePop(i1, timeSample = "unif"),
                       "Subjects", fixed = TRUE)
    }
})



test_that("exercising sampling code, customSize", {

    cs <-  data.frame(parent = c(rep("Root", 4), "a", "b", "d", "e", "c"),
                                child = c("a", "b", "d", "e", "c", "c", rep("g", 3)),
                      s = 0.1,
                                sh = -0.9,
                      typeDep = "MN")
    
    cbn1 <- allFitnessEffects(cs)

    
    o1 <- oncoSimulIndiv(cbn1, detectionSize = 1e4,
                         onlyCancer = TRUE,
                         max.num.tries = 5000,
                         sampleEvery = 0.03, keepEvery = 1)
    o4 <- oncoSimulPop(4,
                       cbn1, 
                       detectionSize = 1e4,
                       onlyCancer = TRUE,
                       mc.cores = 2,
                       max.num.tries = 5000,
                       sampleEvery = 0.03, keepEvery = 1)

    expect_message(samplePop(o1, popSizeSample = 2100),
                   "Subjects by Genes matrix of 1 subjects and 6 genes")

    expect_warning(samplePop(o1, popSizeSample = 21000),
                   "Pop size never", fixed = TRUE)

    expect_message(samplePop(o4, popSizeSample = c(6100)),
                   "Subjects by Genes matrix of 4 subjects and 6 genes")

    expect_message(samplePop(o4, popSizeSample = c(6100, 6000, 8999, 8030)),
                   "Subjects by Genes matrix of 4 subjects and 6 genes")

    expect_message(samplePop(o4, popSizeSample = c(6100, 6000)),
                   "Subjects by Genes matrix of 4 subjects and 6 genes")
    ## this was fixed to not give warnings
    ## expect_warning(samplePop(o4, popSizeSample = c(6100, 6000)),
    ##                "length popSizeSample != number of subjects")

    expect_warning(samplePop(o4, popSizeSample = 21000),
                   "Pop size never >= requested size", fixed = TRUE)
    
    expect_message(samplePop(o4, typeSample = "single",
                             popSizeSample = c(6100, 0, 5000, 9000)),
                   "Subjects by Genes matrix of 4 subjects and 6 genes")

    expect_message(samplePop(o4, typeSample = "single-nowt",
                             popSizeSample = c(9000, 9000, 8500, 9000)),
                   "Subjects by Genes matrix of 4 subjects and 6 genes")

    ## these are specific for the nowt code to exercise the dealing with
    ## borderline cases and the second to deal with a clearly non-borderline
    o41 <- oncoSimulPop(4,
                        cbn1,
                        initSize = 2e3,
                       detectionSize = 1e3,
                       onlyCancer = TRUE,
                       mc.cores = 2,
                       max.num.tries = 5000,
                       sampleEvery = 0.03, keepEvery = 1)
    expect_message(samplePop(o41, typeSample = "single-nowt",
                             popSizeSample = c(900, 800, 850, 900)),
                   "Subjects by Genes matrix of 4 subjects and 6 genes")
    expect_warning(samplePop(o41, typeSample = "single-nowt",
                             popSizeSample = c(900, 800, 850, 900)),
                   "No non-WT clone with required popSize or at required time")
    o91 <- oncoSimulPop(4,
                       cbn1, 
                       detectionSize = 1e5,
                       onlyCancer = TRUE,
                       mc.cores = 2,
                       max.num.tries = 5000,
                       sampleEvery = 0.03, keepEvery = 1)
     expect_message(samplePop(o91, typeSample = "single-nowt",
                             popSizeSample = c(9000, 9000, 8500, 9000)),
                   "Subjects by Genes matrix of 4 subjects and 6 genes")
})



test_that("exercising sampling code, customSize", {
    cs <-  data.frame(parent = c(rep("Root", 4), "a", "b", "d", "e", "c"),
                                child = c("a", "b", "d", "e", "c", "c", rep("g", 3)),
                      s = 0.1,
                                sh = -0.9,
                      typeDep = "MN")
    cbn1 <- allFitnessEffects(cs, drvNames = c("a", "b", "c", "d", "e", "g"))
    o1 <- oncoSimulIndiv(cbn1, detectionSize = 1e8,
                         onlyCancer = TRUE,
                         detectionDrivers = 1,
                         max.num.tries = 5000,
                         sampleEvery = 0.03,
                         keepEvery = 1000000)
    expect_message(samplePop(o1, timeSample = "unif"),
                   "Only one sampled time period with mutants", fixed = TRUE)
})






test_that("Exercise sampledGenotypes, option for selecting genes", {
    initS <- 200
    u <- 0.2; i <- -0.02; vi <- 0.6; ui <- uv <- -Inf
    od <- allFitnessEffects(
        epistasis = c("u" = u,  "u:i" = ui,
                      "u:v" = uv, "i" = i,
                      "v:-i" = -Inf, "v:i" = vi))
    ## drvNames = c("u", "v"))
    odm <- allMutatorEffects(noIntGenes = c("i" = 50))
    evalAllGenotypesFitAndMut(od, odm, addwt = TRUE)
    op <- oncoSimulPop(20, od, muEF = odm, model = "McFL",
                        mu = 1e-4, 
                        onlyCancer = TRUE, finalTime = 15000, detectionSize = NA, detectionProb = NA,
                        initSize = initS, 
                        keepEvery = NA,
                        fixation = c("u", "v"),
                        mc.cores = 2
                        )
    sp <- samplePop(op)
    expect_output(print(sampledGenotypes(sp, genes = c("u", "i"))))
    expect_output(print(summary(sampledGenotypes(sp, genes = c("u", "i")))))
    expect_output(print(sampledGenotypes(sp)))
    expect_output(print(summary(sampledGenotypes(sp))))
})





test_that("exercising sampling code, propError", {
    cs <-  data.frame(parent = c(rep("Root", 4), "a", "b", "d", "e", "c"),
                                child = c("a", "b", "d", "e", "c", "c", rep("g", 3)),
                      s = 0.1,
                                sh = -0.9,
                      typeDep = "MN")
    cbn1 <- allFitnessEffects(cs)
    o4 <- oncoSimulPop(4,
                       cbn1, 
                       detectionSize = 1e4,
                       onlyCancer = TRUE,
                       mc.cores = 2,
                       max.num.tries = 5000,
                       sampleEvery = 0.03, keepEvery = 1)
    expect_error(samplePop(o4, propError = 1.4),
                 "Proportion with error cannot be > 1", fixed = TRUE)
    expect_message(samplePop(o4, propError = 0),
                   "Subjects by Genes matrix of 4 subjects",
                   fixed = TRUE)
    expect_message(samplePop(o4, propError = 0.1),
                   "Subjects by Genes matrix of 4 subjects",
                   fixed = TRUE)
    expect_message(samplePop(o4, propError = 0.2),
                   "Subjects by Genes matrix of 4 subjects",
                   fixed = TRUE)
    expect_message(samplePop(o4, propError = 0.3),
                   "Subjects by Genes matrix of 4 subjects",
                   fixed = TRUE)
    expect_message(samplePop(o4, propError = 0.9),
                   "Subjects by Genes matrix of 4 subjects",
                   fixed = TRUE)
    expect_message(samplePop(o4, propError = 0.01),
                   "Subjects by Genes matrix of 4 subjects",
                   fixed = TRUE)
    expect_message(samplePop(o4, propError = 0.05),
                   "Subjects by Genes matrix of 4 subjects",
                   fixed = TRUE)
    expect_message(samplePop(o4, propError = 1),
                   "Subjects by Genes matrix of 4 subjects",
                   fixed = TRUE)
})




test_that("sampledGenotypes deals with NAs", {
    sp <- structure(c(1, NA, 1, 0, 0, NA, 0, 0, 0, NA, 0, 0, 0, NA, 0, 
                      0, 0, NA, 1, 0, 0, NA, 0, 0, 0, NA, 0, 0, 1, NA, 1, 0),
                    .Dim = c(4L, 
                             8L),
                    .Dimnames = list(NULL, c("ZZZ", "APC", "EVC2",
                                             "VVV", 
                                             "KRAS", "PIK3CA", "TCF7L2", "TP53")))
    sg <- sampledGenotypes(sp)
    expect_true(is.na(sg[4, 1]))
    expect_true(sg[1, 1] == "WT")
})

cat(paste("\n Ending samplePop tests", date(), "\n"))





cat(paste("  Took ", round(difftime(Sys.time(), inittime, units = "secs"), 2), "\n\n"))
rm(inittime)
